"""
CIViC client — V2 GraphQL API + nightly TSV bulk download fallback.

CIViC is the primary open weapon:
  - CC0 licensed, fully public, no token required
  - GraphQL API at https://civicdb.org/api/graphql
  - Nightly TSV releases at https://civicdb.org/downloads/nightly/
  - Evidence levels: A (validated), B (clinical), C (case study), D (preclinical), E (inferential)

Evidence level → Tier mapping:
  A → Tier_A
  B → Tier_B
  C, D → Tier_C
  E → Tier_D
"""
from __future__ import annotations

import io
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import httpx
import pandas as pd

from ..models import EvidenceTier, KBEvidence, ResponseType, SourceKB

logger = logging.getLogger(__name__)

CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"
CIVIC_NIGHTLY_EVIDENCE_URL = (
    "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
)
CIVIC_NIGHTLY_ASSERTIONS_URL = (
    "https://civicdb.org/downloads/nightly/nightly-AssertionSummaries.tsv"
)

# Cache path
_CACHE_DIR = Path(".cache/civic")
_CACHE_DIR.mkdir(parents=True, exist_ok=True)
_TSV_CACHE = _CACHE_DIR / "evidence_summaries.parquet"


# ── Level → Tier ──────────────────────────────────────────────────────────────

_CIVIC_LEVEL_MAP: Dict[str, EvidenceTier] = {
    "A": EvidenceTier.A,
    "B": EvidenceTier.B,
    "C": EvidenceTier.C,
    "D": EvidenceTier.C,   # Case study → still clinically observed
    "E": EvidenceTier.D,   # Inferential
}


def _map_tier(level: str) -> EvidenceTier:
    return _CIVIC_LEVEL_MAP.get(str(level).strip().upper(), EvidenceTier.UNKNOWN)


def _map_significance(sig: str) -> ResponseType:
    s = sig.lower()
    if "sensitiv" in s or "response" in s or "benefit" in s:
        return ResponseType.SENSITIVITY
    if "resist" in s:
        return ResponseType.RESISTANCE
    if "reduced" in s:
        return ResponseType.REDUCED_SENSITIVITY
    if "adverse" in s or "toxicit" in s:
        return ResponseType.ADVERSE
    return ResponseType.UNKNOWN


# ── Bulk TSV loader ───────────────────────────────────────────────────────────

def _load_civic_tsv(force_refresh: bool = False) -> pd.DataFrame:
    """Download and cache CIViC nightly evidence summaries TSV."""
    if _TSV_CACHE.exists() and not force_refresh:
        logger.info("CIViC: loading from cache %s", _TSV_CACHE)
        return pd.read_parquet(_TSV_CACHE)

    logger.info("CIViC: downloading nightly TSV …")
    try:
        resp = httpx.get(CIVIC_NIGHTLY_EVIDENCE_URL, follow_redirects=True, timeout=120)
        resp.raise_for_status()
        df = pd.read_csv(io.BytesIO(resp.content), sep="\t", low_memory=False)
        df.to_parquet(_TSV_CACHE)
        logger.info("CIViC: cached %d evidence rows", len(df))
        return df
    except Exception as exc:
        logger.warning("CIViC TSV download failed: %s", exc)
        return pd.DataFrame()


# ── GraphQL query ─────────────────────────────────────────────────────────────

_GENE_EVIDENCE_QUERY = """
query GeneEvidence($geneSymbol: String!) {
  genes(entrezSymbol: $geneSymbol) {
    nodes {
      id
      name
      variants {
        nodes {
          id
          name
          molecularProfiles {
            nodes {
              id
              name
              evidenceItems(status: ACCEPTED) {
                nodes {
                  id
                  evidenceLevel
                  evidenceDirection
                  significanceLabel: significance
                  therapies {
                    id
                    name
                  }
                  disease {
                    name
                    doid
                  }
                  statement
                  source {
                    citation
                    pmid
                    clinicalTrials {
                      nctId
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
"""


def _graphql_query(gene: str) -> Optional[Dict[str, Any]]:
    """Execute CIViC GraphQL query for a gene. Returns raw data dict."""
    try:
        resp = httpx.post(
            CIVIC_GRAPHQL_URL,
            json={"query": _GENE_EVIDENCE_QUERY, "variables": {"geneSymbol": gene}},
            timeout=30,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        return resp.json().get("data")
    except Exception as exc:
        logger.warning("CIViC GraphQL failed for %s: %s", gene, exc)
        return None


# ── TSV-based fallback ────────────────────────────────────────────────────────

def _tsv_query(
    gene: str,
    variant: Optional[str],
    cancer_type: Optional[str],
    df: pd.DataFrame,
) -> List[KBEvidence]:
    """Query the cached TSV for evidence matching gene/variant/cancer_type."""
    if df.empty:
        return []

    # Column name normalisation
    gene_col = next((c for c in df.columns if c.lower() in ("gene", "gene_id", "entrez_name")), None)
    if gene_col is None:
        return []

    mask = df[gene_col].str.upper() == gene.upper()

    if variant:
        var_col = next((c for c in df.columns if c.lower() in ("variant", "variant_name", "variant_origin")), None)
        if var_col:
            mask &= df[var_col].str.contains(variant, case=False, na=False)

    if cancer_type:
        ct_col = next((c for c in df.columns if c.lower() in ("disease", "disease_name", "doid_name")), None)
        if ct_col:
            mask &= df[ct_col].str.contains(cancer_type, case=False, na=False)

    rows = df[mask]
    evidence: List[KBEvidence] = []

    therapy_col = next((c for c in df.columns if c.lower() in ("therapies", "therapy_name", "drugs")), None)
    level_col = next((c for c in df.columns if c.lower() in ("evidence_level", "evidencelevel")), None)
    sig_col = next((c for c in df.columns if c.lower() in ("significance", "significancelabel", "clinical_significance")), None)
    stmt_col = next((c for c in df.columns if c.lower() in ("statement", "description")), None)
    pmid_col = next((c for c in df.columns if c.lower() in ("pmid", "source_pmid", "source_citation_id")), None)
    var_col2 = next((c for c in df.columns if c.lower() in ("variant", "variant_name")), None)
    ct_col2 = next((c for c in df.columns if c.lower() in ("disease", "disease_name", "doid_name")), None)
    eid_col = next((c for c in df.columns if c.lower() in ("id", "evidence_id", "eid")), None)

    for _, row in rows.iterrows():
        drug_raw = str(row[therapy_col]) if therapy_col else "Unknown"
        if not drug_raw or drug_raw.lower() in ("nan", "none", ""):
            continue

        level_raw = str(row[level_col]) if level_col else "E"
        sig_raw = str(row[sig_col]) if sig_col else ""
        stmt = str(row[stmt_col]) if stmt_col else None
        pmid = str(row[pmid_col]) if pmid_col else None
        var_name = str(row[var_col2]) if var_col2 else None
        disease_name = str(row[ct_col2]) if ct_col2 else None
        eid = str(row[eid_col]) if eid_col else None

        evidence.append(
            KBEvidence(
                source_kb=SourceKB.CIVIC,
                source_id=f"EID{eid}" if eid else None,
                gene=gene,
                variant=var_name,
                disease=disease_name,
                drug=drug_raw,
                response_type=_map_significance(sig_raw),
                evidence_level_raw=level_raw,
                tier=_map_tier(level_raw),
                evidence_statement=stmt,
                pmids=[pmid] if pmid and pmid.lower() not in ("nan", "none") else [],
                source_url=f"https://civicdb.org/evidence/{eid}" if eid else "https://civicdb.org",
            )
        )

    return evidence


# ── GraphQL → KBEvidence conversion ──────────────────────────────────────────

def _graphql_to_evidence(gene: str, data: Dict[str, Any]) -> List[KBEvidence]:
    """Convert raw CIViC GraphQL response to KBEvidence list."""
    evidence: List[KBEvidence] = []
    genes = (data.get("genes") or {}).get("nodes") or []

    for g in genes:
        variants = (g.get("variants") or {}).get("nodes") or []
        for v in variants:
            var_name = v.get("name", "")
            mps = (v.get("molecularProfiles") or {}).get("nodes") or []
            for mp in mps:
                items = (mp.get("evidenceItems") or {}).get("nodes") or []
                for item in items:
                    therapies = item.get("therapies") or []
                    for therapy in therapies:
                        drug_name = therapy.get("name", "")
                        if not drug_name:
                            continue

                        source = item.get("source") or {}
                        pmid = str(source.get("pmid", "")) or None
                        trials = [
                            ct.get("nctId", "")
                            for ct in (source.get("clinicalTrials") or [])
                        ]
                        disease = item.get("disease") or {}

                        evidence.append(
                            KBEvidence(
                                source_kb=SourceKB.CIVIC,
                                source_id=f"EID{item.get('id', '')}",
                                gene=gene,
                                variant=var_name or None,
                                disease=disease.get("name"),
                                disease_ontology_id=f"DOID:{disease.get('doid', '')}"
                                if disease.get("doid")
                                else None,
                                drug=drug_name,
                                response_type=_map_significance(
                                    item.get("significanceLabel", "")
                                ),
                                evidence_level_raw=item.get("evidenceLevel", "E"),
                                tier=_map_tier(item.get("evidenceLevel", "E")),
                                evidence_statement=item.get("statement"),
                                pmids=[pmid] if pmid and pmid != "None" else [],
                                clinical_trials=trials,
                                source_raw=item,
                                source_url=f"https://civicdb.org/evidence/{item.get('id', '')}",
                            )
                        )
    return evidence


# ── Public interface ──────────────────────────────────────────────────────────

def query_civic(
    gene: str,
    variant: Optional[str] = None,
    cancer_type: Optional[str] = None,
    use_bulk: bool = True,
) -> List[KBEvidence]:
    """
    Query CIViC for evidence records matching gene/variant/cancer_type.

    Strategy:
      1. Try GraphQL (live, accurate, rate-limited)
      2. Fall back to nightly TSV cache (fast, offline-capable)
    """
    # Try GraphQL first (gives richer data)
    data = _graphql_query(gene)
    if data:
        results = _graphql_to_evidence(gene, data)
        if results:
            # Apply variant/cancer_type filters post-fetch
            if variant:
                results = [
                    r for r in results
                    if not r.variant or variant.upper() in (r.variant or "").upper()
                ]
            if cancer_type:
                results = [
                    r for r in results
                    if not r.disease or cancer_type.lower() in (r.disease or "").lower()
                ]
            return results

    # Fallback: bulk TSV
    if use_bulk:
        df = _load_civic_tsv()
        return _tsv_query(gene, variant, cancer_type, df)

    return []
