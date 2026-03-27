"""
Cancer Genome Interpreter (CGI) Biomarkers client.

Data source: CC0 bulk TSV from https://www.cancergenomeinterpreter.org/data/
Latest file: cgi_biomarkers_latest.tsv

CGI evidence levels:
  "FDA approved"  → Tier A
  "European EMA"  → Tier A
  "NCCN guidelines" → Tier A
  "Clinical trials Phase III" → Tier B
  "Clinical trials Phase II" → Tier B
  "Clinical trials Phase I" → Tier C
  "Pre-clinical" → Tier D
  "Case report" → Tier C  (single patient but observed in humans)
  "inferential" → Tier D

CGI columns (from the 2022 release):
  Alteration, Alteration type, Assay type, Biomarker,
  COAD, Combination, Disease, Drug family, Drug full name,
  Evidence level, Gene, Source, Targeting, Response type
"""
from __future__ import annotations

import io
import logging
from pathlib import Path
from typing import List, Optional

import httpx
import pandas as pd

from ..models import EvidenceTier, KBEvidence, ResponseType, SourceKB

logger = logging.getLogger(__name__)

CGI_BIOMARKERS_URL = "https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_latest.tsv"

_CACHE_DIR = Path(".cache/cgi")
_CACHE_DIR.mkdir(parents=True, exist_ok=True)
_CACHE_FILE = _CACHE_DIR / "biomarkers.parquet"


# ── Level → Tier ──────────────────────────────────────────────────────────────

def _map_tier(level: str) -> EvidenceTier:
    l = level.lower().strip()
    if any(k in l for k in ("fda", "ema", "nccn", "approved")):
        return EvidenceTier.A
    if "phase iii" in l or "phase 3" in l:
        return EvidenceTier.B
    if "phase ii" in l or "phase 2" in l:
        return EvidenceTier.B
    if "phase i" in l or "phase 1" in l:
        return EvidenceTier.C
    if "case report" in l:
        return EvidenceTier.C
    if "pre-clinical" in l or "preclinical" in l:
        return EvidenceTier.D
    if "inferential" in l:
        return EvidenceTier.D
    return EvidenceTier.UNKNOWN


def _map_response(rtype: str) -> ResponseType:
    r = rtype.lower()
    if any(k in r for k in ("sensitive", "response", "responsive", "benefit", "activity")):
        return ResponseType.SENSITIVITY
    if "resist" in r:
        return ResponseType.RESISTANCE
    if "reduced" in r:
        return ResponseType.REDUCED_SENSITIVITY
    return ResponseType.UNKNOWN


# ── Bulk TSV loader ───────────────────────────────────────────────────────────

def _load_cgi_df(force_refresh: bool = False) -> pd.DataFrame:
    if _CACHE_FILE.exists() and not force_refresh:
        logger.info("CGI: loading from cache")
        return pd.read_parquet(_CACHE_FILE)

    logger.info("CGI: downloading biomarkers TSV …")
    try:
        resp = httpx.get(CGI_BIOMARKERS_URL, follow_redirects=True, timeout=120)
        resp.raise_for_status()
        df = pd.read_csv(io.BytesIO(resp.content), sep="\t", low_memory=False)
        df.to_parquet(_CACHE_FILE)
        logger.info("CGI: cached %d biomarker rows", len(df))
        return df
    except Exception as exc:
        logger.warning("CGI download failed: %s", exc)
        return pd.DataFrame()


# ── Query ─────────────────────────────────────────────────────────────────────

def query_cgi(
    gene: str,
    variant: Optional[str] = None,
    cancer_type: Optional[str] = None,
    force_refresh: bool = False,
) -> List[KBEvidence]:
    """Return CGI biomarker records for a gene/variant/cancer_type."""
    df = _load_cgi_df(force_refresh)
    if df.empty:
        return []

    # Normalise column names (CGI has changed schema across versions)
    cols_lower = {c.lower(): c for c in df.columns}

    gene_col = cols_lower.get("gene") or cols_lower.get("biomarker") or next(
        (c for c in df.columns if "gene" in c.lower()), None
    )
    if gene_col is None:
        logger.warning("CGI: cannot find gene column in %s", list(df.columns))
        return []

    mask = df[gene_col].str.upper() == gene.upper()

    # Variant filter on Alteration or Biomarker column
    for vc in ("alteration", "biomarker"):
        vc_real = cols_lower.get(vc)
        if vc_real and variant:
            mask &= df[vc_real].str.contains(variant, case=False, na=False)
            break

    # Cancer type filter on Disease column
    disease_col = cols_lower.get("disease") or cols_lower.get("tumor type")
    if disease_col and cancer_type:
        mask &= df[disease_col].str.contains(cancer_type, case=False, na=False)

    rows = df[mask]
    evidence: List[KBEvidence] = []

    drug_col = cols_lower.get("drug full name") or cols_lower.get("drug") or cols_lower.get("targeting")
    family_col = cols_lower.get("drug family")
    level_col = cols_lower.get("evidence level") or cols_lower.get("evidence_level")
    response_col = cols_lower.get("response type") or cols_lower.get("response_type")
    source_col = cols_lower.get("source")
    alt_col = cols_lower.get("alteration")
    alt_type_col = cols_lower.get("alteration type")

    for _, row in rows.iterrows():
        drug_name = str(row[drug_col]) if drug_col else "Unknown"
        if not drug_name or drug_name.lower() in ("nan", "none", ""):
            continue

        level_raw = str(row[level_col]) if level_col else ""
        resp_raw = str(row[response_col]) if response_col else ""
        disease_name = str(row[disease_col]) if disease_col else None
        drug_family = str(row[family_col]) if family_col else None
        source = str(row[source_col]) if source_col else None
        alt = str(row[alt_col]) if alt_col else None
        alt_type = str(row[alt_type_col]) if alt_type_col else None

        # Extract PMIDs from source field (CGI sometimes lists PMIDs there)
        pmids = []
        if source:
            import re
            pmids = re.findall(r"\b\d{7,8}\b", source)

        evidence.append(
            KBEvidence(
                source_kb=SourceKB.CGI,
                gene=gene,
                variant=alt or variant,
                variant_type=alt_type,
                disease=disease_name if disease_name and disease_name.lower() != "nan" else None,
                drug=drug_name,
                drug_class=drug_family if drug_family and drug_family.lower() != "nan" else None,
                response_type=_map_response(resp_raw),
                evidence_level_raw=level_raw,
                tier=_map_tier(level_raw),
                pmids=pmids,
                source_url="https://www.cancergenomeinterpreter.org/biomarkers",
            )
        )

    return evidence
