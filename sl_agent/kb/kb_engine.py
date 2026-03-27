"""
KB Engine — orchestrates all open KB clients.

This replaces OncoKB in the drug mapper. Same interface, zero tokens required.

Entry points:
  query(gene, variant, cancer_type)  → List[DrugRecommendation]
  explain(gene, drug)                → List[KBEvidence] raw
  oncokb_drop_in(gene)               → List[ExternalDrugEvidence]  (for drug_mapper.py)
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional

from .clients.cgi_client import query_cgi
from .clients.civic_client import query_civic
from .clients.clinvar_client import get_variant_oncogenicity
from .clients.jax_client import query_jax
from .evidence_fuser import fuse_evidence
from .models import (
    DrugRecommendation,
    EvidenceTier,
    KBEvidence,
    KBQueryInput,
    KBQueryResponse,
    SourceKB,
)

# ChEMBL is still used for structural drug info (not therapy evidence)
# so we keep it as a supplementary source

logger = logging.getLogger(__name__)

_MAX_WORKERS = 4


# ── Internal: drop-in type for drug_mapper.py compatibility ──────────────────

class OpenKBEvidence:
    """
    Thin shim that mimics the ExternalDrugEvidence interface from the old
    OncoKB layer so drug_mapper.py needs zero changes.
    """
    def __init__(self, rec: DrugRecommendation):
        self.source = "+".join(kb.value for kb in rec.supporting_kbs)
        self.drug_name = rec.drug
        self.drug_chembl_id = None
        self.mechanism_of_action = rec.drug_class
        self.drug_class = _map_drug_class(rec.drug_class or "")
        self.primary_targets = [rec.gene]
        self.max_phase = _tier_to_phase(rec.tier)
        self.oncokb_level = None   # No OncoKB level — use tier instead
        self.mars_tier = rec.tier.value
        self.confidence = rec.confidence
        self.supporting_kbs = [kb.value for kb in rec.supporting_kbs]
        self.pmids = rec.pmids
        self.clinical_trials = rec.clinical_trials


def _tier_to_phase(tier: EvidenceTier) -> int:
    """Convert Mars tier to a max_phase integer for backward compat with rank scoring."""
    return {
        EvidenceTier.A: 4,
        EvidenceTier.B: 3,
        EvidenceTier.C: 2,
        EvidenceTier.D: 0,
        EvidenceTier.UNKNOWN: 0,
    }.get(tier, 0)


def _map_drug_class(class_str: str):
    """Re-import the same inference function from drug_mapper to keep consistent classes."""
    try:
        from ..core.drug_mapper import _infer_drug_class
        return _infer_drug_class(class_str or "")
    except ImportError:
        return None


# ── Parallel KB query ─────────────────────────────────────────────────────────

def _run_all_sources(
    gene: str,
    variant: Optional[str],
    cancer_type: Optional[str],
) -> List[KBEvidence]:
    """Run CIViC, CGI, and JAX in parallel, merge results."""
    all_evidence: List[KBEvidence] = []
    tasks = {
        "civic": lambda: query_civic(gene, variant, cancer_type),
        "cgi": lambda: query_cgi(gene, variant, cancer_type),
        "jax": lambda: query_jax(gene, variant, cancer_type),
    }

    with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as pool:
        futures = {pool.submit(fn): name for name, fn in tasks.items()}
        for future in as_completed(futures):
            name = futures[future]
            try:
                results = future.result(timeout=30)
                logger.info("KB %s returned %d records for %s", name, len(results), gene)
                all_evidence.extend(results)
            except Exception as exc:
                logger.warning("KB %s failed for %s: %s", name, gene, exc)

    return all_evidence


# ── Public API ────────────────────────────────────────────────────────────────

def query(
    gene: str,
    variant: Optional[str] = None,
    cancer_type: Optional[str] = None,
    min_tier: EvidenceTier = EvidenceTier.D,
    include_resistance: bool = True,
    top_n: int = 20,
) -> KBQueryResponse:
    """
    Full open-KB query. Replaces OncoKB.

    Returns fused, tiered drug recommendations from CIViC + CGI + JAX.
    """
    q = KBQueryInput(
        gene=gene,
        variant=variant,
        cancer_type=cancer_type,
        min_tier=min_tier,
        include_resistance=include_resistance,
        top_n=top_n,
    )

    raw = _run_all_sources(gene, variant, cancer_type)
    recommendations = fuse_evidence(raw, q)

    return KBQueryResponse(
        query=q,
        recommendations=recommendations,
        sources_queried=[SourceKB.CIVIC, SourceKB.CGI, SourceKB.JAX],
        total_evidence_items=len(raw),
        notes=f"Open KB stack: CIViC + CGI + JAX-CKB. No OncoKB token required.",
    )


def explain(gene: str, drug: str) -> List[KBEvidence]:
    """
    Return raw evidence records from all KBs for a specific gene/drug pair.
    Used by /kb/explain endpoint.
    """
    raw = _run_all_sources(gene, None, None)
    drug_lower = drug.lower()
    return [e for e in raw if drug_lower in e.drug.lower()]


def oncokb_drop_in(gene: str, variant: Optional[str] = None) -> List[OpenKBEvidence]:
    """
    Drop-in replacement for the oncokb_gene_drugs() function in drug_mapper.py.
    Returns a list of OpenKBEvidence objects with the same interface as ExternalDrugEvidence.

    Usage in drug_mapper.py:
        # OLD: oncokb_evidence = oncokb_gene_drugs(gene)
        # NEW:
        from ..kb.kb_engine import oncokb_drop_in
        oncokb_evidence = oncokb_drop_in(gene, variant)
    """
    response = query(gene, variant=variant, min_tier=EvidenceTier.C)
    return [OpenKBEvidence(rec) for rec in response.recommendations]
