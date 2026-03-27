"""
KB API routes — open clinical knowledge base endpoints.

GET  /kb/query   — fused drug recommendations for gene/variant/cancer_type
GET  /kb/explain — raw evidence records for a gene+drug pair
GET  /kb/sources — list supported source KBs and their status
"""
from __future__ import annotations

import logging
from typing import Optional

from fastapi import APIRouter, HTTPException, status

from ..kb.kb_engine import explain, query
from ..kb.models import (
    EvidenceTier,
    KBExplainResponse,
    KBQueryInput,
    KBQueryResponse,
    SourceKB,
)

logger = logging.getLogger(__name__)
kb_router = APIRouter(prefix="/kb", tags=["Knowledge Base"])


@kb_router.get("/sources")
async def list_sources():
    """List supported open KB sources and their access status."""
    return {
        "sources": [
            {
                "name": "CIViC",
                "url": "https://civicdb.org",
                "license": "CC0",
                "access": "GraphQL API + nightly TSV",
                "token_required": False,
                "coverage": "3200+ variants, 470+ genes, crowdsourced clinical evidence",
            },
            {
                "name": "CGI",
                "url": "https://www.cancergenomeinterpreter.org",
                "license": "CC0",
                "access": "Bulk TSV download",
                "token_required": False,
                "coverage": "Cancer biomarkers database, FDA + clinical trial evidence",
            },
            {
                "name": "JAX-CKB",
                "url": "https://ckbhome.jax.org",
                "license": "Public core",
                "access": "Curated seed (top 15 Mars-priority genes)",
                "token_required": False,
                "coverage": "High-specificity Tier A/B for BRCA1/2, EGFR, ALK, BRAF, KRAS, MBD4, etc.",
            },
            {
                "name": "ChEMBL",
                "url": "https://www.ebi.ac.uk/chembl",
                "license": "CC BY-SA 3.0",
                "access": "REST API",
                "token_required": False,
                "coverage": "Bioactivity, MoA, clinical phase for all small molecules",
            },
            {
                "name": "ClinVar",
                "url": "https://www.ncbi.nlm.nih.gov/clinvar",
                "license": "Public domain",
                "access": "eutils API (filter layer)",
                "token_required": False,
                "coverage": "Variant oncogenicity/pathogenicity classifications",
            },
        ],
        "onkokb_status": "NOT_USED — replaced by open stack",
    }


@kb_router.get("/query", response_model=KBQueryResponse)
async def kb_query(
    gene: str,
    variant: Optional[str] = None,
    cancer_type: Optional[str] = None,
    min_tier: EvidenceTier = EvidenceTier.D,
    include_resistance: bool = True,
    top_n: int = 20,
):
    """
    Query the open KB stack for drug recommendations.

    - **gene**: Hugo gene symbol (e.g. BRCA1, MBD4, KRAS)
    - **variant**: Optional variant context (e.g. V600E, LOF, amplification)
    - **cancer_type**: Optional cancer type filter (e.g. "Ovarian Cancer")
    - **min_tier**: Minimum evidence tier to return (A/B/C/D)
    - **include_resistance**: Include resistance evidence
    - **top_n**: Max recommendations to return

    Sources: CIViC + CGI + JAX-CKB (no OncoKB token required)
    """
    try:
        response = query(
            gene=gene,
            variant=variant,
            cancer_type=cancer_type,
            min_tier=min_tier,
            include_resistance=include_resistance,
            top_n=top_n,
        )
        return response
    except Exception as exc:
        logger.exception("KB query failed for gene=%s", gene)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(exc),
        )


@kb_router.get("/explain", response_model=KBExplainResponse)
async def kb_explain(gene: str, drug: str):
    """
    Return raw evidence records from all KBs for a specific gene + drug pair.

    Shows the receipts — every evidence statement, PMID, trial ID,
    and source KB that supports or contradicts this combination.
    """
    try:
        raw = explain(gene=gene, drug=drug)
        # Build a summary recommendation if we have evidence
        rec = None
        if raw:
            from ..kb.evidence_fuser import fuse_evidence
            from ..kb.models import KBQueryInput
            q = KBQueryInput(gene=gene, top_n=1)
            recs = fuse_evidence(raw, q)
            rec = recs[0] if recs else None

        return KBExplainResponse(
            gene=gene,
            drug=drug,
            raw_evidence=raw,
            recommendation=rec,
        )
    except Exception as exc:
        logger.exception("KB explain failed for gene=%s drug=%s", gene, drug)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(exc),
        )
