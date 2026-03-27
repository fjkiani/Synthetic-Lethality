"""
Multi-modal analysis routes.

POST /analyze/multimodal  — full multi-modal evidence matrix + narrative
GET  /analyze/calibration — cytidine-analog gold-standard calibration narrative
GET  /evidence/matrix     — matrix for a gene without running CRISPR analysis
"""
from __future__ import annotations

import asyncio
import logging
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, status

from ..core.orchestrator import DataStore
from ..multimodal.matrix_builder import run_multimodal_analysis
from ..multimodal.models import MultiModalQueryInput, MultiModalResult, CandidateAxis
from ..multimodal.literature_receipts import (
    get_calibration_narrative,
    list_receipts_for_gene,
)

logger = logging.getLogger(__name__)
mm_router = APIRouter(prefix="/analyze", tags=["Multi-Modal Analysis"])
_executor = ThreadPoolExecutor(max_workers=4)


# ── Dependency ────────────────────────────────────────────────────────────────

async def verify_data_loaded():
    try:
        DataStore.ensure_loaded(require_prism=True)
    except Exception as exc:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=f"Data not yet loaded: {exc}",
        )


# ── Multi-modal analysis endpoint ─────────────────────────────────────────────

@mm_router.post(
    "/multimodal",
    response_model=MultiModalResult,
    summary="Full multi-modal SL evidence matrix + recommendation narrative",
)
async def analyze_multimodal(
    query: MultiModalQueryInput,
    _: None = Depends(verify_data_loaded),
):
    """
    Build a complete multi-modal evidence matrix for the query gene,
    integrating CRISPR, pharmacologic screens (PRISM/GDSC), literature
    receipts, and pathway mechanism evidence.

    Returns:
    - `evidence_matrix`: Full evidence grid (rows = candidate axes, columns = modalities)
    - `agreement_report`: Where CRISPR and pharmacologic data agree / disagree
    - `updated_narrative`: Tumor-board-grade interpretation using all modalities
    - `recommendation_map`: Tier mapping per axis
    - `guardrails_applied`: Doctrine guardrails triggered during analysis

    Doctrine enforced:
    - CRISPR negative + pharmacologic not examined → "Mechanistic candidate only"
      (NOT "not supported / negative")
    - Cytidine analogs for MBD4 = calibration gold standard
    - Recommendation tiers driven by multi-modal matrix, not CRISPR alone
    """
    t0 = time.perf_counter()

    # Run SL engine on the same gene to get CRISPR partner stats
    from ..core.models import SLQueryInput, MutationType
    from ..core.orchestrator import run_sl_analysis
    from ..data.depmap_loader import get_mutant_wt_lines

    loop = asyncio.get_event_loop()

    # We need CRISPR partners + raw cell line stratification
    def _run():
        # Step A: stratify cell lines
        crispr = DataStore.crispr()
        mutations = DataStore.mutations()
        sample_info = DataStore.sample_info()
        cna = DataStore.cna()
        prism = DataStore.prism()
        prism_meta = DataStore.prism_meta()

        if mutations is None:
            mutations_df = __import__("pandas").DataFrame()
        else:
            mutations_df = mutations

        try:
            mut_type = MutationType(query.mutation_type)
        except ValueError:
            mut_type = MutationType.ANY

        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=query.gene,
            mutation_df=mutations_df,
            sample_info=sample_info,
            cancer_type=query.cancer_type,
            mutation_type=mut_type.value,
            cna_df=cna,
        )

        # Step B: get CRISPR SL partners (full engine run for CRISPR modality)
        sl_partners = []
        try:
            from ..core.sl_engine import SLEngine
            from ..core.config import get_settings
            cfg = get_settings()
            engine = SLEngine(
                crispr_df=crispr,
                sample_info_df=sample_info,
                mutation_df=mutations_df,
                cna_df=cna,
                depmap_release=cfg.depmap_release,
            )
            sl_q = SLQueryInput(
                gene=query.gene,
                mutation_type=mut_type,
                cancer_type=query.cancer_type,
                top_n_partners=50,
                fdr_cutoff=0.25,
                delta_dep_cutoff=0.05,  # looser cutoff for matrix — we report all signals
            )
            _, partners, _ = engine.compute_sl_partners(
                query=sl_q,
                mutant_ids=mutant_ids,
                wt_ids=wt_ids,
            )
            sl_partners = partners
        except Exception as e:
            logger.warning("CRISPR engine failed during multimodal analysis: %s", e)

        # Step C: run multi-modal analysis
        from ..core.config import get_settings
        cfg = get_settings()

        result = run_multimodal_analysis(
            query=query,
            sl_partners=sl_partners,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
            prism_df=prism,
            prism_meta=prism_meta,
            gdsc_df=None,  # GDSC optional — wire if loaded
            sample_info=sample_info,
            depmap_release=cfg.depmap_release,
        )
        return result

    try:
        result: MultiModalResult = await loop.run_in_executor(_executor, _run)
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc))
    except RuntimeError as exc:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(exc))
    except Exception as exc:
        logger.exception("Multi-modal analysis failed")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Internal error: {exc}",
        )

    elapsed_ms = round((time.perf_counter() - t0) * 1000, 1)
    logger.info(
        "Multimodal analysis for gene=%s completed in %.1f ms", query.gene, elapsed_ms
    )
    return result


# ── Calibration endpoint ──────────────────────────────────────────────────────

@mm_router.get(
    "/calibration",
    summary="Gold-standard cytidine-analog calibration narrative for MBD4",
    response_model=None,
)
async def get_calibration(gene: str = "MBD4"):
    """
    Return the multi-modal calibration narrative for the gold-standard axis.
    Defines what 'high evidence SL' looks like and sets the bar for all
    other candidate axes.

    Default: MBD4 + cytidine analogs (npj Precis Oncol 2022, PMID 35428381).
    """
    narrative = get_calibration_narrative(gene)
    receipts = list_receipts_for_gene(gene)
    return {
        "gene": gene,
        "calibration_axis": "cytidine_analogs",
        "narrative": narrative,
        "positive_modalities_by_axis": receipts,
        "pmids": ["35428381"],
    }


# ── Evidence matrix endpoint (no live CRISPR run) ────────────────────────────

@mm_router.get(
    "/evidence_matrix",
    summary="Evidence matrix from frozen literature receipts (no live data required)",
)
async def get_evidence_matrix(
    gene: str = "MBD4",
    cancer_type: Optional[str] = None,
):
    """
    Return the evidence matrix populated from frozen literature receipts only,
    without running a live CRISPR or PRISM analysis.

    Useful for:
    - Quick audit / review without DepMap data loaded
    - Exporting the curated evidence state for a gene
    """
    query = MultiModalQueryInput(
        gene=gene,
        cancer_type=cancer_type,
        include_pharmacologic_stratification=False,
        include_literature_receipts=True,
        include_calibration_narrative=True,
    )

    from ..multimodal.matrix_builder import build_evidence_matrix
    matrix = build_evidence_matrix(query=query)

    return {
        "gene": gene,
        "cancer_type": cancer_type,
        "axes": [
            {
                "axis": row.axis.value,
                "label": row.axis_label,
                "mechanism": row.mechanism,
                "recommendation_tier": row.recommendation_tier,
                "overall_evidence_level": row.overall_evidence_level,
                "positive_modalities": row.positive_modalities(),
                "missing_modalities": row.missing_modalities(),
                "interpretation": row.interpretation,
            }
            for row in matrix.rows
        ],
        "gold_standard_summary": matrix.gold_standard_summary,
        "notes": matrix.notes,
    }
