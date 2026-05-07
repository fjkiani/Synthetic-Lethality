"""
BrM targetability matrix builder — 3-section partitioned output.

Sections:
  CALIBRATION: receipt-backed, literature-validated (BACE1, MMP9)
  NOVEL:       DepMap-discovered or receipt-only without DepMap signal
  NEGATIVE:    tested, no signal at FDR<0.25

BACE1 is always rows[0] (gold standard invariant).
Within each section, rows sorted by confidence_score descending.
"""

from __future__ import annotations

import datetime
import json
from typing import Dict, List, Optional

import pandas as pd

from sl_agent.brm.models import (
    BrMGeneRow,
    BrMQueryInput,
    BrMRowClass,
    BrMTargetabilityMatrix,
    ExploitMode,
)
from sl_agent.brm.brm_receipts import get_brm_receipt
from sl_agent.brm.exploit_router import route
from sl_agent.brm.expression_stratifier import run_expression_stratified_sl
from sl_agent.brm.depmap_brm_adapter import query_expression, query_codependency


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

def build_brm_matrix(
    query: BrMQueryInput,
    expression_df: Optional[pd.DataFrame] = None,
    crispr_df: Optional[pd.DataFrame] = None,
    sample_info: Optional[pd.DataFrame] = None,
) -> BrMTargetabilityMatrix:
    """
    Build the BrM targetability matrix for the given query.

    Parameters
    ----------
    query : BrMQueryInput
    expression_df : pd.DataFrame, optional
        DepMap expression (rows=cell lines, cols=genes). None → non-blocking.
    crispr_df : pd.DataFrame, optional
        DepMap CRISPR scores. None → non-blocking.
    sample_info : pd.DataFrame, optional
        DepMap sample info with OncotreeLineage. None → non-blocking.
    """
    now = datetime.datetime.utcnow().isoformat() + "+00:00"
    rows: List[BrMGeneRow] = []

    for gene in query.genes:
        # Layer 1: expression-stratified SL
        strat = None
        if (
            query.include_depmap_expression
            and expression_df is not None
            and crispr_df is not None
            and sample_info is not None
        ):
            strat = run_expression_stratified_sl(
                gene=gene,
                expression_df=expression_df,
                crispr_df=crispr_df,
                sample_info=sample_info,
            )

        # Layer 2: co-dependency
        codep = None
        if (
            query.include_depmap_expression
            and crispr_df is not None
            and sample_info is not None
        ):
            codep = query_codependency(
                gene=gene,
                crispr_df=crispr_df,
                sample_info=sample_info,
            )

        # DepMap expression note
        depmap_note = None
        if (
            query.include_depmap_expression
            and expression_df is not None
            and sample_info is not None
        ):
            note_obj = query_expression(
                gene=gene,
                expression_df=expression_df,
                sample_info=sample_info,
            )
            if note_obj is not None:
                depmap_note = note_obj.summary

        # Live literature (non-blocking)
        live_lit_note = None
        if query.include_live_literature:
            try:
                from sl_agent.brm.live_literature import query_brm_literature
                live_lit_note = query_brm_literature(gene, query.cancer_type)
            except Exception:
                pass

        # Route through exploit engine
        row = route(
            gene=gene,
            strat=strat,
            codep=codep,
            depmap_note=depmap_note,
            live_lit_note=live_lit_note,
        )
        rows.append(row)

    # --- 3-section partitioning ---
    calibration_rows = [r for r in rows if r.row_class == BrMRowClass.CALIBRATION]
    novel_rows       = [r for r in rows if r.row_class == BrMRowClass.NOVEL]
    negative_rows    = [r for r in rows if r.row_class == BrMRowClass.NEGATIVE]

    # sort within sections by confidence descending
    calibration_rows.sort(key=lambda r: r.confidence_score, reverse=True)
    novel_rows.sort(key=lambda r: r.confidence_score, reverse=True)
    negative_rows.sort(key=lambda r: r.confidence_score, reverse=True)

    # BACE1 invariant: always first in calibration (and thus rows[0])
    bace1_rows = [r for r in calibration_rows if r.gene == "BACE1"]
    other_cal  = [r for r in calibration_rows if r.gene != "BACE1"]
    calibration_rows = bace1_rows + other_cal

    # full ordered list: CALIBRATION → NOVEL → NEGATIVE
    ordered_rows = calibration_rows + novel_rows + negative_rows

    # --- summaries ---
    exploit_summary: Dict[str, int] = {}
    for r in ordered_rows:
        k = r.primary_exploit_mode.value
        exploit_summary[k] = exploit_summary.get(k, 0) + 1

    cascade_summary: Dict[str, int] = {}
    for r in ordered_rows:
        k = r.cascade_step.value
        cascade_summary[k] = cascade_summary.get(k, 0) + 1

    row_class_summary: Dict[str, int] = {
        "calibration": len(calibration_rows),
        "novel":       len(novel_rows),
        "negative":    len(negative_rows),
    }

    return BrMTargetabilityMatrix(
        query_context=f"{query.cancer_type} {query.context.replace('_', ' ')}",
        cancer_type=query.cancer_type,
        rows=ordered_rows,
        calibration_rows=calibration_rows,
        novel_rows=novel_rows,
        negative_rows=negative_rows,
        exploit_mode_summary=exploit_summary,
        cascade_step_summary=cascade_summary,
        row_class_summary=row_class_summary,
        frozen_at=now,
    )


# ---------------------------------------------------------------------------
# Freeze artifact
# ---------------------------------------------------------------------------

def freeze_matrix(matrix: BrMTargetabilityMatrix, path: str) -> None:
    """Serialize matrix to JSON artifact (manuscript discipline)."""
    data = matrix.model_dump()
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
