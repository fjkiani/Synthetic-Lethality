"""
BrM extension for the multimodal layer.

Thin adapter: imports BrMQueryInput and build_brm_matrix from sl_agent.brm.
Exposes run_brm_analysis() for use by the API layer.

No modifications to existing multimodal files.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from sl_agent.brm import BrMQueryInput, BrMTargetabilityMatrix, build_brm_matrix


def run_brm_analysis(
    query: Optional[BrMQueryInput] = None,
    expression_df: Optional[pd.DataFrame] = None,
    crispr_df: Optional[pd.DataFrame] = None,
    sample_info: Optional[pd.DataFrame] = None,
) -> BrMTargetabilityMatrix:
    """
    Run BrM targetability analysis.

    Parameters
    ----------
    query : BrMQueryInput, optional
        Query parameters. Defaults to BrMQueryInput() (7 panel genes, NSCLC).
    expression_df : pd.DataFrame, optional
        DepMap expression matrix. None → non-blocking (receipt-only mode).
    crispr_df : pd.DataFrame, optional
        DepMap CRISPR scores. None → non-blocking.
    sample_info : pd.DataFrame, optional
        DepMap sample info. None → non-blocking.

    Returns
    -------
    BrMTargetabilityMatrix
        3-section partitioned matrix (CALIBRATION / NOVEL / NEGATIVE).
        Always includes RUO envelope and disclaimer.
    """
    if query is None:
        query = BrMQueryInput()

    return build_brm_matrix(
        query=query,
        expression_df=expression_df,
        crispr_df=crispr_df,
        sample_info=sample_info,
    )
