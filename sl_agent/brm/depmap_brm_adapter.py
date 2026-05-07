"""
DepMap BrM adapter — non-blocking expression and co-dependency queries.

Returns None gracefully if DepMap data is not loaded.
Lineage filter: OncotreeLineage column, "Non-Small Cell Lung Cancer".
Falls back to pan-cancer if fewer than 5 lung lines available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class BrMDepmapNote:
    gene: str
    lineage: str
    mean_expression: Optional[float]
    expression_percentile_in_lineage: Optional[float]
    n_lines_in_lineage: int
    summary: str


def query_expression(
    gene: str,
    expression_df: Optional[pd.DataFrame],
    sample_info: Optional[pd.DataFrame],
    lineage: str = "lung",
    lineage_col: str = "OncotreeLineage",
    nsclc_value: str = "Non-Small Cell Lung Cancer",
) -> Optional[BrMDepmapNote]:
    """
    Query expression for a gene in NSCLC lines.

    Returns None if expression_df or sample_info is None (non-blocking).
    """
    if expression_df is None or sample_info is None:
        return None

    if gene not in expression_df.columns:
        return BrMDepmapNote(
            gene=gene,
            lineage=lineage,
            mean_expression=None,
            expression_percentile_in_lineage=None,
            n_lines_in_lineage=0,
            summary=f"{gene} not found in expression matrix.",
        )

    # NSCLC lines
    nsclc_ids = sample_info[
        sample_info[lineage_col] == nsclc_value
    ].index.tolist()
    common_ids = list(set(nsclc_ids) & set(expression_df.index))

    # fallback to pan-cancer if too few
    if len(common_ids) < 5:
        common_ids = list(expression_df.index)
        lineage_used = "pan-cancer (fallback)"
    else:
        lineage_used = "NSCLC"

    expr = expression_df.loc[common_ids, gene].dropna()
    if len(expr) == 0:
        return BrMDepmapNote(
            gene=gene,
            lineage=lineage_used,
            mean_expression=None,
            expression_percentile_in_lineage=None,
            n_lines_in_lineage=len(common_ids),
            summary=f"{gene}: no expression values in {lineage_used} lines.",
        )

    mean_expr = float(expr.mean())
    # percentile within lineage
    all_expr = expression_df.loc[common_ids].values.flatten()
    all_expr = all_expr[~np.isnan(all_expr)]
    pct = float(np.mean(all_expr <= mean_expr) * 100) if len(all_expr) > 0 else None

    summary = (
        f"{gene} mean expression in {lineage_used} (n={len(expr)}): "
        f"{mean_expr:.2f} log2TPM+1"
        + (f", {pct:.0f}th percentile" if pct is not None else "")
        + "."
    )

    return BrMDepmapNote(
        gene=gene,
        lineage=lineage_used,
        mean_expression=mean_expr,
        expression_percentile_in_lineage=pct,
        n_lines_in_lineage=len(common_ids),
        summary=summary,
    )


def query_codependency(
    gene: str,
    crispr_df: Optional[pd.DataFrame],
    sample_info: Optional[pd.DataFrame],
    lineage_col: str = "OncotreeLineage",
    nsclc_value: str = "Non-Small Cell Lung Cancer",
    top_n: int = 10,
) -> Optional[dict]:
    """
    Compute Pearson co-dependency correlations for a gene in NSCLC lines.

    Returns dict {partner_gene: r_value} for top_n partners, or None if unavailable.
    """
    if crispr_df is None or sample_info is None:
        return None

    if gene not in crispr_df.columns:
        return None

    nsclc_ids = sample_info[
        sample_info[lineage_col] == nsclc_value
    ].index.tolist()
    common_ids = list(set(nsclc_ids) & set(crispr_df.index))

    if len(common_ids) < 10:
        return None

    sub = crispr_df.loc[common_ids].dropna(axis=1, thresh=int(len(common_ids) * 0.5))
    if gene not in sub.columns:
        return None

    gene_dep = sub[gene]
    other_genes = [g for g in sub.columns if g != gene]

    corrs = {}
    for g in other_genes:
        partner_dep = sub[g]
        mask = gene_dep.notna() & partner_dep.notna()
        if mask.sum() < 10:
            continue
        r = float(np.corrcoef(gene_dep[mask], partner_dep[mask])[0, 1])
        if not np.isnan(r):
            corrs[g] = r

    # return top_n by |r|
    top = sorted(corrs.items(), key=lambda x: abs(x[1]), reverse=True)[:top_n]
    return dict(top)
