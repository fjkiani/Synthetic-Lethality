"""
Expression-stratified SL discovery for BrM universe genes.

Method:
  1. For each BrM gene, split NSCLC lines into Q75 (high) vs Q25 (low) expression quartiles
  2. Mann-Whitney U test on CRISPR dependency scores of SL candidate genes
  3. BH correction across all tests per BrM gene
  4. FDR-weighted novelty score (v2 formula, POL-002)

v2 novelty score formula (POL-002):
  weighted_n = min((n_fdr10 * 2 + n_fdr10_25 * 1) / 50, 1.0)
  score = weighted_n * 0.40 + best_delta_score * 0.20 + codep_score * 0.30 + fdr_quality * 0.10

Returns SLPartner objects for integration into BrMGeneRow.sl_partners.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from sl_agent.brm.models import SLPartner


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FDR_STRONG   = 0.10   # strong tier threshold
FDR_MODERATE = 0.25   # moderate tier threshold
DELTA_STRONG = 0.30   # |delta| threshold for strong tier
DELTA_MIN    = 0.10   # minimum |delta| to report
Q_HIGH       = 0.75   # upper quartile cutoff
Q_LOW        = 0.25   # lower quartile cutoff
MIN_LINES    = 5      # minimum lines per group


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class StratificationResult:
    """Full result for one BrM gene."""
    gene: str
    n_nsclc_lines: int
    n_high: int
    n_low: int
    partners: List[SLPartner] = field(default_factory=list)
    novelty_score: float = 0.0
    novelty_score_v1: float = 0.0   # for sensitivity comparison
    novelty_score_v3: float = 0.0   # delta-tiered variant
    codep_score: float = 0.0
    best_delta: float = 0.0
    n_fdr10: int = 0
    n_fdr10_25: int = 0
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Core stratification
# ---------------------------------------------------------------------------

def run_expression_stratified_sl(
    gene: str,
    expression_df: pd.DataFrame,
    crispr_df: pd.DataFrame,
    sample_info: pd.DataFrame,
    lineage_col: str = "OncotreeLineage",
    lineage_value: str = "Non-Small Cell Lung Cancer",
    codep_scores: Optional[Dict[str, float]] = None,
) -> StratificationResult:
    """
    Run expression-stratified SL for one BrM gene.

    Parameters
    ----------
    gene : str
        BrM universe gene to stratify on.
    expression_df : pd.DataFrame
        Rows = cell lines (DepMap IDs), columns = gene symbols.
    crispr_df : pd.DataFrame
        Rows = cell lines (DepMap IDs), columns = gene symbols (CRISPR scores).
    sample_info : pd.DataFrame
        Must contain lineage_col column; index = DepMap IDs.
    codep_scores : dict, optional
        Pre-computed co-dependency scores {partner_gene: r_value}.
    """
    result = StratificationResult(gene=gene, n_nsclc_lines=0, n_high=0, n_low=0)

    # --- filter to NSCLC lines ---
    nsclc_ids = sample_info[
        sample_info[lineage_col] == lineage_value
    ].index.tolist()

    # intersect with expression and crispr
    expr_ids = set(expression_df.index)
    crispr_ids = set(crispr_df.index)
    common_ids = list(set(nsclc_ids) & expr_ids & crispr_ids)

    result.n_nsclc_lines = len(common_ids)

    if gene not in expression_df.columns:
        result.error = f"{gene} not in expression matrix"
        return result

    if len(common_ids) < MIN_LINES * 2:
        result.error = f"Insufficient NSCLC lines: {len(common_ids)}"
        return result

    expr = expression_df.loc[common_ids, gene]
    q_high_val = expr.quantile(Q_HIGH)
    q_low_val  = expr.quantile(Q_LOW)

    high_ids = expr[expr >= q_high_val].index.tolist()
    low_ids  = expr[expr <= q_low_val].index.tolist()

    result.n_high = len(high_ids)
    result.n_low  = len(low_ids)

    if result.n_high < MIN_LINES or result.n_low < MIN_LINES:
        result.error = (
            f"Insufficient lines after quartile split: "
            f"high={result.n_high}, low={result.n_low}"
        )
        return result

    # --- test all CRISPR genes ---
    crispr_sub = crispr_df.loc[common_ids]
    candidate_genes = [g for g in crispr_sub.columns if g != gene]

    pvals: List[float] = []
    deltas: List[float] = []
    tested_genes: List[str] = []

    for cg in candidate_genes:
        dep_high = crispr_sub.loc[high_ids, cg].dropna()
        dep_low  = crispr_sub.loc[low_ids,  cg].dropna()
        if len(dep_high) < MIN_LINES or len(dep_low) < MIN_LINES:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, p = stats.mannwhitneyu(dep_high, dep_low, alternative="two-sided")
        delta = dep_high.mean() - dep_low.mean()
        pvals.append(p)
        deltas.append(delta)
        tested_genes.append(cg)

    if not tested_genes:
        result.error = "No testable CRISPR genes"
        return result

    # BH correction
    _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")

    # collect significant partners
    partners: List[SLPartner] = []
    for cg, delta, fdr in zip(tested_genes, deltas, fdrs):
        if fdr <= FDR_MODERATE and abs(delta) >= DELTA_MIN and delta < 0:
            tier = (
                "strong"   if fdr <= FDR_STRONG and abs(delta) >= DELTA_STRONG else
                "moderate" if fdr <= FDR_MODERATE and abs(delta) >= DELTA_MIN else
                "weak"
            )
            partners.append(SLPartner(
                partner_gene=cg,
                delta_dep=round(delta, 4),
                fdr=round(float(fdr), 6),
                n_high=len(high_ids),
                n_low=len(low_ids),
                partner_quality_tier=tier,
            ))

    partners.sort(key=lambda p: (p.fdr, p.delta_dep))
    result.partners = partners

    # --- novelty score components ---
    n_fdr10    = sum(1 for p in partners if p.fdr <= FDR_STRONG)
    n_fdr10_25 = sum(1 for p in partners if FDR_STRONG < p.fdr <= FDR_MODERATE)
    result.n_fdr10    = n_fdr10
    result.n_fdr10_25 = n_fdr10_25

    best_delta = max((abs(p.delta_dep) for p in partners), default=0.0)
    result.best_delta = best_delta

    # co-dependency score (from pre-computed scores or 0)
    codep = 0.0
    if codep_scores and gene in codep_scores:
        r = codep_scores[gene]
        codep = max(0.0, min(1.0, abs(r)))
    result.codep_score = codep

    # v2 formula (POL-002)
    weighted_n      = min((n_fdr10 * 2 + n_fdr10_25 * 1) / 50.0, 1.0)
    best_delta_score = min(best_delta / 0.5, 1.0)
    fdr_quality      = min(n_fdr10 / 5.0, 1.0)
    result.novelty_score = (
        weighted_n       * 0.40 +
        best_delta_score * 0.20 +
        codep            * 0.30 +
        fdr_quality      * 0.10
    )

    # v1 formula (for sensitivity comparison)
    n_partners = len(partners)
    result.novelty_score_v1 = (
        min(n_partners / 10.0, 1.0) * 0.40 +
        best_delta_score             * 0.20 +
        codep                        * 0.30 +
        fdr_quality                  * 0.10
    )

    # v3 formula (delta-tiered, POL-003 sensitivity)
    n_strong   = sum(1 for p in partners if p.partner_quality_tier == "strong")
    n_moderate = sum(1 for p in partners if p.partner_quality_tier == "moderate")
    weighted_n_v3 = min((n_strong * 3 + n_moderate * 1) / 50.0, 1.0)
    result.novelty_score_v3 = (
        weighted_n_v3    * 0.40 +
        best_delta_score * 0.20 +
        codep            * 0.30 +
        fdr_quality      * 0.10
    )

    return result


def run_batch_stratification(
    genes: List[str],
    expression_df: pd.DataFrame,
    crispr_df: pd.DataFrame,
    sample_info: pd.DataFrame,
    lineage_col: str = "OncotreeLineage",
    lineage_value: str = "Non-Small Cell Lung Cancer",
    codep_scores: Optional[Dict[str, float]] = None,
) -> Dict[str, StratificationResult]:
    """Run stratification for a list of BrM genes. Returns dict keyed by gene."""
    results = {}
    for gene in genes:
        results[gene] = run_expression_stratified_sl(
            gene=gene,
            expression_df=expression_df,
            crispr_df=crispr_df,
            sample_info=sample_info,
            lineage_col=lineage_col,
            lineage_value=lineage_value,
            codep_scores=codep_scores,
        )
    return results
