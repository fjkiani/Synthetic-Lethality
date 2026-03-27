"""
Synthetic Lethality computation engine.

Steps performed here:
  1. Pan-essential gene filter   (noisy universally-essential genes removed)
  2. Group stratification        (mutant vs WT cell lines)
  3. Statistical testing         (Wilcoxon rank-sum + optional t-test)
  4. Multiple-testing correction (Benjamini–Hochberg FDR)
  5. Effect size                 (Cohen's d, delta-dependency)
  6. Co-dependency analysis      (Pearson r of dependency profiles)
  7. Ranking & filtering
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .models import (
    AnalysisScope,
    FrameworkSupport,
    InputContext,
    MutationType,
    SLPartner,
    SLQueryInput,
)

logger = logging.getLogger(__name__)


# ── Pan-essential blacklist (common core fitness genes) ───────────────────────
# These genes are essential in nearly ALL cell lines and produce confounding hits.
PAN_ESSENTIAL_GENES = {
    "RPL5", "RPL11", "RPS14", "RPS19",
    "POLR2A", "POLR2B", "CDC42", "RAC1",
    "ACTB", "GAPDH", "PCNA", "MCM2", "MCM7",
    "SF3B1", "SRSF1", "U2AF1",
    "CDK1", "CDK2", "CDK4",
    "KPNB1", "XPO1",
    "VCP", "HSP90AA1",
    "RB1",  # often essential via CDK4/6 confound
}

# Minimum cell lines per group for a valid test
MIN_GROUP_SIZE = 5


# ── Main engine ───────────────────────────────────────────────────────────────

class SLEngine:
    """
    Performs SL inference given pre-loaded DepMap matrices.
    All matrices must be indexed by ACH model IDs on rows.
    Gene effect matrix must have genes as columns.
    """

    def __init__(
        self,
        crispr_df: pd.DataFrame,
        sample_info_df: pd.DataFrame,
        mutation_df: Optional[pd.DataFrame] = None,
        cna_df: Optional[pd.DataFrame] = None,
        expression_df: Optional[pd.DataFrame] = None,
        depmap_release: str = "unknown",
    ):
        self.crispr = crispr_df
        self.sample_info = sample_info_df
        self.mutations = mutation_df
        self.cna = cna_df
        self.expression = expression_df
        self.release = depmap_release

        # Pre-compute pan-essential gene set from data (top 5% most essential)
        self._data_pan_essential: set = self._compute_pan_essential()

    # ── Pre-processing ─────────────────────────────────────────────────────

    def _compute_pan_essential(self, percentile: float = 1.0) -> set:
        """
        Genes where median Chronos score is in the bottom `percentile`
        across all lines → universally essential, exclude from SL.

        Uses mean of per-line medians to resist outlier-group pull.
        A gene qualifies only if its frequency of dependency
        (score < -0.5) exceeds 90% of all lines.
        """
        # Fraction of lines where gene is dependent (Chronos < -0.5)
        dep_freq = (self.crispr < -0.5).mean(axis=0)
        return set(dep_freq.index[dep_freq >= 0.90].tolist())

    def _all_pan_essential(self) -> set:
        return PAN_ESSENTIAL_GENES | self._data_pan_essential

    # ── Stratification ─────────────────────────────────────────────────────

    def _stratify(
        self,
        gene: str,
        mutant_ids: List[str],
        wt_ids: List[str],
        cancer_type: Optional[str],
    ) -> Tuple[pd.DataFrame, pd.DataFrame, AnalysisScope, str]:
        """
        Align mutant / WT IDs against available CRISPR lines.
        Falls back to pan-cancer if cancer_type subset is too small.
        Returns: (mut_dep, wt_dep, scope, note)
        """
        available = set(self.crispr.index)
        mut_avail = [m for m in mutant_ids if m in available]
        wt_avail = [w for w in wt_ids if w in available]

        scope = AnalysisScope.CANCER_SPECIFIC
        note = ""

        if len(mut_avail) < MIN_GROUP_SIZE or len(wt_avail) < MIN_GROUP_SIZE:
            if cancer_type:
                logger.warning(
                    "Insufficient lines for %s cancer-specific analysis "
                    "(mut=%d, wt=%d). Falling back to pan-cancer.",
                    cancer_type, len(mut_avail), len(wt_avail),
                )
                # Pan-cancer fallback: use all available mutant/WT across all lines
                note = (
                    f"Fell back to pan-cancer analysis: only {len(mut_avail)} mutant "
                    f"and {len(wt_avail)} WT lines available for {cancer_type}."
                )
                mut_avail = mutant_ids  # caller must provide pan-cancer lists already
                wt_avail = wt_ids
                scope = AnalysisScope.PAN_CANCER

        return (
            self.crispr.loc[mut_avail],
            self.crispr.loc[wt_avail],
            scope,
            note,
        )

    # ── Statistical tests ──────────────────────────────────────────────────

    @staticmethod
    def _cohend(a: np.ndarray, b: np.ndarray) -> float:
        """Pooled Cohen's d."""
        n1, n2 = len(a), len(b)
        if n1 < 2 or n2 < 2:
            return float("nan")
        pooled_sd = np.sqrt(
            ((n1 - 1) * a.std(ddof=1) ** 2 + (n2 - 1) * b.std(ddof=1) ** 2)
            / (n1 + n2 - 2)
        )
        if pooled_sd == 0:
            return 0.0
        return float((a.mean() - b.mean()) / pooled_sd)

    def _test_gene(
        self, gene: str, mut_dep: pd.Series, wt_dep: pd.Series
    ) -> Dict:
        """
        Run Wilcoxon rank-sum test + compute effect sizes for one partner gene.
        """
        mut_vals = mut_dep.dropna().values
        wt_vals = wt_dep.dropna().values

        if len(mut_vals) < MIN_GROUP_SIZE or len(wt_vals) < MIN_GROUP_SIZE:
            return None

        stat, pval = stats.mannwhitneyu(mut_vals, wt_vals, alternative="less")
        delta = float(mut_vals.mean() - wt_vals.mean())
        d = self._cohend(mut_vals, wt_vals)

        return {
            "gene": gene,
            "delta_dependency": round(delta, 4),
            "p_value": float(pval),
            "n_mut": len(mut_vals),
            "n_wt": len(wt_vals),
            "test_type": "wilcoxon_ranksum_one_sided",
            "effect_size_cohend": round(d, 4) if not np.isnan(d) else None,
        }

    def _compute_codependency(
        self,
        query_gene: str,
        partner_genes: List[str],
        all_lines: List[str],
    ) -> Dict[str, float]:
        """
        Pearson r between the dependency profile of query_gene and each partner
        across all available cell lines.
        """
        if query_gene not in self.crispr.columns:
            return {}
        q_profile = self.crispr.loc[all_lines, query_gene].dropna()
        shared = q_profile.index.tolist()
        codeps: Dict[str, float] = {}
        for pg in partner_genes:
            if pg not in self.crispr.columns:
                continue
            p_profile = self.crispr.loc[shared, pg].dropna()
            common = q_profile.index.intersection(p_profile.index)
            if len(common) < 10:
                continue
            r, _ = stats.pearsonr(q_profile.loc[common], p_profile.loc[common])
            codeps[pg] = round(float(r), 4)
        return codeps

    # ── Main run ───────────────────────────────────────────────────────────

    def compute_sl_partners(
        self,
        query: SLQueryInput,
        mutant_ids: List[str],
        wt_ids: List[str],
    ) -> Tuple[InputContext, List[SLPartner], str]:
        """
        Core SL computation. Returns (input_context, ranked_partners, warnings).
        """
        warnings: List[str] = []
        pan_essential = self._all_pan_essential()

        # Stratify
        mut_dep, wt_dep, scope, note = self._stratify(
            query.gene, mutant_ids, wt_ids, query.cancer_type
        )
        if note:
            warnings.append(note)

        if len(mut_dep) < MIN_GROUP_SIZE:
            raise ValueError(
                f"Not enough mutant lines ({len(mut_dep)}) even in pan-cancer mode. "
                "Consider relaxing mutation_type or checking gene name."
            )

        # Candidate genes: columns present in CRISPR, not the query gene itself
        candidates = [
            g for g in self.crispr.columns
            if g != query.gene and g not in pan_essential
        ]

        # Run per-gene tests
        raw_results = []
        for gene in candidates:
            if gene not in mut_dep.columns:
                continue
            res = self._test_gene(gene, mut_dep[gene], wt_dep[gene])
            if res is not None:
                raw_results.append(res)

        if not raw_results:
            raise ValueError("No candidate genes passed minimum group-size threshold.")

        # Multiple testing correction (BH)
        p_vals = [r["p_value"] for r in raw_results]
        _, fdr_vals, _, _ = multipletests(p_vals, method="fdr_bh")
        for r, fdr in zip(raw_results, fdr_vals):
            r["fdr"] = round(float(fdr), 6)

        # Apply filters
        filtered = [
            r for r in raw_results
            if r["fdr"] <= query.fdr_cutoff
            and abs(r["delta_dependency"]) >= query.delta_dep_cutoff
            and r["delta_dependency"] < 0   # mutant MORE dependent = SL signal
        ]

        # Sort by FDR, then by delta magnitude
        filtered.sort(key=lambda x: (x["fdr"], x["delta_dependency"]))

        # Limit
        top = filtered[: query.top_n_partners]
        top_genes = [r["gene"] for r in top]

        # Co-dependency (optional)
        codeps: Dict[str, float] = {}
        if query.include_codependency and top_genes:
            all_lines = list(set(mut_dep.index.tolist() + wt_dep.index.tolist()))
            codeps = self._compute_codependency(query.gene, top_genes, all_lines)

        # Build SLPartner objects
        partners: List[SLPartner] = []
        for r in top:
            is_essential = r["gene"] in PAN_ESSENTIAL_GENES
            partners.append(
                SLPartner(
                    gene=r["gene"],
                    delta_dependency=r["delta_dependency"],
                    p_value=r["p_value"],
                    fdr=r["fdr"],
                    n_mut=r["n_mut"],
                    n_wt=r["n_wt"],
                    test_type=r["test_type"],
                    effect_size_cohend=r.get("effect_size_cohend"),
                    codependency_r=codeps.get(r["gene"]),
                    pathway=None,  # enriched downstream
                    supporting_frameworks=FrameworkSupport(),
                    is_essential_in_wt=is_essential,
                )
            )

        # Build InputContext
        mut_freq = (
            round(len(mutant_ids) / (len(mutant_ids) + len(wt_ids)), 3)
            if (mutant_ids or wt_ids)
            else None
        )
        ctx = InputContext(
            query_gene=query.gene,
            mutation_type=query.mutation_type,
            cancer_type=query.cancer_type,
            scope=scope,
            n_total_lines=len(mut_dep) + len(wt_dep),
            n_mut_lines=len(mut_dep),
            n_wt_lines=len(wt_dep),
            depmap_release=self.release,
            mutation_frequency_depmap=mut_freq,
            notes=note or None,
        )

        return ctx, partners, "\n".join(warnings)
