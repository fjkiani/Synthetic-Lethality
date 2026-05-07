"""Tests for expression-stratified SL discovery."""

import numpy as np
import pandas as pd
import pytest

from sl_agent.brm.expression_stratifier import (
    run_expression_stratified_sl,
    run_batch_stratification,
    FDR_MODERATE,
    FDR_STRONG,
    MIN_LINES,
)
from sl_agent.brm.models import SLPartner


def _make_synthetic_data(n_lines: int = 50, seed: int = 42):
    """
    Create synthetic expression + CRISPR data with a planted SL signal.

    ZEB1 high-expression lines have lower ITGAV dependency (delta ~ -0.5).
    """
    rng = np.random.default_rng(seed)
    genes = ["ZEB1", "VIM", "SPARC", "TGFB1", "ITGB1"]
    crispr_genes = ["ITGAV", "FERMT2", "CDS2", "LMNA", "PTK2", "RANDOM1", "RANDOM2"]
    line_ids = [f"ACH-{i:06d}" for i in range(n_lines)]

    # Expression: ZEB1 has bimodal distribution
    expr_data = rng.normal(5, 1, size=(n_lines, len(genes)))
    expr_df = pd.DataFrame(expr_data, index=line_ids, columns=genes)
    # Make top 25% of ZEB1 clearly high
    zeb1_vals = np.sort(rng.normal(5, 1, n_lines))
    zeb1_vals[-int(n_lines * 0.25):] += 3  # high group
    expr_df["ZEB1"] = zeb1_vals

    # CRISPR: ITGAV is more essential in ZEB1-high lines
    crispr_data = rng.normal(-0.1, 0.2, size=(n_lines, len(crispr_genes)))
    crispr_df = pd.DataFrame(crispr_data, index=line_ids, columns=crispr_genes)
    # Plant signal: ZEB1-high lines have lower ITGAV dependency
    high_mask = expr_df["ZEB1"] >= expr_df["ZEB1"].quantile(0.75)
    crispr_df.loc[high_mask, "ITGAV"] -= 0.5

    # Sample info
    sample_info = pd.DataFrame(
        {"OncotreeLineage": ["Non-Small Cell Lung Cancer"] * n_lines},
        index=line_ids,
    )

    return expr_df, crispr_df, sample_info


class TestRunExpressionStratifiedSL:
    def setup_method(self):
        self.expr_df, self.crispr_df, self.sample_info = _make_synthetic_data(n_lines=60)

    def test_returns_stratification_result(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        assert result.gene == "ZEB1"
        assert result.error is None

    def test_detects_planted_sl_signal(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        partner_genes = [p.partner_gene for p in result.partners]
        assert "ITGAV" in partner_genes, "Planted ITGAV signal not detected"

    def test_partners_have_negative_delta(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        for p in result.partners:
            assert p.delta_dep < 0, f"Partner {p.partner_gene} has positive delta"

    def test_partners_fdr_within_threshold(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        for p in result.partners:
            assert p.fdr <= FDR_MODERATE

    def test_missing_gene_returns_error(self):
        result = run_expression_stratified_sl(
            "NOTREAL", self.expr_df, self.crispr_df, self.sample_info
        )
        assert result.error is not None

    def test_novelty_score_in_range(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        assert 0.0 <= result.novelty_score <= 1.0

    def test_v2_v1_v3_scores_computed(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        assert result.novelty_score >= 0.0
        assert result.novelty_score_v1 >= 0.0
        assert result.novelty_score_v3 >= 0.0

    def test_partner_quality_tiers_valid(self):
        result = run_expression_stratified_sl(
            "ZEB1", self.expr_df, self.crispr_df, self.sample_info
        )
        for p in result.partners:
            assert p.partner_quality_tier in ("strong", "moderate", "weak")


class TestBatchStratification:
    def test_returns_dict_keyed_by_gene(self):
        expr_df, crispr_df, sample_info = _make_synthetic_data(n_lines=60)
        results = run_batch_stratification(
            ["ZEB1", "VIM"], expr_df, crispr_df, sample_info
        )
        assert "ZEB1" in results
        assert "VIM" in results

    def test_missing_gene_in_batch_has_error(self):
        expr_df, crispr_df, sample_info = _make_synthetic_data(n_lines=60)
        results = run_batch_stratification(
            ["NOTREAL"], expr_df, crispr_df, sample_info
        )
        assert results["NOTREAL"].error is not None

    def test_insufficient_lines_returns_error(self):
        expr_df, crispr_df, sample_info = _make_synthetic_data(n_lines=8)
        result = run_expression_stratified_sl(
            "ZEB1", expr_df, crispr_df, sample_info
        )
        # With only 8 lines, quartile split gives < MIN_LINES per group
        assert result.error is not None
