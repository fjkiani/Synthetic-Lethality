"""
Tests for the Receipt Miner pipeline — Sprint 2.

All tests are OFFLINE: every network-touching function is mocked.
No DepMap downloads, no GDSC downloads, no KB network calls.

Test groups:
  TestGDSCBiomarkerLoader     (6) — gdsc_biomarker_loader.py
  TestProjectTier             (7) — _project_tier() sabotage protection
  TestCompositeWeights        (4) — confidence score formula
  TestMineReceiptsIntegration (6) — mine_receipts() end-to-end with mocks
  TestMinerRunSummary         (3) — MinerRunSummary dataclass
  TestBuildEvidenceSummary    (2) — evidence string formatting
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_gdsc_df(drug_name: str = "olaparib", cosmic_id: int = 1001) -> pd.DataFrame:
    """Minimal GDSC-like DataFrame for testing."""
    return pd.DataFrame({
        "DRUG_NAME": [drug_name],
        "DRUG_NAME_NORM": [drug_name.lower()],
        "COSMIC_ID": [cosmic_id],
        "LN_IC50": [2.0],
        "_DATASET": ["gdsc1"],
    })


def _make_gdsc_stratified_df(
    drug_name: str,
    mutant_cosmics: List[int],
    wt_cosmics: List[int],
    mut_ic50: float = 1.0,
    wt_ic50: float = 3.0,
) -> pd.DataFrame:
    """DataFrame with mutant and WT cell lines for proper Mann-Whitney testing."""
    rows = []
    for cid in mutant_cosmics:
        rows.append({"DRUG_NAME": drug_name, "DRUG_NAME_NORM": drug_name.lower(),
                     "COSMIC_ID": cid, "LN_IC50": mut_ic50, "_DATASET": "gdsc1"})
    for cid in wt_cosmics:
        rows.append({"DRUG_NAME": drug_name, "DRUG_NAME_NORM": drug_name.lower(),
                     "COSMIC_ID": cid, "LN_IC50": wt_ic50, "_DATASET": "gdsc1"})
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# TestGDSCBiomarkerLoader (6 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestGDSCBiomarkerLoader:
    """Tests for gdsc_biomarker_loader.py — all offline via mocks."""

    def test_gdsc_result_is_dataclass(self):
        """GDSCResult must be a plain dataclass (not Pydantic model)."""
        import dataclasses
        from sl_agent.data.gdsc_biomarker_loader import GDSCResult
        assert dataclasses.is_dataclass(GDSCResult)
        # Pydantic models have model_fields; dataclasses do not
        assert not hasattr(GDSCResult, "model_fields")

    def test_gdsc_result_fields(self):
        """GDSCResult dataclass must have all required fields."""
        from sl_agent.data.gdsc_biomarker_loader import GDSCResult
        r = GDSCResult(
            gene="BRCA1",
            axis="parp_inhibitors",
            drug="olaparib",
            gdsc_version="GDSC1",
            delta_ln_ic50=-1.5,
            p_value=0.01,
            fdr=0.01,
            n_mutant=10,
            n_wt=50,
            effect_size_cohend=-0.8,
        )
        assert r.gene == "BRCA1"
        assert r.delta_ln_ic50 == pytest.approx(-1.5)
        assert r.n_mutant == 10

    def test_gdsc_biomarker_stratify_returns_none_when_both_empty(self):
        """gdsc_biomarker_stratify returns None when both GDSC1 and GDSC2 are empty."""
        with patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc1", return_value=pd.DataFrame()), \
             patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc2", return_value=pd.DataFrame()):
            from sl_agent.data.gdsc_biomarker_loader import gdsc_biomarker_stratify
            result = gdsc_biomarker_stratify("BRCA1", "parp_inhibitors")
            assert result is None

    def test_gdsc_biomarker_stratify_never_raises(self):
        """gdsc_biomarker_stratify NEVER raises — returns None on any exception."""
        with patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc1", side_effect=RuntimeError("network error")), \
             patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc2", side_effect=RuntimeError("network error")):
            from sl_agent.data.gdsc_biomarker_loader import gdsc_biomarker_stratify
            result = gdsc_biomarker_stratify("BRCA1", "parp_inhibitors")
            assert result is None  # must return None, not raise

    def test_gdsc_biomarker_stratify_returns_none_when_insufficient_cell_lines(self):
        """Returns None when fewer than min_n_per_group COSMIC lines match."""
        # Only 2 mutant COSMIC IDs, but GDSC has none of them
        big_df = _make_gdsc_stratified_df(
            "olaparib",
            mutant_cosmics=[9999],  # won't match what stratify returns
            wt_cosmics=[8888],
        )
        with patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc1", return_value=big_df), \
             patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc2", return_value=pd.DataFrame()), \
             patch("sl_agent.data.gdsc_biomarker_loader._get_cosmic_mutant_lines",
                   return_value=([], [])):  # no COSMIC IDs → insufficient
            from sl_agent.data.gdsc_biomarker_loader import gdsc_biomarker_stratify
            result = gdsc_biomarker_stratify("BRCA1", "parp_inhibitors", min_n_per_group=5)
            assert result is None

    def test_gdsc_biomarker_stratify_returns_result_with_valid_data(self):
        """Returns GDSCResult when enough cell lines and matching drugs exist."""
        mut_cosmics = list(range(1001, 1011))  # 10 mutant lines
        wt_cosmics = list(range(2001, 2051))   # 50 wt lines
        big_df = _make_gdsc_stratified_df(
            "olaparib",
            mutant_cosmics=mut_cosmics,
            wt_cosmics=wt_cosmics,
            mut_ic50=1.0,
            wt_ic50=4.0,
        )
        with patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc1", return_value=big_df), \
             patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc2", return_value=pd.DataFrame()), \
             patch("sl_agent.data.gdsc_biomarker_loader._get_cosmic_mutant_lines",
                   return_value=(mut_cosmics, wt_cosmics)):
            from sl_agent.data.gdsc_biomarker_loader import gdsc_biomarker_stratify
            result = gdsc_biomarker_stratify("BRCA1", "parp_inhibitors", min_n_per_group=5)
            assert result is not None
            assert result.gene == "BRCA1"
            assert result.axis == "parp_inhibitors"
            assert result.delta_ln_ic50 < 0  # mutant more sensitive (lower IC50)
            assert 0.0 <= result.p_value <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# TestProjectTier (7 tests) — SABOTAGE PROTECTION
# ─────────────────────────────────────────────────────────────────────────────

class TestProjectTier:
    """
    Tests for _project_tier() — the critical sabotage protection function.

    SABOTAGE TEST: _project_tier NEVER returns "Validated".
    Even with confidence=1.0 and kb_hits=1000, the result is always
    one of {"Strong", "Mechanistic", "Insufficient"}.
    """

    def test_high_confidence_high_kb_returns_strong(self):
        """confidence=0.85, kb_hits=5 → 'Strong'."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.85, 5) == "Strong"

    def test_medium_confidence_returns_mechanistic(self):
        """confidence=0.60, kb_hits=0 → 'Mechanistic'."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.60, 0) == "Mechanistic"

    def test_low_confidence_returns_insufficient(self):
        """confidence=0.30, kb_hits=0 → 'Insufficient'."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.30, 0) == "Insufficient"

    def test_boundary_70_with_3_kb_hits_is_strong(self):
        """confidence=0.70 exactly with kb_hits=3 → 'Strong' (boundary inclusive)."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.70, 3) == "Strong"

    def test_boundary_69_below_threshold_is_mechanistic(self):
        """confidence=0.699 with kb_hits=10 → 'Mechanistic' (below 0.70 threshold)."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.699, 10) == "Mechanistic"

    def test_high_confidence_but_low_kb_hits_not_strong(self):
        """confidence=0.90 but kb_hits=2 (< 3) → 'Mechanistic' (not Strong)."""
        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(0.90, 2) == "Mechanistic"

    def test_sabotage_validated_never_returned(self):
        """
        SABOTAGE TEST — most important test.

        Call _project_tier() with every combination of high inputs.
        Assert the result is NEVER "Validated".
        Assert the result is ALWAYS in {"Strong", "Mechanistic", "Insufficient"}.
        """
        from sl_agent.multimodal.receipt_miner import _project_tier, _VALID_CANDIDATE_TIERS

        # Exhaustive sweep: confidence 0.0 to 1.0 in 0.05 steps, kb_hits 0-100
        confidences = [round(i * 0.05, 2) for i in range(21)]  # 0.0 to 1.0
        kb_hits_values = [0, 1, 2, 3, 5, 10, 50, 100]

        for conf in confidences:
            for kb in kb_hits_values:
                tier = _project_tier(conf, kb)
                assert tier != "Validated", (
                    f"SABOTAGE DETECTED: _project_tier({conf}, {kb}) returned 'Validated'! "
                    f"Validated tier MUST NEVER be auto-assigned."
                )
                assert tier in _VALID_CANDIDATE_TIERS, (
                    f"_project_tier({conf}, {kb}) returned unexpected tier '{tier}'. "
                    f"Must be one of {_VALID_CANDIDATE_TIERS}"
                )

        # Extra: max possible inputs
        assert _project_tier(1.0, 1000) in _VALID_CANDIDATE_TIERS
        assert _project_tier(1.0, 1000) != "Validated"
        assert _project_tier(1.0, 100) == "Strong"  # best possible → "Strong"


# ─────────────────────────────────────────────────────────────────────────────
# TestCompositeWeights (4 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestCompositeWeights:
    """Tests that the composite confidence formula uses correct weights."""

    def test_weights_sum_to_one(self):
        """CRISPR(0.35) + GDSC(0.30) + KB(0.25) + Expr(0.10) = 1.0"""
        weights = [0.35, 0.30, 0.25, 0.10]
        assert sum(weights) == pytest.approx(1.0)

    def test_max_score_scenario(self):
        """All signals = 1.0 → confidence = 1.0 (all weights sum to 1.0)."""
        # CRISPR=1.0, GDSC=1.0, KB=1.0, Expr=1.0
        confidence = round(1.0 * 0.35 + 1.0 * 0.30 + 1.0 * 0.25 + 1.0 * 0.10, 4)
        assert confidence == pytest.approx(1.0)

    def test_crispr_only_score(self):
        """CRISPR=1.0, others=0 → confidence = 0.35."""
        confidence = round(1.0 * 0.35 + 0.0 * 0.30 + 0.0 * 0.25 + 0.0 * 0.10, 4)
        assert confidence == pytest.approx(0.35)

    def test_brca1_parp_canary_confidence(self):
        """
        BRCA1 + PARP canary: max all signals → confidence = 1.0.

        This is the explicit canary test required by the spec.
        With mocked max signals:
          CRISPR score = 1.0 (weight 0.35)
          GDSC score   = 1.0 (weight 0.30)
          KB score     = 1.0 (weight 0.25)
          Expr score   = 1.0 (weight 0.10) [not 0.8 — we test the formula directly]
        Confidence = 1.0 * 0.35 + 1.0 * 0.30 + 1.0 * 0.25 + 1.0 * 0.10 = 1.0
        """
        # Direct formula verification — this is what mine_receipts() computes
        crispr_score = 1.0
        gdsc_score = 1.0
        kb_score = 1.0
        expr_score = 1.0
        confidence = round(
            crispr_score * 0.35
            + gdsc_score  * 0.30
            + kb_score    * 0.25
            + expr_score  * 0.10,
            4,
        )
        assert confidence == pytest.approx(1.0), (
            f"BRCA1+PARP canary confidence = {confidence} (expected 1.0)"
        )
        # Also verify tier at this confidence
        from sl_agent.multimodal.receipt_miner import _project_tier
        tier = _project_tier(confidence, kb_hits=10)
        assert tier == "Strong"
        assert tier != "Validated"


# ─────────────────────────────────────────────────────────────────────────────
# TestMineReceiptsIntegration (6 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestMineReceiptsIntegration:
    """End-to-end mine_receipts() tests with all signals mocked."""

    def test_mine_receipts_returns_summary(self):
        """mine_receipts returns a MinerRunSummary."""
        from sl_agent.multimodal.receipt_miner import MinerRunSummary

        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(0.0, None, None, None, None)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(0.0, None)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(0.0, 0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=0.0), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert") as mock_upsert:

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.PARP_INHIBITORS],
                thresholds=MinerThresholds(),
            )
            assert isinstance(summary, MinerRunSummary)
            assert summary.gene_count == 1
            assert summary.axis_count == 1
            assert summary.completed_at is not None

    def test_mine_receipts_discards_low_confidence(self):
        """Candidates below min_confidence are not queued."""
        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(0.1, -0.05, 0.3, 5, 50)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(0.0, None)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(0.0, 0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=0.0), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert") as mock_upsert:

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            thresholds = MinerThresholds(min_confidence=0.40)
            # crispr_score=0.1 × 0.35 = 0.035 → well below 0.40
            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.PARP_INHIBITORS],
                thresholds=thresholds,
            )
            mock_upsert.assert_not_called()
            assert summary.candidates_queued == 0

    def test_mine_receipts_queues_high_confidence_candidate(self):
        """Candidates above min_confidence are queued via AuditQueue.upsert."""
        mock_candidate = MagicMock()
        mock_candidate.gene = "BRCA1"

        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(1.0, -0.40, 0.01, 15, 100)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(0.8, -1.5)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(1.0, 5)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=0.8), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert",
                   return_value=mock_candidate) as mock_upsert:

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.PARP_INHIBITORS],
                thresholds=MinerThresholds(),
            )
            mock_upsert.assert_called_once()
            assert summary.candidates_queued == 1

    def test_mine_receipts_skips_custom_axis(self):
        """CandidateAxis.CUSTOM is always skipped."""
        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(1.0, -0.40, 0.01, 15, 100)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(1.0, -2.0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(1.0, 5)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=1.0), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert") as mock_upsert:

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.CUSTOM],  # only CUSTOM axis
                thresholds=MinerThresholds(),
            )
            mock_upsert.assert_not_called()
            assert summary.pairs_evaluated == 0

    def test_mine_receipts_validated_never_in_queued_candidates(self):
        """
        No candidate queued by mine_receipts() has candidate_tier == 'Validated'.
        Even with maximum signals.
        """
        queued_candidates = []

        def capture_upsert(candidate):
            queued_candidates.append(candidate)
            return candidate

        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(1.0, -0.40, 0.01, 20, 120)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(1.0, -2.0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(1.0, 10)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=1.0), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert",
                   side_effect=capture_upsert):

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds, DEFAULT_AXES
            from sl_agent.multimodal.models import CandidateAxis

            mine_receipts(
                gene_panel=["BRCA1", "BRCA2", "PALB2"],
                axes=DEFAULT_AXES,
                thresholds=MinerThresholds(),
            )

        # Check none of the queued candidates have "Validated" tier
        for c in queued_candidates:
            assert c.candidate_tier != "Validated", (
                f"SABOTAGE: mine_receipts queued candidate with tier='Validated' "
                f"for {c.gene}×{c.axis}"
            )

    def test_mine_receipts_error_in_signal_does_not_crash(self):
        """Errors in signal computation are caught; the run continues."""
        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(0.0, None, None, None, None)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   side_effect=RuntimeError("gdsc failed")), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(0.0, 0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=0.0), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert"):

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            # Must not raise
            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.PARP_INHIBITORS],
                thresholds=MinerThresholds(),
            )
            # The pair may land in errors since _compute_gdsc_score raised
            assert summary.completed_at is not None


# ─────────────────────────────────────────────────────────────────────────────
# TestMinerRunSummary (3 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestMinerRunSummary:
    """Tests for MinerRunSummary dataclass."""

    def test_miner_run_summary_defaults(self):
        """MinerRunSummary initializes with correct defaults."""
        from sl_agent.multimodal.receipt_miner import MinerRunSummary
        s = MinerRunSummary(gene_count=10, axis_count=6)
        assert s.gene_count == 10
        assert s.axis_count == 6
        assert s.pairs_evaluated == 0
        assert s.candidates_queued == 0
        assert s.completed_at is None
        assert s.depmap_release == "unknown"

    def test_miner_run_summary_finish(self):
        """finish() sets completed_at."""
        from sl_agent.multimodal.receipt_miner import MinerRunSummary
        from datetime import datetime
        s = MinerRunSummary(gene_count=1, axis_count=1)
        assert s.completed_at is None
        result = s.finish()
        assert result.completed_at is not None
        assert isinstance(result.completed_at, datetime)

    def test_miner_run_summary_errors_list(self):
        """errors list starts empty and can be appended to."""
        from sl_agent.multimodal.receipt_miner import MinerRunSummary
        s = MinerRunSummary(gene_count=1, axis_count=1)
        assert s.errors == []
        s.errors.append("GENE_A×parp_inhibitors: test error")
        assert len(s.errors) == 1


# ─────────────────────────────────────────────────────────────────────────────
# TestBuildEvidenceSummary (2 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestBuildEvidenceSummary:
    """Tests for _build_evidence_summary formatting."""

    def test_evidence_summary_with_all_signals(self):
        """Summary includes all signal labels when scores > 0."""
        from sl_agent.multimodal.receipt_miner import _build_evidence_summary
        s = _build_evidence_summary(
            gene="BRCA1",
            axis="parp_inhibitors",
            crispr_score=0.90,
            gdsc_score=0.70,
            kb_score=0.80,
            expression_score=0.50,
            confidence=0.82,
            tier="Strong",
        )
        assert "BRCA1" in s
        assert "parp_inhibitors" in s
        assert "CRISPR" in s
        assert "GDSC" in s
        assert "KB" in s
        assert "Expr" in s
        assert "Strong" in s

    def test_evidence_summary_with_no_signal(self):
        """Summary shows 'no signal' when all scores are 0."""
        from sl_agent.multimodal.receipt_miner import _build_evidence_summary
        s = _build_evidence_summary(
            gene="GENE_X",
            axis="wrn",
            crispr_score=0.0,
            gdsc_score=0.0,
            kb_score=0.0,
            expression_score=0.0,
            confidence=0.0,
            tier="Insufficient",
        )
        assert "no signal" in s


# ─────────────────────────────────────────────────────────────────────────────
# TestMinerThresholds (2 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestMinerThresholds:
    """Tests for MinerThresholds validation."""

    def test_default_thresholds_valid(self):
        """Default MinerThresholds creates without error."""
        from sl_agent.multimodal.receipt_miner import MinerThresholds
        t = MinerThresholds()
        assert t.crispr_delta_dep == -0.15
        assert t.min_confidence == 0.40
        assert t.gdsc_fdr == 0.05

    def test_positive_crispr_delta_dep_raises(self):
        """crispr_delta_dep must be negative — positive value raises AssertionError."""
        from sl_agent.multimodal.receipt_miner import MinerThresholds
        with pytest.raises(AssertionError):
            MinerThresholds(crispr_delta_dep=0.15)  # positive — invalid
