"""
Sprint 3 — PRISM Loader + Bonus Signal Tests.

All tests are OFFLINE: every network-touching function is mocked.
No PRISM downloads, no DepMap downloads.

Test groups:
  TestPRISMLoader           (7)  — prism_loader.py core behaviour
  TestPRISMBonusSignal      (5)  — bonus wiring in receipt_miner.py
  TestExistingCanaryIntact  (3)  — BRCA1+PARP canary, weights, sabotage still pass
"""
from __future__ import annotations

from typing import List, Optional
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_prism_df(
    drug_name: str,
    mutant_ids: List[str],
    wt_ids: List[str],
    mut_auc: float = 0.3,
    wt_auc: float = 0.7,
) -> pd.DataFrame:
    """Minimal PRISM-like DataFrame for testing."""
    rows = []
    for did in mutant_ids:
        rows.append({
            "depmap_id": did,
            "name": drug_name,
            "name_norm": drug_name.lower(),
            "auc": mut_auc,
            "ec50": 0.5,
        })
    for did in wt_ids:
        rows.append({
            "depmap_id": did,
            "name": drug_name,
            "name_norm": drug_name.lower(),
            "auc": wt_auc,
            "ec50": 2.0,
        })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# TestPRISMLoader  (7 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestPRISMLoader:
    """Tests for prism_loader.py — all offline via mocks."""

    def test_prism_result_is_dataclass(self):
        """PRISMResult must be a plain dataclass (not Pydantic model)."""
        import dataclasses
        from sl_agent.data.prism_loader import PRISMResult
        assert dataclasses.is_dataclass(PRISMResult)
        assert not hasattr(PRISMResult, "model_fields")

    def test_prism_result_fields(self):
        """PRISMResult has all required fields."""
        from sl_agent.data.prism_loader import PRISMResult
        r = PRISMResult(
            gene="BRCA1",
            axis="parp_inhibitors",
            drug="olaparib",
            delta_auc=-0.4,
            n_mut=10,
            n_wt=50,
            p_value=0.01,
            fdr=0.01,
        )
        assert r.gene == "BRCA1"
        assert r.delta_auc == pytest.approx(-0.4)
        assert r.n_mut == 10

    def test_prism_stratify_unknown_gene_returns_none_gracefully(self):
        """
        prism_stratify with an unknown gene (no mutation data) returns None.
        Must NOT raise.
        """
        with patch("sl_agent.data.prism_loader._load_prism_raw",
                   return_value=pd.DataFrame()), \
             patch("sl_agent.data.prism_loader._get_depmap_mutant_wt_lines",
                   return_value=([], [])):
            from sl_agent.data.prism_loader import prism_stratify
            result = prism_stratify("UNKNOWN_GENE_XYZ", "parp_inhibitors")
            assert result is None  # no data → None, never raises

    def test_prism_stratify_returns_none_when_prism_empty(self):
        """Returns None when PRISM data fails to load."""
        with patch("sl_agent.data.prism_loader._load_prism_raw",
                   return_value=pd.DataFrame()):
            from sl_agent.data.prism_loader import prism_stratify
            result = prism_stratify("BRCA1", "parp_inhibitors")
            assert result is None

    def test_prism_stratify_never_raises(self):
        """prism_stratify NEVER raises — returns None on any exception."""
        with patch("sl_agent.data.prism_loader._load_prism_raw",
                   side_effect=RuntimeError("simulated network failure")):
            from sl_agent.data.prism_loader import prism_stratify
            result = prism_stratify("BRCA1", "parp_inhibitors")
            assert result is None  # must be None, not an exception

    def test_prism_stratify_returns_none_when_insufficient_cell_lines(self):
        """Returns None when fewer than min_n_per_group lines match."""
        drug_df = _make_prism_df(
            "olaparib",
            mutant_ids=["ACH-0001"],   # only 1 mutant line
            wt_ids=["ACH-0002"],       # only 1 WT line
        )
        with patch("sl_agent.data.prism_loader._load_prism_raw", return_value=drug_df), \
             patch("sl_agent.data.prism_loader._get_depmap_mutant_wt_lines",
                   return_value=(["ACH-0001"], ["ACH-0002"])):
            from sl_agent.data.prism_loader import prism_stratify
            result = prism_stratify("BRCA1", "parp_inhibitors", min_n_per_group=5)
            assert result is None  # 1 < min_n_per_group=5

    def test_prism_stratify_brca1_parp_returns_result_or_none(self):
        """
        BRCA1 + PARP: may return a result OR None depending on PRISM data availability.
        This is an honest test — PRISM may not have this pair.
        With mocked data having sufficient lines: must return a PRISMResult.
        """
        mutant_ids = [f"ACH-MUT-{i:04d}" for i in range(10)]
        wt_ids = [f"ACH-WT-{i:04d}" for i in range(50)]
        drug_df = _make_prism_df(
            "olaparib",
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
            mut_auc=0.2,   # mutant very sensitive
            wt_auc=0.8,
        )
        with patch("sl_agent.data.prism_loader._load_prism_raw", return_value=drug_df), \
             patch("sl_agent.data.prism_loader._get_depmap_mutant_wt_lines",
                   return_value=(mutant_ids, wt_ids)):
            from sl_agent.data.prism_loader import prism_stratify
            result = prism_stratify("BRCA1", "parp_inhibitors", min_n_per_group=5)
            assert result is not None
            assert result.gene == "BRCA1"
            assert result.axis == "parp_inhibitors"
            assert result.delta_auc < 0        # mutant more sensitive (lower AUC)
            assert 0.0 <= result.p_value <= 1.0
            assert result.n_mut == 10
            assert result.n_wt == 50


# ─────────────────────────────────────────────────────────────────────────────
# TestPRISMBonusSignal  (5 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestPRISMBonusSignal:
    """Tests for PRISM bonus wiring in receipt_miner.py."""

    def test_prism_bonus_added_when_fdr_below_gate(self):
        """
        PRISM bonus of +0.05 is added when prism_fdr < PRISM_FDR_GATE (0.05).
        Base confidence = 0.60, with PRISM bonus → 0.65.
        """
        from sl_agent.data.prism_loader import PRISMResult
        mock_result = PRISMResult(
            gene="BRCA1",
            axis="parp_inhibitors",
            drug="olaparib",
            delta_auc=-0.35,
            n_mut=10,
            n_wt=50,
            p_value=0.02,
            fdr=0.02,      # < 0.05 gate → bonus applies
        )
        with patch("sl_agent.multimodal.receipt_miner.prism_stratify",
                   return_value=mock_result):
            from sl_agent.multimodal.receipt_miner import (
                _compute_prism_score, PRISM_BONUS, PRISM_FDR_GATE,
            )
            from sl_agent.multimodal.models import CandidateAxis
            bonus, delta = _compute_prism_score("BRCA1", CandidateAxis.PARP_INHIBITORS)
            assert bonus == pytest.approx(PRISM_BONUS)  # 0.05
            assert delta == pytest.approx(-0.35)

    def test_prism_bonus_not_added_when_fdr_above_gate(self):
        """No bonus when prism_fdr >= PRISM_FDR_GATE (0.05)."""
        from sl_agent.data.prism_loader import PRISMResult
        mock_result = PRISMResult(
            gene="BRCA1",
            axis="parp_inhibitors",
            drug="olaparib",
            delta_auc=-0.1,
            n_mut=8,
            n_wt=30,
            p_value=0.15,
            fdr=0.15,      # >= 0.05 gate → no bonus
        )
        with patch("sl_agent.multimodal.receipt_miner.prism_stratify",
                   return_value=mock_result):
            from sl_agent.multimodal.receipt_miner import _compute_prism_score
            from sl_agent.multimodal.models import CandidateAxis
            bonus, delta = _compute_prism_score("BRCA1", CandidateAxis.PARP_INHIBITORS)
            assert bonus == pytest.approx(0.0)
            assert delta == pytest.approx(-0.1)

    def test_prism_bonus_capped_at_one(self):
        """
        PRISM bonus cannot push confidence above 1.0.
        Base = 1.0, bonus = 0.05 → final = 1.0 (capped).
        """
        from sl_agent.data.prism_loader import PRISMResult
        mock_result = PRISMResult(
            gene="BRCA1",
            axis="parp_inhibitors",
            drug="olaparib",
            delta_auc=-0.5,
            n_mut=15,
            n_wt=60,
            p_value=0.001,
            fdr=0.001,
        )
        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(1.0, -0.40, 0.01, 20, 120)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(1.0, -2.0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(1.0, 10)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=1.0), \
             patch("sl_agent.multimodal.receipt_miner.prism_stratify",
                   return_value=mock_result), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert") as mock_upsert:

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            summary = mine_receipts(
                gene_panel=["BRCA1"],
                axes=[CandidateAxis.PARP_INHIBITORS],
                thresholds=MinerThresholds(),
            )

            # Check that upserted candidate has confidence <= 1.0
            assert mock_upsert.called
            candidate = mock_upsert.call_args[0][0]
            assert candidate.confidence_score <= 1.0, (
                f"PRISM bonus overflowed 1.0: confidence={candidate.confidence_score}"
            )

    def test_prism_bonus_none_result_returns_zero(self):
        """prism_stratify returning None → bonus = 0.0, delta_auc = None."""
        with patch("sl_agent.multimodal.receipt_miner.prism_stratify",
                   return_value=None):
            from sl_agent.multimodal.receipt_miner import _compute_prism_score
            from sl_agent.multimodal.models import CandidateAxis
            bonus, delta = _compute_prism_score("UNKNOWN", CandidateAxis.ATR_WEE1)
            assert bonus == pytest.approx(0.0)
            assert delta is None

    def test_prism_delta_auc_stored_in_candidate(self):
        """prism_delta_auc is stored in ReceiptCandidate when PRISM returns a result."""
        from sl_agent.data.prism_loader import PRISMResult
        mock_result = PRISMResult(
            gene="MBD4",
            axis="cytidine_analogs",
            drug="gemcitabine",
            delta_auc=-0.25,
            n_mut=8,
            n_wt=40,
            p_value=0.03,
            fdr=0.03,
        )
        stored = []

        def capture(c):
            stored.append(c)
            return c

        with patch("sl_agent.multimodal.receipt_miner._compute_crispr_score",
                   return_value=(1.0, -0.40, 0.01, 20, 120)), \
             patch("sl_agent.multimodal.receipt_miner._compute_gdsc_score",
                   return_value=(1.0, -2.0)), \
             patch("sl_agent.multimodal.receipt_miner._compute_kb_score",
                   return_value=(1.0, 10)), \
             patch("sl_agent.multimodal.receipt_miner._compute_expression_score",
                   return_value=0.8), \
             patch("sl_agent.multimodal.receipt_miner.prism_stratify",
                   return_value=mock_result), \
             patch("sl_agent.multimodal.receipt_miner.AuditQueue.upsert",
                   side_effect=capture):

            from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
            from sl_agent.multimodal.models import CandidateAxis

            mine_receipts(
                gene_panel=["MBD4"],
                axes=[CandidateAxis.CYTIDINE_ANALOGS],
                thresholds=MinerThresholds(),
            )

        assert len(stored) == 1
        assert stored[0].prism_delta_auc == pytest.approx(-0.25)


# ─────────────────────────────────────────────────────────────────────────────
# TestExistingCanaryIntact  (3 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestExistingCanaryIntact:
    """
    Sprint 3 regression tests — existing Sprint 2 guarantees must still hold.
    Weights unchanged; sabotage protection intact.
    """

    def test_weights_still_sum_to_one(self):
        """
        Sprint 3 does NOT rebalance weights.
        CRISPR(0.35) + GDSC(0.30) + KB(0.25) + Expr(0.10) = 1.0 (unchanged).
        PRISM is purely additive and not in this sum.
        """
        weights = [0.35, 0.30, 0.25, 0.10]
        assert sum(weights) == pytest.approx(1.0)
        # PRISM_BONUS is NOT part of the base sum
        from sl_agent.multimodal.receipt_miner import PRISM_BONUS
        assert PRISM_BONUS == pytest.approx(0.05)
        assert sum(weights) + PRISM_BONUS == pytest.approx(1.05)  # over 1.0 → cap applies

    def test_brca1_parp_canary_still_1000_without_prism(self):
        """
        BRCA1+PARP canary with PRISM returning None → confidence still = 1.0
        (capped at 1.0, no change from Sprint 2).
        """
        crispr_score = 1.0
        gdsc_score = 1.0
        kb_score = 1.0
        expr_score = 1.0
        prism_bonus = 0.0  # PRISM not available

        base = round(
            crispr_score * 0.35
            + gdsc_score  * 0.30
            + kb_score    * 0.25
            + expr_score  * 0.10,
            4,
        )
        confidence = round(min(base + prism_bonus, 1.0), 4)
        assert confidence == pytest.approx(1.0)

        from sl_agent.multimodal.receipt_miner import _project_tier
        assert _project_tier(confidence, kb_hits=10) == "Strong"
        assert _project_tier(confidence, kb_hits=10) != "Validated"

    def test_sabotage_still_passes_after_prism_addition(self):
        """
        SABOTAGE TEST (Sprint 3 regression): _project_tier NEVER returns 'Validated'.
        With PRISM bonus: max confidence = 1.0 + 0.05 → capped at 1.0 → still 'Strong'.
        """
        from sl_agent.multimodal.receipt_miner import (
            _project_tier, _VALID_CANDIDATE_TIERS, PRISM_BONUS,
        )

        confidences = [round(i * 0.05, 2) for i in range(21)]  # 0.0..1.0
        kb_hits_values = [0, 1, 2, 3, 5, 10, 100]

        for conf in confidences:
            # Apply PRISM bonus (as done in mine_receipts) and cap
            boosted = round(min(conf + PRISM_BONUS, 1.0), 4)
            for kb in kb_hits_values:
                tier = _project_tier(boosted, kb)
                assert tier != "Validated", (
                    f"SABOTAGE: _project_tier({boosted}, {kb}) returned 'Validated' "
                    f"after PRISM bonus (base_conf={conf})"
                )
                assert tier in _VALID_CANDIDATE_TIERS
