"""Tests for BrM matrix builder (3-section partitioned output)."""

import datetime
import pytest

from sl_agent.brm import build_brm_matrix, BrMQueryInput
from sl_agent.brm.models import BrMRowClass, ExploitMode, CascadeStep


class TestMatrixStructure:
    def setup_method(self):
        self.matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)

    def test_returns_seven_rows(self):
        assert len(self.matrix.rows) == 7

    def test_bace1_is_first_row(self):
        assert self.matrix.rows[0].gene == "BACE1"

    def test_exploit_mode_summary_sums_to_row_count(self):
        total = sum(self.matrix.exploit_mode_summary.values())
        assert total == len(self.matrix.rows)

    def test_row_class_summary_sums_to_row_count(self):
        total = sum(self.matrix.row_class_summary.values())
        assert total == len(self.matrix.rows)

    def test_frozen_at_parseable_iso8601(self):
        ts = self.matrix.frozen_at
        # Should parse without error
        parsed = datetime.datetime.fromisoformat(ts.replace("+00:00", "+00:00"))
        assert parsed is not None

    def test_ruo_true(self):
        assert self.matrix.ruo is True

    def test_disclaimer_non_empty(self):
        assert len(self.matrix.disclaimer) > 0

    def test_version_v2(self):
        assert self.matrix.version == "v2"

    def test_depmap_release_24q4(self):
        assert self.matrix.depmap_release == "24Q4"


class TestMatrixSections:
    def setup_method(self):
        self.matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)

    def test_calibration_rows_non_empty(self):
        assert len(self.matrix.calibration_rows) >= 1

    def test_bace1_in_calibration(self):
        cal_genes = [r.gene for r in self.matrix.calibration_rows]
        assert "BACE1" in cal_genes

    def test_calibration_before_novel_in_rows(self):
        """All calibration rows appear before novel rows in the ordered list."""
        row_classes = [r.row_class for r in self.matrix.rows]
        cal_indices  = [i for i, c in enumerate(row_classes) if c == BrMRowClass.CALIBRATION]
        novel_indices = [i for i, c in enumerate(row_classes) if c == BrMRowClass.NOVEL]
        if cal_indices and novel_indices:
            assert max(cal_indices) < min(novel_indices)

    def test_row_class_summary_keys(self):
        assert "calibration" in self.matrix.row_class_summary
        assert "novel" in self.matrix.row_class_summary
        assert "negative" in self.matrix.row_class_summary


class TestMatrixNonBlocking:
    def test_depmap_note_none_when_no_data(self):
        matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)
        for row in matrix.rows:
            assert row.depmap_expression_note is None

    def test_sl_partners_empty_when_no_data(self):
        matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)
        for row in matrix.rows:
            assert row.sl_partners == []

    def test_novelty_score_none_when_no_data(self):
        matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)
        # CALIBRATION rows should have None novelty_score when no DepMap data
        for row in matrix.calibration_rows:
            assert row.novelty_score is None


class TestMatrixInvariants:
    def setup_method(self):
        self.matrix = build_brm_matrix(BrMQueryInput(), expression_df=None, crispr_df=None, sample_info=None)

    def test_all_confidence_scores_in_range(self):
        for row in self.matrix.rows:
            assert 0.0 <= row.confidence_score <= 1.0, (
                f"{row.gene} confidence_score={row.confidence_score} out of range"
            )

    def test_evidence_tier_never_validated(self):
        for row in self.matrix.rows:
            assert row.evidence_tier != "Validated", f"{row.gene} has 'Validated' tier"

    def test_cascade_step_summary_valid_values(self):
        valid_steps = {e.value for e in CascadeStep}
        for step in self.matrix.cascade_step_summary:
            assert step in valid_steps, f"Invalid cascade step: {step}"

    def test_exploit_mode_summary_valid_values(self):
        valid_modes = {e.value for e in ExploitMode}
        for mode in self.matrix.exploit_mode_summary:
            assert mode in valid_modes, f"Invalid exploit mode: {mode}"

    def test_secondary_never_equals_primary(self):
        for row in self.matrix.rows:
            for mode in row.secondary_exploit_modes:
                assert mode != row.primary_exploit_mode, (
                    f"{row.gene}: secondary {mode} == primary {row.primary_exploit_mode}"
                )

    def test_matrix_has_seven_rows(self):
        assert len(self.matrix.rows) == 7
