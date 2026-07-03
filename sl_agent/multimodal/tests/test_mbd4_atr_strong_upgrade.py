"""
Integration test: MBD4/ATR_WEE1 STRONG tier upgrade via canonical GDSC2 receipt.

Verifies that when matrix_builder runs with include_literature_receipts=True
and no live GDSC data, the frozen GDSC receipt (ceralasertib, p=0.021, d=-0.503)
flows into EvidenceRow.gdsc and triggers "Strong candidate dependency axis" tier.

Canonical evidence source:
  fjkiani/crispro/publications/00-mbd4-manuscript/mbd4_parp_response/
  artifacts/canonical_atr_wee1_rerun_20260405/
  n_LOF=14 vs n_WT=914, Δ LN_IC50=−0.732, p=0.021, Cohen's d=−0.503
"""
import pytest
from sl_agent.multimodal.models import CandidateAxis, ModalityStatus
from sl_agent.multimodal.literature_receipts import get_literature_receipts
from sl_agent.multimodal.modality_fuser import _assign_recommendation_tier, _weighted_score


class TestMBD4ATRStrongUpgrade:
    """Verify MBD4/ATR_WEE1 reaches Strong tier via canonical GDSC2 receipt."""

    def test_mbd4_atr_gdsc_receipt_is_positive(self):
        """GDSC receipt for MBD4/ATR_WEE1 must be POSITIVE with correct stats."""
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)
        assert "gdsc" in receipts, "GDSC receipt must exist for MBD4/ATR_WEE1"
        gdsc = receipts["gdsc"]
        assert gdsc.status == ModalityStatus.POSITIVE, (
            f"GDSC status must be POSITIVE, got {gdsc.status}"
        )
        assert gdsc.p_value == pytest.approx(0.021, abs=1e-4), (
            f"p_value must be 0.021 (ceralasertib primary), got {gdsc.p_value}"
        )
        assert gdsc.effect_size == pytest.approx(-0.503, abs=1e-4), (
            f"Cohen's d must be -0.503, got {gdsc.effect_size}"
        )
        assert gdsc.n_mut == 14, f"n_LOF must be 14, got {gdsc.n_mut}"
        assert gdsc.n_wt == 914, f"n_WT must be 914, got {gdsc.n_wt}"
        assert gdsc.delta_auc == pytest.approx(-0.056, abs=1e-4), (
            f"delta_AUC must be -0.056, got {gdsc.delta_auc}"
        )

    def test_mbd4_atr_expression_receipt_is_positive(self):
        """Expression receipt must be POSITIVE (mechanistic chain)."""
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)
        assert "expression" in receipts
        assert receipts["expression"].status == ModalityStatus.POSITIVE

    def test_mbd4_atr_tier_logic_fires_strong(self):
        """
        With gdsc=POSITIVE and expression=POSITIVE, tier logic must return
        'Strong candidate dependency axis'.

        Tier rule: (crispr_pos OR pharma_pos) AND n_positive >= 2 AND NOT crispr_neg
        pharma_pos = gdsc.POSITIVE → True
        n_positive = 2 (gdsc + expression)
        crispr_neg = False (crispr is MISSING)
        → "Strong candidate dependency axis"
        """
        from sl_agent.multimodal.models import (
            EvidenceRow, ModalityEvidence, Modality, ModalityStatus
        )
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)

        # Build EvidenceRow with receipts applied (simulating matrix_builder Step 4)
        row = EvidenceRow(
            axis=CandidateAxis.ATR_WEE1,
            axis_label="ATR/WEE1 Inhibitors",
            mechanism="BER stress → ssDNA at stalled forks → ATR/WEE1 checkpoint dependency",
        )
        # Apply receipts (as matrix_builder Step 4 does)
        if "gdsc" in receipts and row.gdsc.status == ModalityStatus.MISSING:
            row.gdsc = receipts["gdsc"]
        if "expression" in receipts and row.expression.status == ModalityStatus.MISSING:
            row.expression = receipts["expression"]
        if "in_vitro" in receipts:
            row.in_vitro = receipts["in_vitro"]
        if "in_vivo" in receipts:
            row.in_vivo = receipts["in_vivo"]
        if "clinical" in receipts:
            row.clinical = receipts["clinical"]

        # Verify cell states
        cells = row.cells()
        assert cells["gdsc"].status == ModalityStatus.POSITIVE, "gdsc must be POSITIVE after receipt merge"
        assert cells["expression"].status == ModalityStatus.POSITIVE, "expression must be POSITIVE"
        assert cells["crispr"].status == ModalityStatus.MISSING, "crispr must be MISSING (not NEGATIVE)"

        # Count positives
        n_positive = sum(1 for v in cells.values() if v.status == ModalityStatus.POSITIVE)
        assert n_positive >= 2, f"Need >= 2 POSITIVE modalities, got {n_positive}"

        # Compute tier
        score = _weighted_score(row)
        tier = _assign_recommendation_tier(row, score, n_positive, rs_score=None)

        assert tier == "Strong candidate dependency axis", (
            f"MBD4/ATR_WEE1 must reach Strong tier with GDSC receipt, got: '{tier}'\n"
            f"  gdsc={cells['gdsc'].status}, expression={cells['expression'].status}, "
            f"  crispr={cells['crispr'].status}, n_positive={n_positive}, score={score:.2f}"
        )

    def test_mbd4_atr_in_vitro_still_missing(self):
        """in_vitro must remain MISSING — no isogenic data exists yet."""
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)
        assert receipts["in_vitro"].status == ModalityStatus.MISSING

    def test_mbd4_atr_clinical_still_missing(self):
        """clinical must remain MISSING — no patient-level ATRi data for MBD4-LOF."""
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)
        assert receipts["clinical"].status == ModalityStatus.MISSING

    def test_mbd4_atr_does_not_reach_validated(self):
        """Must NOT reach Validated tier — no clinical/in_vivo POSITIVE."""
        from sl_agent.multimodal.models import (
            EvidenceRow, ModalityStatus
        )
        receipts = get_literature_receipts("MBD4", CandidateAxis.ATR_WEE1)
        row = EvidenceRow(
            axis=CandidateAxis.ATR_WEE1,
            axis_label="ATR/WEE1 Inhibitors",
            mechanism="BER stress → ATR/WEE1 checkpoint dependency",
        )
        for key in ("gdsc", "expression", "in_vitro", "in_vivo", "clinical"):
            if key in receipts and getattr(row, key).status == ModalityStatus.MISSING:
                setattr(row, key, receipts[key])

        cells = row.cells()
        n_positive = sum(1 for v in cells.values() if v.status == ModalityStatus.POSITIVE)
        score = _weighted_score(row)
        tier = _assign_recommendation_tier(row, score, n_positive, rs_score=None)

        assert tier != "Validated SL therapeutic lever", (
            "MBD4/ATR_WEE1 must NOT reach Validated — no clinical/in_vivo POSITIVE data"
        )


class TestPRMT5MTAPAxis:
    """Verify PRMT5_MTAP axis is correctly registered and receipts are populated."""

    def test_prmt5_mtap_in_candidate_axis_enum(self):
        """PRMT5_MTAP must be in CandidateAxis enum."""
        assert hasattr(CandidateAxis, "PRMT5_MTAP")
        assert CandidateAxis.PRMT5_MTAP.value == "prmt5_mtap"

    def test_mtap_prmt5_receipt_crispr_positive(self):
        """MTAP/PRMT5_MTAP CRISPR receipt must be POSITIVE."""
        receipts = get_literature_receipts("MTAP", CandidateAxis.PRMT5_MTAP)
        assert "crispr" in receipts
        assert receipts["crispr"].status == ModalityStatus.POSITIVE

    def test_mtap_prmt5_receipt_in_vitro_positive(self):
        """MTAP/PRMT5_MTAP in_vitro receipt must be POSITIVE."""
        receipts = get_literature_receipts("MTAP", CandidateAxis.PRMT5_MTAP)
        assert "in_vitro" in receipts
        assert receipts["in_vitro"].status == ModalityStatus.POSITIVE

    def test_mtap_prmt5_receipt_clinical_mixed(self):
        """MTAP/PRMT5_MTAP clinical receipt must be MIXED (Phase I/II ongoing)."""
        receipts = get_literature_receipts("MTAP", CandidateAxis.PRMT5_MTAP)
        assert "clinical" in receipts
        assert receipts["clinical"].status == ModalityStatus.MIXED
