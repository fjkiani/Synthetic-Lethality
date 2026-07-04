"""
Breast Cancer SL Framework Expansion Tests — 2026-07-03

Tests the 5-modality evidence engine for breast cancer gene-axis pairs.
Covers:
  - VALIDATED tier: BRCA1/2 (PARP), PIK3CA (PI3K_AKT), PALB2 (PARP)
  - STRONG tier: ARID1A (ATR_WEE1)
  - MECHANISTIC tier: RB1, FBXW7, CDH1, PIK3R1

Framework audit:
  Modality weights: clinical=4, in_vivo=3, in_vitro=2, prism=2, gdsc=2, expression=1, crispr=1
  VALIDATED: strong_clinical AND n_positive >= 3
  STRONG: (crispr_pos OR pharma_pos) AND n_positive >= 2 AND NOT crispr_neg
  MECHANISTIC: expression/pathway only
  RS promotion: only ATR_WEE1 + PKMYT1; only Mechanistic→Strong; only when RS=high AND independent evidence
"""
import pytest
from sl_agent.multimodal.models import (
    CandidateAxis, ModalityStatus, EvidenceRow
)
from sl_agent.multimodal.literature_receipts import (
    get_literature_receipts, list_receipts_for_gene
)
from sl_agent.multimodal.modality_fuser import _assign_recommendation_tier, _weighted_score


# ── Helper ────────────────────────────────────────────────────────────────────

def build_row_with_receipts(gene: str, axis: CandidateAxis) -> EvidenceRow:
    """Build EvidenceRow with frozen receipts applied (simulates matrix_builder Step 4)."""
    receipts = get_literature_receipts(gene, axis)
    row = EvidenceRow(
        axis=axis,
        axis_label=axis.value,
        mechanism=f"{gene} x {axis.value}",
    )
    for key in ("gdsc", "crispr", "in_vitro", "in_vivo", "clinical"):
        if key in receipts and getattr(row, key).status == ModalityStatus.MISSING:
            setattr(row, key, receipts[key])
    if "expression" in receipts and row.expression.status == ModalityStatus.MISSING:
        row.expression = receipts["expression"]
    return row


def get_tier(gene: str, axis: CandidateAxis) -> tuple:
    row = build_row_with_receipts(gene, axis)
    cells = row.cells()
    n_positive = sum(1 for v in cells.values() if v.status == ModalityStatus.POSITIVE)
    score = _weighted_score(row)
    tier = _assign_recommendation_tier(row, score, n_positive, rs_score=None)
    return tier, n_positive, score


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATED TIER TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestValidatedTierBreastCancer:
    """BRCA1/2 and PIK3CA/PALB2 must reach Validated tier."""

    def test_brca1_parp_validated(self):
        """BRCA1 + PARP_INHIBITORS must be Validated (OlympiAD Phase III + isogenic + CRISPR)."""
        tier, n_pos, score = get_tier("BRCA1", CandidateAxis.PARP_INHIBITORS)
        assert tier == "Validated SL therapeutic lever", (
            f"BRCA1/PARP must be Validated, got: {tier} (n_pos={n_pos}, score={score:.1f})"
        )
        assert n_pos >= 3, f"Need >= 3 POSITIVE modalities, got {n_pos}"
        assert score >= 10.0, f"Score must be >= 10 for Validated, got {score:.1f}"

    def test_brca1_parp_has_clinical_positive(self):
        """BRCA1/PARP clinical receipt must be POSITIVE (OlympiAD)."""
        receipts = get_literature_receipts("BRCA1", CandidateAxis.PARP_INHIBITORS)
        assert receipts["clinical"].status == ModalityStatus.POSITIVE
        assert "OlympiAD" in receipts["clinical"].summary or "28578601" in receipts["clinical"].pmids

    def test_brca1_parp_has_in_vivo_positive(self):
        """BRCA1/PARP in_vivo receipt must be POSITIVE (PDX data)."""
        receipts = get_literature_receipts("BRCA1", CandidateAxis.PARP_INHIBITORS)
        assert receipts["in_vivo"].status == ModalityStatus.POSITIVE

    def test_brca2_parp_validated(self):
        """BRCA2 + PARP_INHIBITORS must be Validated."""
        tier, n_pos, score = get_tier("BRCA2", CandidateAxis.PARP_INHIBITORS)
        assert tier == "Validated SL therapeutic lever", (
            f"BRCA2/PARP must be Validated, got: {tier}"
        )
        assert n_pos >= 3

    def test_pik3ca_pi3k_akt_validated(self):
        """PIK3CA + PI3K_AKT must be Validated (SOLAR-1 Phase III + in vitro + CRISPR)."""
        tier, n_pos, score = get_tier("PIK3CA", CandidateAxis.PI3K_AKT)
        assert tier == "Validated SL therapeutic lever", (
            f"PIK3CA/PI3K_AKT must be Validated, got: {tier} (n_pos={n_pos}, score={score:.1f})"
        )
        assert n_pos >= 3

    def test_pik3ca_pi3k_akt_has_clinical_positive(self):
        """PIK3CA/PI3K_AKT clinical receipt must be POSITIVE (SOLAR-1)."""
        receipts = get_literature_receipts("PIK3CA", CandidateAxis.PI3K_AKT)
        assert receipts["clinical"].status == ModalityStatus.POSITIVE
        assert "SOLAR-1" in receipts["clinical"].summary or "31091374" in receipts["clinical"].pmids

    def test_palb2_parp_validated(self):
        """PALB2 + PARP_INHIBITORS must be Validated (TBCRC 048 Phase II + in vitro)."""
        tier, n_pos, score = get_tier("PALB2", CandidateAxis.PARP_INHIBITORS)
        assert tier == "Validated SL therapeutic lever", (
            f"PALB2/PARP must be Validated, got: {tier} (n_pos={n_pos}, score={score:.1f})"
        )
        assert n_pos >= 3

    def test_palb2_parp_has_clinical_positive(self):
        """PALB2/PARP clinical receipt must be POSITIVE (TBCRC 048)."""
        receipts = get_literature_receipts("PALB2", CandidateAxis.PARP_INHIBITORS)
        assert receipts["clinical"].status == ModalityStatus.POSITIVE
        assert "33119476" in receipts["clinical"].pmids  # Tung et al. JCO 2020 (DOI: 10.1200/jco.20.02151)


# ══════════════════════════════════════════════════════════════════════════════
# STRONG TIER TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestStrongTierBreastCancer:
    """ARID1A/ATR_WEE1 must reach Strong tier."""

    def test_arid1a_atr_strong(self):
        """ARID1A + ATR_WEE1 must be Strong (CRISPR + in_vitro + expression)."""
        tier, n_pos, score = get_tier("ARID1A", CandidateAxis.ATR_WEE1)
        assert tier == "Strong candidate dependency axis", (
            f"ARID1A/ATR_WEE1 must be Strong, got: {tier} (n_pos={n_pos}, score={score:.1f})"
        )
        assert n_pos >= 2

    def test_arid1a_atr_has_crispr_positive(self):
        """ARID1A/ATR_WEE1 CRISPR receipt must be POSITIVE (Williamson 2016)."""
        receipts = get_literature_receipts("ARID1A", CandidateAxis.ATR_WEE1)
        assert receipts["crispr"].status == ModalityStatus.POSITIVE
        assert "27958275" in receipts["crispr"].pmids  # Williamson et al. Nat Commun 2016 (DOI: 10.1038/ncomms13837)

    def test_arid1a_atr_has_in_vitro_positive(self):
        """ARID1A/ATR_WEE1 in_vitro receipt must be POSITIVE."""
        receipts = get_literature_receipts("ARID1A", CandidateAxis.ATR_WEE1)
        assert receipts["in_vitro"].status == ModalityStatus.POSITIVE

    def test_arid1a_atr_does_not_reach_validated(self):
        """ARID1A/ATR_WEE1 must NOT reach Validated — no clinical/in_vivo POSITIVE."""
        tier, _, _ = get_tier("ARID1A", CandidateAxis.ATR_WEE1)
        assert tier != "Validated SL therapeutic lever"


# ══════════════════════════════════════════════════════════════════════════════
# MECHANISTIC TIER TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestMechanisticTierBreastCancer:
    """Novel hypotheses must reach Mechanistic tier (expression only)."""

    def test_rb1_atr_mechanistic(self):
        """RB1 + ATR_WEE1 must be Mechanistic (expression only — no CRISPR/in_vitro yet)."""
        tier, n_pos, score = get_tier("RB1", CandidateAxis.ATR_WEE1)
        assert tier == "Mechanistic candidate only", (
            f"RB1/ATR_WEE1 must be Mechanistic, got: {tier}"
        )
        assert n_pos == 1  # expression only

    def test_fbxw7_pkmyt1_mechanistic(self):
        """FBXW7 + PKMYT1 must be Mechanistic (cyclin E stabilization hypothesis)."""
        tier, n_pos, score = get_tier("FBXW7", CandidateAxis.PKMYT1)
        assert tier == "Mechanistic candidate only", (
            f"FBXW7/PKMYT1 must be Mechanistic, got: {tier}"
        )

    def test_cdh1_atr_mechanistic(self):
        """CDH1 + ATR_WEE1 must be Mechanistic (EMT-driven RS hypothesis)."""
        tier, n_pos, score = get_tier("CDH1", CandidateAxis.ATR_WEE1)
        assert tier == "Mechanistic candidate only", (
            f"CDH1/ATR_WEE1 must be Mechanistic, got: {tier}"
        )

    def test_pik3r1_pi3k_akt_mechanistic(self):
        """PIK3R1 + PI3K_AKT must be Mechanistic (regulatory subunit LOF hypothesis)."""
        tier, n_pos, score = get_tier("PIK3R1", CandidateAxis.PI3K_AKT)
        assert tier == "Mechanistic candidate only", (
            f"PIK3R1/PI3K_AKT must be Mechanistic, got: {tier}"
        )

    def test_rb1_in_vitro_missing(self):
        """RB1/ATR_WEE1 in_vitro must be MISSING — no isogenic data yet."""
        receipts = get_literature_receipts("RB1", CandidateAxis.ATR_WEE1)
        assert receipts["in_vitro"].status == ModalityStatus.MISSING

    def test_fbxw7_in_vitro_missing(self):
        """FBXW7/PKMYT1 in_vitro must be MISSING — no RP-6306 data yet."""
        receipts = get_literature_receipts("FBXW7", CandidateAxis.PKMYT1)
        assert receipts["in_vitro"].status == ModalityStatus.MISSING


# ══════════════════════════════════════════════════════════════════════════════
# PI3K_AKT AXIS TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestPI3KAKTAxis:
    """Verify PI3K_AKT axis is correctly registered."""

    def test_pi3k_akt_in_candidate_axis_enum(self):
        """PI3K_AKT must be in CandidateAxis enum."""
        assert hasattr(CandidateAxis, "PI3K_AKT")
        assert CandidateAxis.PI3K_AKT.value == "pi3k_akt"

    def test_pi3k_akt_distinct_from_parp_inhibitors(self):
        """PI3K_AKT must be a distinct axis from PARP_INHIBITORS."""
        assert CandidateAxis.PI3K_AKT != CandidateAxis.PARP_INHIBITORS

    def test_pik3r1_pi3k_akt_has_expression_positive(self):
        """PIK3R1/PI3K_AKT expression receipt must be POSITIVE."""
        receipts = get_literature_receipts("PIK3R1", CandidateAxis.PI3K_AKT)
        assert "expression" in receipts
        assert receipts["expression"].status == ModalityStatus.POSITIVE


# ══════════════════════════════════════════════════════════════════════════════
# FRAMEWORK AUDIT TESTS
# ══════════════════════════════════════════════════════════════════════════════

class TestFrameworkAudit:
    """Audit the 5-modality promotion framework invariants."""

    def test_validated_requires_clinical_or_in_vivo_positive(self):
        """Validated tier requires strong_clinical (clinical OR in_vivo POSITIVE)."""
        # ARID1A has CRISPR + in_vitro + expression but no clinical/in_vivo → Strong, not Validated
        tier, _, _ = get_tier("ARID1A", CandidateAxis.ATR_WEE1)
        assert tier == "Strong candidate dependency axis"
        assert tier != "Validated SL therapeutic lever"

    def test_strong_requires_crispr_or_pharma_positive(self):
        """Strong tier requires CRISPR or pharmacologic (GDSC/PRISM) POSITIVE."""
        # RB1 has only expression → Mechanistic, not Strong
        tier, _, _ = get_tier("RB1", CandidateAxis.ATR_WEE1)
        assert tier == "Mechanistic candidate only"

    def test_mechanistic_does_not_require_crispr(self):
        """Mechanistic tier fires on expression alone (no CRISPR required)."""
        tier, n_pos, _ = get_tier("CDH1", CandidateAxis.ATR_WEE1)
        assert tier == "Mechanistic candidate only"
        assert n_pos == 1  # expression only

    def test_total_breast_cancer_receipts_count(self):
        """Total receipt count must include all breast cancer additions."""
        from sl_agent.multimodal.literature_receipts import _FROZEN_RECEIPTS
        bc_genes = {"BRCA1", "BRCA2", "PIK3CA", "PALB2", "ARID1A", "RB1", "FBXW7", "CDH1", "PIK3R1"}
        bc_receipts = [(g, a) for (g, a) in _FROZEN_RECEIPTS.keys() if g in bc_genes]
        assert len(bc_receipts) == 9, (
            f"Expected 9 breast cancer receipts, got {len(bc_receipts)}: {bc_receipts}"
        )
