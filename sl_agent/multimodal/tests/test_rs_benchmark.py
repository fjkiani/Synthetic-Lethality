"""
RS (Replication Stress) benchmark tests.

Tests:
  - RS scoring unit tests (compute_rs_score with various feature combos)
  - RS tier promotion (fuse_matrix with RS scores and evidence rows)
  - Regression guards (non-RS axes unaffected, RS cannot reach Validated)

Total: 12 tests
"""
from __future__ import annotations

import pytest

from sl_agent.multimodal.replication_stress import RSFeatureSet, RSScore, compute_rs_score
from sl_agent.multimodal.models import (
    CandidateAxis,
    EvidenceMatrix,
    EvidenceRow,
    ModalityEvidence,
    ModalityStatus,
    Modality,
)
from sl_agent.multimodal.modality_fuser import fuse_matrix


# ─────────────────────────────────────────────────────────────────────────────
# Class 1: RS Scoring Unit Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestRSScoring:

    def test_mbd4_alone(self):
        """MBD4 alone: rs=0.05, level='low'"""
        features = RSFeatureSet(mbd4_lof=True)
        result = compute_rs_score(features)
        assert abs(result.score - 0.05) < 1e-9
        assert result.level == "low"

    def test_mbd4_tp53(self):
        """MBD4+TP53: 0.05+0.10=0.15, level='low'"""
        features = RSFeatureSet(mbd4_lof=True, tp53_lof=True)
        result = compute_rs_score(features)
        assert abs(result.score - 0.15) < 1e-9
        assert result.level == "low"

    def test_mbd4_tp53_arid1a(self):
        """MBD4+TP53+ARID1A: 0.05+0.10+0.20=0.35, level='moderate'"""
        features = RSFeatureSet(mbd4_lof=True, tp53_lof=True, arid1a_lof=True)
        result = compute_rs_score(features)
        assert abs(result.score - 0.35) < 1e-9
        assert result.level == "moderate"

    def test_ccne1_tp53(self):
        """CCNE1+TP53: 0.30+0.10=0.40, level='high'"""
        features = RSFeatureSet(ccne1_amplified=True, tp53_lof=True)
        result = compute_rs_score(features)
        assert abs(result.score - 0.40) < 1e-9
        assert result.level == "high"

    def test_ccne1_myc_arid1a(self):
        """CCNE1+MYC+ARID1A: 0.30+0.20+0.20=0.70, level='high'"""
        features = RSFeatureSet(ccne1_amplified=True, myc_amplified=True, arid1a_lof=True)
        result = compute_rs_score(features)
        assert result.score <= 1.0
        assert result.level == "high"

    def test_no_features(self):
        """No features: score=0.0, level='none'"""
        features = RSFeatureSet()
        result = compute_rs_score(features)
        assert result.score == 0.0
        assert result.level == "none"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers — build minimal rows for tier-promotion tests
# ─────────────────────────────────────────────────────────────────────────────

def _make_row(
    axis: CandidateAxis,
    label: str,
    crispr_status: ModalityStatus = ModalityStatus.MISSING,
    prism_status: ModalityStatus = ModalityStatus.MISSING,
    gdsc_status: ModalityStatus = ModalityStatus.MISSING,
    expression_status: ModalityStatus = ModalityStatus.POSITIVE,
    in_vitro_status: ModalityStatus = ModalityStatus.MISSING,
    in_vivo_status: ModalityStatus = ModalityStatus.MISSING,
    clinical_status: ModalityStatus = ModalityStatus.MISSING,
) -> EvidenceRow:
    """Build a minimal EvidenceRow with specified modality statuses."""
    row = EvidenceRow(
        axis=axis,
        axis_label=label,
        mechanism="test mechanism",
    )
    row.crispr = ModalityEvidence(
        modality=Modality.CRISPR_DEPENDENCY,
        status=crispr_status,
        delta_dep=-0.15 if crispr_status == ModalityStatus.POSITIVE else 0.03,
    )
    row.prism = ModalityEvidence(modality=Modality.PRISM_PHARMACOLOGIC, status=prism_status)
    row.gdsc = ModalityEvidence(modality=Modality.GDSC_PHARMACOLOGIC, status=gdsc_status)
    row.expression = ModalityEvidence(modality=Modality.EXPRESSION_ASSOC, status=expression_status)
    row.in_vitro = ModalityEvidence(modality=Modality.IN_VITRO_FUNCTIONAL, status=in_vitro_status)
    row.in_vivo = ModalityEvidence(modality=Modality.IN_VIVO_PDX, status=in_vivo_status)
    row.clinical = ModalityEvidence(modality=Modality.CLINICAL, status=clinical_status)
    return row


# ─────────────────────────────────────────────────────────────────────────────
# Class 2: RS Tier Promotion Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestRSTierPromotion:

    def test_ccne1_pkmyt1_promotes_to_strong(self):
        """CCNE1 with PKMYT1 axis, high RS + CRISPR receipt → Strong"""
        # Build row: PKMYT1 axis, CRISPR positive, expression positive
        # Other modalities MISSING → base tier would be Mechanistic or Strong
        # With CRISPR=POSITIVE + expression=POSITIVE → n_pos=2 → Strong already
        # But let's test with CRISPR=POSITIVE only (n_pos=1, expression=MISSING)
        # to verify RS promotion from Mechanistic → Strong
        row = _make_row(
            axis=CandidateAxis.PKMYT1,
            label="PKMYT1",
            crispr_status=ModalityStatus.POSITIVE,
            expression_status=ModalityStatus.MISSING,
        )
        matrix = EvidenceMatrix(query_gene="CCNE1", rows=[row])

        # High RS: CCNE1 + TP53 → 0.40 → "high"
        rs = compute_rs_score(RSFeatureSet(ccne1_amplified=True, tp53_lof=True))
        assert rs.level == "high"

        fused = fuse_matrix(matrix, rs_score=rs)
        pkmyt1_row = fused.rows[0]
        # CRISPR positive + high RS → should be promoted to Strong
        assert pkmyt1_row.recommendation_tier == "Strong candidate dependency axis"

    def test_mbd4_alone_stays_mechanistic(self):
        """MBD4 alone: rs=low → ATR/WEE1 stays Mechanistic (no promotion)"""
        row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="MBD4", rows=[row])

        # Low RS: MBD4 alone → 0.05 → "low"
        rs = compute_rs_score(RSFeatureSet(mbd4_lof=True))
        assert rs.level == "low"

        fused = fuse_matrix(matrix, rs_score=rs)
        atr_row = fused.rows[0]
        assert atr_row.recommendation_tier == "Mechanistic candidate only"

    def test_moderate_rs_no_pharma_stays_mechanistic(self):
        """MBD4+TP53+ARID1A: rs=moderate but no pharma/CRISPR → stays Mechanistic"""
        row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="MBD4", rows=[row])

        # Moderate RS: MBD4 + TP53 + ARID1A → 0.35 → "moderate"
        rs = compute_rs_score(RSFeatureSet(mbd4_lof=True, tp53_lof=True, arid1a_lof=True))
        assert rs.level == "moderate"

        fused = fuse_matrix(matrix, rs_score=rs)
        atr_row = fused.rows[0]
        # Moderate RS, no CRISPR/pharma → stays Mechanistic
        assert atr_row.recommendation_tier == "Mechanistic candidate only"

    def test_high_rs_no_evidence_cannot_promote(self):
        """Sabotage: high RS (CCNE1+TP53) but all modalities MISSING → no promotion"""
        row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.MISSING,
        )
        matrix = EvidenceMatrix(query_gene="MBD4", rows=[row])

        # High RS: CCNE1 + TP53 → 0.40 → "high"
        rs = compute_rs_score(RSFeatureSet(ccne1_amplified=True, tp53_lof=True))
        assert rs.level == "high"

        fused = fuse_matrix(matrix, rs_score=rs)
        atr_row = fused.rows[0]
        # High RS but NO CRISPR/pharma evidence → cannot promote
        assert atr_row.recommendation_tier != "Strong candidate dependency axis"


# ─────────────────────────────────────────────────────────────────────────────
# Class 3: RS Regression Guards
# ─────────────────────────────────────────────────────────────────────────────

class TestRSRegressionGuard:

    def test_non_rs_axis_unaffected(self):
        """PARP axis should not be affected by RS score"""
        row = _make_row(
            axis=CandidateAxis.PARP_INHIBITORS,
            label="PARP Inhibitors",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        # Fuse without RS
        matrix_no_rs = EvidenceMatrix(query_gene="MBD4", rows=[row])
        fused_no_rs = fuse_matrix(matrix_no_rs)
        tier_no_rs = fused_no_rs.rows[0].recommendation_tier

        # Build fresh row (fuse_matrix mutates in place)
        row2 = _make_row(
            axis=CandidateAxis.PARP_INHIBITORS,
            label="PARP Inhibitors",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        # Fuse with high RS
        rs = compute_rs_score(RSFeatureSet(ccne1_amplified=True, tp53_lof=True))
        matrix_with_rs = EvidenceMatrix(query_gene="MBD4", rows=[row2])
        fused_with_rs = fuse_matrix(matrix_with_rs, rs_score=rs)
        tier_with_rs = fused_with_rs.rows[0].recommendation_tier

        # PARP axis tier must be identical with or without RS
        assert tier_no_rs == tier_with_rs

    def test_rs_cannot_reach_validated(self):
        """RS modulation must never promote directly to Validated"""
        # Build a PKMYT1 row with CRISPR positive + pharma positive + expression
        # But NO clinical/in_vivo → cannot reach Validated via base logic
        # High RS should promote to Strong at most, never Validated
        row = _make_row(
            axis=CandidateAxis.PKMYT1,
            label="PKMYT1",
            crispr_status=ModalityStatus.POSITIVE,
            prism_status=ModalityStatus.POSITIVE,
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="CCNE1", rows=[row])

        # Max RS
        rs = compute_rs_score(RSFeatureSet(
            ccne1_amplified=True, myc_amplified=True, arid1a_lof=True,
            tp53_lof=True, msi_h=True, high_ploidy=True, mbd4_lof=True,
        ))
        assert rs.level == "high"

        fused = fuse_matrix(matrix, rs_score=rs)
        tier = fused.rows[0].recommendation_tier
        assert tier != "Validated SL therapeutic lever", (
            f"RS must NEVER promote to Validated, but got: {tier}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Class 4: Combinatorial Vulnerability Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestCombinatorialVulnerability:

    def test_mbd4_combinatorial_identified_no_crispr(self):
        """
        Fixture 12: MBD4 BER × checkpoint — no ATR CRISPR data.
        assess_combinatorial_vulnerability returns a CombinationVulnerability
        with combined_tier=None and non-empty data_gaps.
        """
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability
        from sl_agent.multimodal.models import CombinationVulnerability

        # Build matrix with ATR/WEE1 row, CRISPR MISSING
        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="MBD4", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet(mbd4_lof=True, tp53_lof=True))
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("MBD4", rs, fused)

        assert combo is not None, "Should identify combinatorial vulnerability for MBD4"
        assert isinstance(combo, CombinationVulnerability)
        assert combo.pathway_a_gene == "MBD4"
        assert combo.pathway_b == "checkpoint_inhibition"
        assert "ATR" in combo.pathway_b_targets
        assert combo.combined_tier is None, "No combined tier without CRISPR data"
        assert len(combo.data_gaps) > 0, "Should have data gaps"
        assert combo.rs_score is not None
        assert combo.convergence_mechanism == "mitotic_catastrophe"

    def test_arid1a_combinatorial_realized(self):
        """
        Fixture 13: ARID1A combinatorial — has CRISPR POSITIVE → combined_tier realized.
        assess_combinatorial_vulnerability returns non-None combined_tier.
        """
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability
        from sl_agent.multimodal.models import CombinationVulnerability

        # ARID1A has ATR CRISPR POSITIVE
        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.POSITIVE,
            expression_status=ModalityStatus.POSITIVE,
            in_vitro_status=ModalityStatus.POSITIVE,
            in_vivo_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="ARID1A", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet(arid1a_lof=True, tp53_lof=True))
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("ARID1A", rs, fused)

        assert combo is not None
        assert combo.pathway_a_gene == "ARID1A"
        assert combo.combined_tier is not None, "CRISPR data present → tier should be realized"
        assert len(combo.data_gaps) == 0, "No gaps when all evidence present"

    def test_combinatorial_stored_on_matrix(self):
        """
        build_evidence_matrix wires combinatorial assessment onto matrix.combination_vulnerabilities
        when rs_features are provided and ATR/WEE1 axis is present.
        """
        from sl_agent.multimodal.matrix_builder import build_evidence_matrix
        from sl_agent.multimodal.models import MultiModalQueryInput, CandidateAxis
        from sl_agent.multimodal.replication_stress import RSFeatureSet

        query = MultiModalQueryInput(
            gene="MBD4",
            mutation_type="loss_of_function",
            axes=[CandidateAxis.ATR_WEE1, CandidateAxis.CYTIDINE_ANALOGS],
            include_pharmacologic_stratification=False,
            rs_features=RSFeatureSet(mbd4_lof=True, tp53_lof=True),
        )
        matrix = build_evidence_matrix(query)

        # RS score should be stored
        assert matrix.rs_score is not None
        assert matrix.rs_score.level in ("low", "moderate", "high")

        # Combinatorial vulnerabilities should be assessed
        assert isinstance(matrix.combination_vulnerabilities, list)
        # With MBD4+TP53 RS (low), should still identify the combinatorial structure
        assert len(matrix.combination_vulnerabilities) >= 1

    def test_zero_rs_no_combinatorial(self):
        """Gene with zero RS features should not produce combinatorial assessment."""
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability

        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="UNKNOWNGENE", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet())  # No features → score = 0.0
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("UNKNOWNGENE", rs, fused)
        assert combo is None, "Zero RS score → no combinatorial vulnerability"
