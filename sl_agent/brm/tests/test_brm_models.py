"""Tests for BrM v2 Pydantic models."""

import pytest
from pydantic import ValidationError

from sl_agent.brm.models import (
    BrMEvidenceClass,
    BrMGeneRow,
    BrMQueryInput,
    BrMRowClass,
    BrMTargetabilityMatrix,
    CascadeStep,
    ContextSpecificity,
    ExploitMode,
    RUO_DISCLAIMER,
    ScoreBasis,
    SLPartner,
    TargetNodeRole,
)


# ---------------------------------------------------------------------------
# Enum membership
# ---------------------------------------------------------------------------

class TestEnums:
    def test_exploit_mode_members(self):
        assert set(ExploitMode) == {
            ExploitMode.DIRECT,
            ExploitMode.PATHWAY,
            ExploitMode.DELIVERY,
            ExploitMode.NO_PATH,
        }

    def test_brm_evidence_class_members(self):
        expected = {
            "in_vivo_crispr_screen", "patient_expression", "genetic_suppression",
            "pharmacologic", "pathway_inference", "bbb_biology",
        }
        assert {e.value for e in BrMEvidenceClass} == expected

    def test_cascade_step_members(self):
        expected = {
            "intravasation", "circulation", "bbb_transit", "colonization",
            "immune_evasion", "dormancy_escape", "proliferation", "multi_step",
        }
        assert {e.value for e in CascadeStep} == expected

    def test_target_node_role_members(self):
        expected = {"self", "partner", "pathway_bottleneck", "microenvironment", "delivery_interface"}
        assert {e.value for e in TargetNodeRole} == expected

    def test_context_specificity_members(self):
        expected = {"brm_specific", "cns_relevant", "pan_cancer"}
        assert {e.value for e in ContextSpecificity} == expected

    def test_brm_row_class_members(self):
        assert set(BrMRowClass) == {BrMRowClass.CALIBRATION, BrMRowClass.NOVEL, BrMRowClass.NEGATIVE}

    def test_score_basis_members(self):
        assert set(ScoreBasis) == {ScoreBasis.RECEIPT, ScoreBasis.DEPMAP, ScoreBasis.HYBRID}


# ---------------------------------------------------------------------------
# SLPartner
# ---------------------------------------------------------------------------

class TestSLPartner:
    def _make(self, **kwargs):
        defaults = dict(
            partner_gene="ITGAV",
            delta_dep=-0.45,
            fdr=0.05,
            n_high=24,
            n_low=24,
            partner_quality_tier="strong",
        )
        defaults.update(kwargs)
        return SLPartner(**defaults)

    def test_valid_partner(self):
        p = self._make()
        assert p.partner_gene == "ITGAV"
        assert p.delta_dep == -0.45

    def test_fdr_out_of_range_raises(self):
        with pytest.raises(ValidationError):
            self._make(fdr=1.5)

    def test_invalid_tier_raises(self):
        with pytest.raises(ValidationError):
            self._make(partner_quality_tier="excellent")

    def test_valid_tiers(self):
        for tier in ("strong", "moderate", "weak"):
            p = self._make(partner_quality_tier=tier)
            assert p.partner_quality_tier == tier


# ---------------------------------------------------------------------------
# BrMGeneRow
# ---------------------------------------------------------------------------

def _make_row(**kwargs):
    defaults = dict(
        row_class=BrMRowClass.CALIBRATION,
        score_basis=ScoreBasis.RECEIPT,
        gene="BACE1",
        visibility_axis="absent_from_standard_DDR_panels",
        brm_evidence_class=BrMEvidenceClass.IN_VIVO_CRISPR_SCREEN,
        primary_evidence_modality="in_vivo_CRISPRa_screen_NSCLC_PDX",
        cascade_step=CascadeStep.COLONIZATION,
        cascade_step_label="Brain parenchyma colonization",
        primary_exploit_mode=ExploitMode.DIRECT,
        primary_vulnerability_hypothesis="BACE1 cleaves EGFR to drive colonization.",
        target_node="BACE1",
        target_node_role=TargetNodeRole.SELF,
        drug_class="BACE1 inhibitor",
        cns_bbb_feasibility="established",
        cns_bbb_rationale="Verubecestat is CNS-penetrant.",
        context_specificity=ContextSpecificity.BRM_SPECIFIC,
        evidence_tier="Strong",
        recommendation_status="Prioritize",
        confidence_score=0.82,
        frozen_at="2026-05-07T00:00:00+00:00",
    )
    defaults.update(kwargs)
    return BrMGeneRow(**defaults)


class TestBrMGeneRow:
    def test_valid_row(self):
        row = _make_row()
        assert row.gene == "BACE1"
        assert row.evidence_tier == "Strong"

    def test_confidence_out_of_range_raises(self):
        with pytest.raises(ValidationError):
            _make_row(confidence_score=1.5)

    def test_confidence_negative_raises(self):
        with pytest.raises(ValidationError):
            _make_row(confidence_score=-0.1)

    def test_validated_tier_raises(self):
        with pytest.raises(ValidationError):
            _make_row(evidence_tier="Validated")

    def test_secondary_exploit_modes_defaults_empty(self):
        row = _make_row()
        assert row.secondary_exploit_modes == []

    def test_emt_cluster_defaults_false(self):
        row = _make_row()
        assert row.emt_cluster is False

    def test_broad_state_marker_defaults_false(self):
        row = _make_row()
        assert row.broad_state_marker is False

    def test_sl_partners_defaults_empty(self):
        row = _make_row()
        assert row.sl_partners == []

    def test_novelty_score_optional(self):
        row = _make_row(novelty_score=0.72)
        assert row.novelty_score == 0.72

    def test_novelty_score_none_by_default(self):
        row = _make_row()
        assert row.novelty_score is None


# ---------------------------------------------------------------------------
# BrMTargetabilityMatrix
# ---------------------------------------------------------------------------

class TestBrMTargetabilityMatrix:
    def _make_matrix(self, **kwargs):
        row = _make_row()
        defaults = dict(
            query_context="NSCLC brain metastasis",
            cancer_type="NSCLC",
            rows=[row],
            frozen_at="2026-05-07T00:00:00+00:00",
        )
        defaults.update(kwargs)
        return BrMTargetabilityMatrix(**defaults)

    def test_ruo_true_by_default(self):
        m = self._make_matrix()
        assert m.ruo is True

    def test_disclaimer_non_empty(self):
        m = self._make_matrix()
        assert len(m.disclaimer) > 0
        assert "Research use only" in m.disclaimer

    def test_gold_standard_gene_default(self):
        m = self._make_matrix()
        assert m.gold_standard_gene == "BACE1"

    def test_version_v2(self):
        m = self._make_matrix()
        assert m.version == "v2"

    def test_depmap_release_default(self):
        m = self._make_matrix()
        assert m.depmap_release == "24Q4"


# ---------------------------------------------------------------------------
# BrMQueryInput
# ---------------------------------------------------------------------------

class TestBrMQueryInput:
    def test_default_genes(self):
        q = BrMQueryInput()
        assert "BACE1" in q.genes
        assert len(q.genes) == 7

    def test_default_cancer_type(self):
        q = BrMQueryInput()
        assert q.cancer_type == "NSCLC"

    def test_custom_genes(self):
        q = BrMQueryInput(genes=["ZEB1", "VIM"])
        assert q.genes == ["ZEB1", "VIM"]
