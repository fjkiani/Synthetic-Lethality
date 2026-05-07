"""Tests for frozen BrM literature receipts."""

import pytest

from sl_agent.brm.brm_receipts import (
    all_receipts,
    get_brm_receipt,
    list_receipt_genes,
)


PANEL_GENES = ["BACE1", "MMP9", "CLDN5", "CCL2", "ICAM1", "MMP2", "TWIST1"]


class TestReceiptCoverage:
    def test_all_panel_genes_have_receipts(self):
        for gene in PANEL_GENES:
            r = get_brm_receipt(gene)
            assert r is not None, f"No receipt for {gene}"

    def test_case_insensitive_lookup(self):
        assert get_brm_receipt("bace1") is not None
        assert get_brm_receipt("Bace1") is not None
        assert get_brm_receipt("BACE1") is not None

    def test_unknown_gene_returns_none(self):
        assert get_brm_receipt("NOTAREALGENE99") is None

    def test_list_receipt_genes_contains_all_panel(self):
        genes = list_receipt_genes()
        for gene in PANEL_GENES:
            assert gene in genes

    def test_all_receipts_returns_dict(self):
        receipts = all_receipts()
        assert isinstance(receipts, dict)
        assert len(receipts) >= len(PANEL_GENES)


class TestBACE1Receipt:
    def setup_method(self):
        self.r = get_brm_receipt("BACE1")

    def test_chafe_2025_pmid_present(self):
        assert "10.1126/scitranslmed.adu2459" in self.r.supporting_pmids

    def test_primary_exploit_mode_direct(self):
        assert self.r.primary_exploit_mode == "direct"

    def test_delivery_in_secondary(self):
        assert "delivery" in self.r.secondary_exploit_modes

    def test_cascade_step_colonization(self):
        assert self.r.cascade_step == "colonization"

    def test_target_node_self(self):
        assert self.r.target_node == "BACE1"
        assert self.r.target_node_role == "self"

    def test_cns_bbb_established(self):
        assert self.r.cns_bbb_feasibility == "established"

    def test_evidence_tier_not_validated(self):
        assert self.r.evidence_tier != "Validated"

    def test_evidence_tier_strong(self):
        assert self.r.evidence_tier == "Strong"

    def test_data_gaps_non_empty(self):
        assert len(self.r.data_gaps) > 0

    def test_context_specificity_brm(self):
        assert self.r.context_specificity == "brm_specific"


class TestMMP9Receipt:
    def setup_method(self):
        self.r = get_brm_receipt("MMP9")

    def test_cascade_step_bbb_transit(self):
        assert self.r.cascade_step == "bbb_transit"

    def test_primary_exploit_delivery(self):
        assert self.r.primary_exploit_mode == "delivery"

    def test_has_pmids(self):
        assert len(self.r.supporting_pmids) >= 1

    def test_evidence_tier_not_validated(self):
        assert self.r.evidence_tier != "Validated"


class TestTWIST1Receipt:
    def setup_method(self):
        self.r = get_brm_receipt("TWIST1")

    def test_target_node_brd4(self):
        assert self.r.target_node == "BRD4"

    def test_primary_exploit_pathway(self):
        assert self.r.primary_exploit_mode == "pathway"

    def test_evidence_tier_not_validated(self):
        assert self.r.evidence_tier != "Validated"


class TestICAM1Receipt:
    def setup_method(self):
        self.r = get_brm_receipt("ICAM1")

    def test_broad_state_marker_flag(self):
        assert self.r.broad_state_marker is True

    def test_evidence_tier_not_validated(self):
        assert self.r.evidence_tier != "Validated"


class TestAllReceiptsInvariants:
    def test_no_receipt_has_validated_tier(self):
        for gene, r in all_receipts().items():
            assert r.evidence_tier != "Validated", f"{gene} has 'Validated' tier"

    def test_all_receipts_have_pmids(self):
        for gene, r in all_receipts().items():
            assert len(r.supporting_pmids) >= 1, f"{gene} has no PMIDs"

    def test_all_receipts_have_context_specificity(self):
        for gene, r in all_receipts().items():
            assert r.context_specificity, f"{gene} missing context_specificity"

    def test_non_strong_receipts_have_data_gaps(self):
        for gene, r in all_receipts().items():
            if r.evidence_tier != "Strong":
                assert len(r.data_gaps) > 0, f"{gene} (tier={r.evidence_tier}) has no data_gaps"
