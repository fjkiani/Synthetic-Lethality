"""
Tests for the biomarker subsystem.

Covers the load-bearing honesty gates carried from the PI3K/mTOR audit:
  1. Registry integrity  — nothing is 'validated' without a passing external
     replication receipt on an INDEPENDENT cohort.
  2. Non-overlap guard    — same-cohort reprocessing is not replication.
  3. Grounding guard       — every number in an interpretation traces to a receipt.
  4. Honest status         — the OV prognostic signature ships discovery_only
     because replication across independent cohorts is INCONSISTENT.
  5. Advisory-only wiring — biomarker_context never changes drug routing.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from sl_agent.biomarkers import (
    BiomarkerType,
    BiomarkerStatus,
    ValidationReceipt,
    BiomarkerModel,
    BiomarkerContext,
    BiomarkerResult,
    BiomarkerRegistry,
)
from sl_agent.biomarkers.datasets import CohortSplit, CohortOverlapError
from sl_agent.biomarkers.grounding_guard import (
    check_interpretation_grounded,
    assert_grounded,
)


# ── fixtures: minimal receipts ────────────────────────────────────────────────

def _independent_external():
    return ValidationReceipt(
        kind="external_replication", cohort_id="GSE9891", n=268, n_events=113,
        is_independent_of_training=True, metric_name="c_index", metric_value=0.625,
        ci95_low=0.569, ci95_high=0.681, p_value=0.00001,
    )


def _same_cohort_sensitivity():
    return ValidationReceipt(
        kind="sensitivity_check", cohort_id="TCGAOVARIAN", n=521, n_events=281,
        is_independent_of_training=False, metric_name="c_index", metric_value=0.572,
    )


# ── 1. registry integrity ─────────────────────────────────────────────────────

def test_registry_refuses_validated_without_external_replication():
    reg = BiomarkerRegistry()
    m = BiomarkerModel(
        cancer="ovarian_hgsoc", biomarker_type=BiomarkerType.PROGNOSTIC,
        method_id="fake_validated", status=BiomarkerStatus.VALIDATED,
        receipts=[_same_cohort_sensitivity()],  # NOT independent
        interpretation="",
    )
    with pytest.raises(ValueError):
        reg.register(m)


def test_registry_get_excludes_discovery_by_default():
    reg = BiomarkerRegistry()
    disc = BiomarkerModel(
        cancer="ovarian_hgsoc", biomarker_type=BiomarkerType.PROGNOSTIC,
        method_id="disc_only", status=BiomarkerStatus.DISCOVERY_ONLY,
        receipts=[_independent_external()], interpretation="",
    )
    reg.register(disc)
    assert reg.get("ovarian_hgsoc", BiomarkerType.PROGNOSTIC) == []
    got = reg.get("ovarian_hgsoc", BiomarkerType.PROGNOSTIC, include_discovery=True)
    assert len(got) == 1 and got[0].method_id == "disc_only"


def test_validated_model_requires_independent_receipt_flag():
    # A properly-validated model (independent external replication present).
    m = BiomarkerModel(
        cancer="x", biomarker_type=BiomarkerType.PROGNOSTIC, method_id="ok",
        status=BiomarkerStatus.VALIDATED, receipts=[_independent_external()],
        interpretation="",
    )
    assert m.is_validated() is True
    # Strip independence -> no longer validated.
    m2 = BiomarkerModel(
        cancer="x", biomarker_type=BiomarkerType.PROGNOSTIC, method_id="bad",
        status=BiomarkerStatus.VALIDATED, receipts=[_same_cohort_sensitivity()],
        interpretation="",
    )
    assert m2.is_validated() is False


# ── 2. non-overlap guard ──────────────────────────────────────────────────────

def test_non_overlap_guard_raises_on_identical_cohort():
    with pytest.raises(CohortOverlapError):
        CohortSplit("ov_tcga", "ov_tcga").assert_non_overlap()


def test_non_overlap_guard_raises_on_shared_tcga_patient():
    # Same patient, two sample barcodes -> must be caught (patient-level norm).
    split = CohortSplit(
        "trainA", "valB",
        train_ids=["TCGA-13-0720-01A"],
        val_ids=["TCGA-13-0720-02B"],
    )
    with pytest.raises(CohortOverlapError):
        split.assert_non_overlap()


def test_non_overlap_guard_passes_for_distinct_studies():
    split = CohortSplit(
        "TCGAOVARIAN", "GSE9891",
        train_ids=["TCGA-13-0720"], val_ids=["GSM249012"],
    )
    assert split.assert_non_overlap().is_independent is True


# ── 3. grounding guard ────────────────────────────────────────────────────────

def test_grounding_guard_flags_fabricated_number():
    m = BiomarkerModel(
        cancer="x", biomarker_type=BiomarkerType.PROGNOSTIC, method_id="g",
        receipts=[_independent_external()],
        interpretation="The hazard ratio was 1.30 with p=0.023.",  # not in receipts
    )
    viol = check_interpretation_grounded(m)
    assert "1.30" in viol and "0.023" in viol


def test_grounding_guard_accepts_receipt_backed_numbers():
    m = BiomarkerModel(
        cancer="x", biomarker_type=BiomarkerType.PROGNOSTIC, method_id="g",
        receipts=[_independent_external()],
        interpretation="C-index 0.625 in GSE9891 (n=268, 113 events).",
    )
    assert check_interpretation_grounded(m) == []
    assert_grounded(m)  # must not raise


# ── 4. honest status from real receipts (inconsistent replication) ────────────

RECEIPT_DIR = Path("/mnt/shared-workspace/shared/biomarker_val")


@pytest.mark.skipif(
    not (RECEIPT_DIR / "w2_prognostic_replication_results.json").exists(),
    reason="W2 receipts not present in this environment",
)
def test_ov_prognostic_ships_discovery_only_due_to_inconsistent_replication():
    from sl_agent.biomarkers.assemble import build_prognostic_model
    m = build_prognostic_model(
        RECEIPT_DIR / "w1_prognostic_results.json",
        RECEIPT_DIR / "w2_prognostic_replication_results.json",
    )
    # Signal in GSE9891 but null in GSE32062 -> NOT consistently replicated.
    assert m.status == BiomarkerStatus.DISCOVERY_ONLY
    assert m.is_validated() is False
    # Interpretation must be honest and grounded.
    assert_grounded(m)
    low = m.interpretation.lower()
    assert "inconsistent" in low
    assert "not a validated prognostic" in low
    # Must NOT overstate two same-source cohorts as independent replication:
    # the only 'independent' claims must refer to the GEO cohorts, and the
    # TCGA re-run must be labelled a sensitivity check, not replication.
    assert "sensitivity check" in low


# ── 5. advisory-only integration ──────────────────────────────────────────────

def test_biomarker_context_is_advisory_and_never_changes_routing():
    from sl_agent.multimodal.models import EvidenceMatrix, EvidenceRow, CandidateAxis

    row = EvidenceRow(axis=CandidateAxis.PARP_INHIBITORS, axis_label="PARP inhibitors",
                      mechanism="synthetic lethality with HR deficiency",
                      recommendation_tier="Strong")
    em = EvidenceMatrix(query_gene="BRCA1", cancer_type="ovarian_hgsoc", rows=[row])
    baseline = em.recommendation_summary()

    ctx = BiomarkerContext(
        cancer="ovarian_hgsoc",
        results=[BiomarkerResult(
            model_ref="riester2014_prognostic_ov",
            biomarker_type=BiomarkerType.PROGNOSTIC,
            status=BiomarkerStatus.DISCOVERY_ONLY,
            interpretation="discovery-only prognostic context",
        )],
    )
    em.biomarker_context = ctx  # attach advisory context

    # Attaching biomarker context must not change any recommendation tier.
    assert em.recommendation_summary() == baseline
    assert em.biomarker_context.advisory_only is True
    assert em.biomarker_context.results[0].status == BiomarkerStatus.DISCOVERY_ONLY


def test_evidence_matrix_biomarker_context_defaults_none():
    from sl_agent.multimodal.models import EvidenceMatrix
    em = EvidenceMatrix(query_gene="TP53")
    assert em.biomarker_context is None  # additive, back-compatible
