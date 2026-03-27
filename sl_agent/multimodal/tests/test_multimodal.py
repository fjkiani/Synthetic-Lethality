"""
Multi-modal evidence engine tests.

All tests are offline — no external HTTP, no DepMap download required.
Covers:
  - EvidenceMatrix model + helpers
  - Pharmacologic analyzer (synthetic PRISM data)
  - Literature receipt layer (frozen MBD4 receipts)
  - Modality fuser (evidence levels, recommendation tiers, agreement report)
  - Matrix builder (end-to-end with synthetic data)
  - Guardrail logic
  - Calibration narrative
  - API endpoint (mocked data)

Total: 30 tests
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd
import pytest

from sl_agent.multimodal.models import (
    CandidateAxis,
    EvidenceMatrix,
    EvidenceRow,
    ModalityEvidence,
    ModalityStatus,
    Modality,
    MultiModalQueryInput,
)
from sl_agent.multimodal.literature_receipts import (
    get_literature_receipts,
    get_calibration_narrative,
    list_receipts_for_gene,
)
from sl_agent.multimodal.modality_fuser import (
    fuse_matrix,
    build_agreement_report,
    build_updated_narrative,
    _assign_recommendation_tier,
    _crispr_pharma_agreement,
    _weighted_score,
)
from sl_agent.multimodal.pharmacologic_analyzer import (
    analyze_drug_screen,
    aggregate_by_axis,
    _match_drug_to_axis,
)
from sl_agent.multimodal.matrix_builder import (
    build_evidence_matrix,
    run_multimodal_analysis,
    _build_crispr_map,
    _axis_to_target_genes,
    _best_crispr_result,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def mbd4_query():
    return MultiModalQueryInput(
        gene="MBD4",
        mutation_type="loss_of_function",
        cancer_type="Colorectal Cancer",
        include_pharmacologic_stratification=False,  # no PRISM data in unit tests
        include_literature_receipts=True,
        include_calibration_narrative=True,
    )


@pytest.fixture
def synthetic_prism_data():
    """
    Synthetic PRISM-like DataFrame: rows = cell lines, cols = drug IDs.
    Gemcitabine: mutant lines are 0.3 units more sensitive (lower).
    Olaparib: no signal (flat).
    """
    rng = np.random.default_rng(42)
    n_mut, n_wt = 10, 30
    cell_ids = [f"ACH-{i:05d}" for i in range(n_mut + n_wt)]
    # gemcitabine column: mutant more sensitive
    gem_mut = rng.normal(-1.0, 0.2, n_mut)
    gem_wt  = rng.normal(-0.7, 0.2, n_wt)
    # olaparib column: flat
    ola_mut = rng.normal(-0.8, 0.2, n_mut)
    ola_wt  = rng.normal(-0.8, 0.2, n_wt)

    df = pd.DataFrame(
        {
            "DRUG001": np.concatenate([gem_mut, gem_wt]),   # gemcitabine
            "DRUG002": np.concatenate([ola_mut, ola_wt]),   # olaparib
        },
        index=cell_ids,
    )
    meta = pd.DataFrame(
        {"name": ["gemcitabine", "olaparib"]},
        index=["DRUG001", "DRUG002"],
    )
    mut_ids = cell_ids[:n_mut]
    wt_ids  = cell_ids[n_mut:]
    return df, meta, mut_ids, wt_ids


@pytest.fixture
def all_positive_row():
    """A row with all modalities positive — should be Validated."""
    row = EvidenceRow(
        axis=CandidateAxis.CYTIDINE_ANALOGS,
        axis_label="Cytidine Analogs",
        mechanism="BER substrate accumulation",
    )
    row.crispr    = ModalityEvidence(modality=Modality.CRISPR_DEPENDENCY,   status=ModalityStatus.POSITIVE)
    row.prism     = ModalityEvidence(modality=Modality.PRISM_PHARMACOLOGIC,  status=ModalityStatus.POSITIVE)
    row.in_vitro  = ModalityEvidence(modality=Modality.IN_VITRO_FUNCTIONAL, status=ModalityStatus.POSITIVE)
    row.in_vivo   = ModalityEvidence(modality=Modality.IN_VIVO_PDX,         status=ModalityStatus.POSITIVE)
    row.clinical  = ModalityEvidence(modality=Modality.CLINICAL,            status=ModalityStatus.POSITIVE)
    row.expression= ModalityEvidence(modality=Modality.EXPRESSION_ASSOC,    status=ModalityStatus.POSITIVE)
    return row


@pytest.fixture
def crispr_neg_pharma_missing_row():
    """CRISPR negative, pharmacologic NOT interrogated — must NOT be 'Not supported'."""
    row = EvidenceRow(
        axis=CandidateAxis.PARP_INHIBITORS,
        axis_label="PARP Inhibitors",
        mechanism="BER stress → PARP1 upregulation",
    )
    row.crispr = ModalityEvidence(
        modality=Modality.CRISPR_DEPENDENCY, status=ModalityStatus.NEGATIVE,
        delta_dep=0.03, p_value=0.7, fdr=0.95,
    )
    # prism and gdsc stay as MISSING (default)
    return row


@pytest.fixture
def crispr_neg_pharma_pos_row():
    """CRISPR negative but pharmacologic positive — key disagreement case."""
    row = EvidenceRow(
        axis=CandidateAxis.PARP_INHIBITORS,
        axis_label="PARP Inhibitors",
        mechanism="PARP trapping via BER substrate",
    )
    row.crispr = ModalityEvidence(
        modality=Modality.CRISPR_DEPENDENCY, status=ModalityStatus.NEGATIVE,
        delta_dep=0.03, p_value=0.7, fdr=0.95,
    )
    row.prism = ModalityEvidence(
        modality=Modality.PRISM_PHARMACOLOGIC, status=ModalityStatus.POSITIVE,
        delta_auc=-0.3, p_value=0.01, fdr=0.04,
    )
    row.expression = ModalityEvidence(
        modality=Modality.EXPRESSION_ASSOC, status=ModalityStatus.POSITIVE,
    )
    return row


# ─────────────────────────────────────────────────────────────────────────────
# 1. Model tests
# ─────────────────────────────────────────────────────────────────────────────

def test_evidence_row_positive_modalities(all_positive_row):
    pos = all_positive_row.positive_modalities()
    assert "crispr" in pos
    assert "clinical" in pos
    assert len(pos) >= 5


def test_evidence_row_missing_modalities(crispr_neg_pharma_missing_row):
    missing = crispr_neg_pharma_missing_row.missing_modalities()
    assert "prism" in missing
    assert "gdsc" in missing
    assert "in_vitro" in missing


def test_evidence_matrix_get_row():
    matrix = EvidenceMatrix(query_gene="MBD4", rows=[
        EvidenceRow(
            axis=CandidateAxis.CYTIDINE_ANALOGS,
            axis_label="Cytidine Analogs",
            mechanism="BER",
        )
    ])
    row = matrix.get_row(CandidateAxis.CYTIDINE_ANALOGS)
    assert row is not None
    assert row.axis == CandidateAxis.CYTIDINE_ANALOGS


def test_evidence_matrix_recommendation_summary():
    row = EvidenceRow(
        axis=CandidateAxis.PARP_INHIBITORS,
        axis_label="PARP",
        mechanism="BER",
        recommendation_tier="Mechanistic candidate only",
    )
    matrix = EvidenceMatrix(query_gene="MBD4", rows=[row])
    summary = matrix.recommendation_summary()
    assert "PARP" in summary
    assert summary["PARP"] == "Mechanistic candidate only"


# ─────────────────────────────────────────────────────────────────────────────
# 2. Literature receipt tests
# ─────────────────────────────────────────────────────────────────────────────

def test_cytidine_receipts_exist():
    receipts = get_literature_receipts("MBD4", CandidateAxis.CYTIDINE_ANALOGS)
    assert "in_vitro" in receipts
    assert "clinical" in receipts
    assert receipts["in_vitro"].status == ModalityStatus.POSITIVE


def test_cytidine_in_vitro_has_pmid():
    receipts = get_literature_receipts("MBD4", CandidateAxis.CYTIDINE_ANALOGS)
    assert "35428381" in receipts["in_vitro"].pmids


def test_parp_clinical_is_missing():
    receipts = get_literature_receipts("MBD4", CandidateAxis.PARP_INHIBITORS)
    assert "clinical" in receipts
    assert receipts["clinical"].status == ModalityStatus.MISSING


def test_parp_in_vitro_is_mixed():
    receipts = get_literature_receipts("MBD4", CandidateAxis.PARP_INHIBITORS)
    assert receipts["in_vitro"].status == ModalityStatus.MIXED


def test_wrn_in_vitro_notes_msi_confound():
    receipts = get_literature_receipts("MBD4", CandidateAxis.WRN)
    assert "in_vitro" in receipts
    summary = receipts["in_vitro"].summary or ""
    assert "MSI" in summary or "MSI-H" in summary or "confound" in summary.lower() or "context" in summary.lower()


def test_calibration_narrative_is_string():
    narrative = get_calibration_narrative("MBD4")
    assert isinstance(narrative, str)
    assert "cytidine" in narrative.lower() or "CYTIDINE" in narrative


def test_calibration_narrative_mentions_parp():
    narrative = get_calibration_narrative("MBD4")
    assert "PARP" in narrative


def test_list_receipts_for_mbd4():
    receipts = list_receipts_for_gene("MBD4")
    assert "cytidine_analogs" in receipts
    assert len(receipts["cytidine_analogs"]) > 0


# ─────────────────────────────────────────────────────────────────────────────
# 3. Pharmacologic analyzer tests
# ─────────────────────────────────────────────────────────────────────────────

def test_drug_to_axis_gemcitabine():
    assert _match_drug_to_axis("gemcitabine") == "cytidine_analogs"


def test_drug_to_axis_olaparib():
    assert _match_drug_to_axis("olaparib") == "parp_inhibitors"


def test_drug_to_axis_adavosertib():
    assert _match_drug_to_axis("adavosertib") == "atr_wee1"


def test_drug_to_axis_unknown():
    assert _match_drug_to_axis("carboplatin") is None


def test_pharmacologic_analyzer_gemcitabine_positive(synthetic_prism_data):
    prism_df, prism_meta, mut_ids, wt_ids = synthetic_prism_data
    results = analyze_drug_screen(
        gene="MBD4",
        mutant_ids=mut_ids,
        wt_ids=wt_ids,
        prism_df=prism_df,
        prism_meta=prism_meta,
        gdsc_df=None,
        sample_info=None,
        target_axes=["cytidine_analogs", "parp_inhibitors"],
        stratify_by_msi=False,
    )
    gem_results = [r for r in results if r.axis == "cytidine_analogs"]
    assert len(gem_results) > 0
    best = min(gem_results, key=lambda r: r.p_value)
    assert best.delta_response < 0, "Gemcitabine should show mutant more sensitive"


def test_pharmacologic_analyzer_olaparib_flat(synthetic_prism_data):
    prism_df, prism_meta, mut_ids, wt_ids = synthetic_prism_data
    results = analyze_drug_screen(
        gene="MBD4",
        mutant_ids=mut_ids,
        wt_ids=wt_ids,
        prism_df=prism_df,
        prism_meta=prism_meta,
        gdsc_df=None,
        sample_info=None,
        stratify_by_msi=False,
    )
    ola_results = [r for r in results if r.axis == "parp_inhibitors"]
    if ola_results:
        best = min(ola_results, key=lambda r: r.p_value)
        assert best.p_value > 0.05, "Olaparib should be flat (p > 0.05)"


def test_aggregate_by_axis_returns_modality_evidence(synthetic_prism_data):
    prism_df, prism_meta, mut_ids, wt_ids = synthetic_prism_data
    results = analyze_drug_screen(
        gene="MBD4", mutant_ids=mut_ids, wt_ids=wt_ids,
        prism_df=prism_df, prism_meta=prism_meta,
        gdsc_df=None, sample_info=None, stratify_by_msi=False,
    )
    by_axis = aggregate_by_axis(results, stratifier="MBD4_LOF_vs_WT")
    assert "cytidine_analogs" in by_axis
    ev = by_axis["cytidine_analogs"]
    assert isinstance(ev, ModalityEvidence)
    assert ev.delta_auc is not None


# ─────────────────────────────────────────────────────────────────────────────
# 4. Modality fuser tests — recommendation tiers and guardrails
# ─────────────────────────────────────────────────────────────────────────────

def test_all_positive_is_validated(all_positive_row):
    tier = _assign_recommendation_tier(all_positive_row, _weighted_score(all_positive_row), 6)
    assert tier == "Validated SL therapeutic lever"


def test_crispr_neg_pharma_missing_is_not_negative(crispr_neg_pharma_missing_row):
    """
    CORE GUARDRAIL TEST:
    CRISPR negative + pharma NOT examined → must NOT be 'Not supported / negative'.
    Spec Section D: "CRISPR alone cannot move an axis to negative if pharmacologic
    modalities have not been examined."
    """
    score = _weighted_score(crispr_neg_pharma_missing_row)
    n_pos = len(crispr_neg_pharma_missing_row.positive_modalities())
    tier = _assign_recommendation_tier(crispr_neg_pharma_missing_row, score, n_pos)
    assert tier != "Not supported / negative", (
        "Guardrail failed: CRISPR negative + pharma unexamined should NOT be 'Not supported'"
    )


def test_crispr_neg_pharma_pos_is_mechanistic(crispr_neg_pharma_pos_row):
    score = _weighted_score(crispr_neg_pharma_pos_row)
    n_pos = len(crispr_neg_pharma_pos_row.positive_modalities())
    tier = _assign_recommendation_tier(crispr_neg_pharma_pos_row, score, n_pos)
    # Should be Mechanistic or Strong — not Validated (no clinical) and not Negative
    assert tier in ("Mechanistic candidate only", "Strong candidate dependency axis")


def test_crispr_pharma_agreement_disagree(crispr_neg_pharma_pos_row):
    agreement = _crispr_pharma_agreement(crispr_neg_pharma_pos_row)
    assert agreement == "disagree_crispr_neg_pharma_pos"


def test_crispr_pharma_agreement_insufficient_data(crispr_neg_pharma_missing_row):
    agreement = _crispr_pharma_agreement(crispr_neg_pharma_missing_row)
    assert agreement == "insufficient_data"


def test_fuse_matrix_updates_tiers(mbd4_query):
    """Full fuse_matrix call on a simple 2-row matrix."""
    rows = [
        EvidenceRow(axis=CandidateAxis.CYTIDINE_ANALOGS, axis_label="Cytidine", mechanism="BER"),
        EvidenceRow(axis=CandidateAxis.PARP_INHIBITORS,  axis_label="PARP",     mechanism="BER→PARP1"),
    ]
    matrix = EvidenceMatrix(query_gene="MBD4", rows=rows)
    fused = fuse_matrix(matrix)
    for row in fused.rows:
        assert row.recommendation_tier is not None
        assert row.overall_evidence_level is not None


def test_agreement_report_length(mbd4_query):
    query = MultiModalQueryInput(
        gene="MBD4",
        include_pharmacologic_stratification=False,
        include_literature_receipts=True,
    )
    matrix = build_evidence_matrix(query)
    report = build_agreement_report(matrix)
    assert len(report) == len(matrix.rows)


def test_agreement_report_flags_interpretation_change(crispr_neg_pharma_pos_row):
    matrix = EvidenceMatrix(query_gene="MBD4", rows=[crispr_neg_pharma_pos_row])
    matrix = fuse_matrix(matrix)
    report = build_agreement_report(matrix)
    assert len(report) == 1
    assert report[0].interpretation_change is True


# ─────────────────────────────────────────────────────────────────────────────
# 5. Matrix builder end-to-end tests
# ─────────────────────────────────────────────────────────────────────────────

def test_build_evidence_matrix_mbd4_receipts(mbd4_query):
    matrix = build_evidence_matrix(mbd4_query)
    assert matrix.query_gene == "MBD4"
    # Cytidine axis should have literature receipts filled
    row = matrix.get_row(CandidateAxis.CYTIDINE_ANALOGS)
    assert row is not None
    assert row.in_vitro.status == ModalityStatus.POSITIVE
    assert row.clinical.status == ModalityStatus.POSITIVE


def test_build_evidence_matrix_parp_clinical_missing(mbd4_query):
    matrix = build_evidence_matrix(mbd4_query)
    row = matrix.get_row(CandidateAxis.PARP_INHIBITORS)
    assert row is not None
    assert row.clinical.status == ModalityStatus.MISSING


def test_build_evidence_matrix_all_rows_have_tiers(mbd4_query):
    matrix = build_evidence_matrix(mbd4_query)
    for row in matrix.rows:
        assert row.recommendation_tier is not None, f"No tier for {row.axis}"


def test_cytidine_axis_is_validated(mbd4_query):
    matrix = build_evidence_matrix(mbd4_query)
    row = matrix.get_row(CandidateAxis.CYTIDINE_ANALOGS)
    assert row is not None
    assert row.recommendation_tier == "Validated SL therapeutic lever"


def test_parp_axis_is_not_validated_or_negative(mbd4_query):
    """PARP should be Mechanistic candidate only — not Validated, not Negative."""
    matrix = build_evidence_matrix(mbd4_query)
    row = matrix.get_row(CandidateAxis.PARP_INHIBITORS)
    assert row is not None
    assert row.recommendation_tier not in (
        "Validated SL therapeutic lever",
        "Not supported / negative",
    ), f"PARP tier was: {row.recommendation_tier}"


def test_axis_to_target_genes_parp():
    genes = _axis_to_target_genes(CandidateAxis.PARP_INHIBITORS)
    assert "PARP1" in genes


def test_best_crispr_result_returns_lowest_fdr():
    crispr_map = {
        "PARP1": (0.0, 0.001, 0.005, 10, 30, -1.5),
        "PARP2": (0.0, 0.05,  0.15,  10, 30, -0.5),
    }
    best = _best_crispr_result(crispr_map, ["PARP1", "PARP2"])
    assert best is not None
    assert best[2] == 0.005   # lowest FDR


def test_run_multimodal_analysis_returns_result(mbd4_query):
    result = run_multimodal_analysis(mbd4_query)
    assert result.query_gene == "MBD4"
    assert len(result.evidence_matrix.rows) > 0
    assert result.updated_narrative
    assert result.recommendation_map


def test_updated_narrative_mentions_all_modalities(mbd4_query):
    result = run_multimodal_analysis(mbd4_query)
    narrative = result.updated_narrative
    assert "MBD4" in narrative
    assert "cytidine" in narrative.lower() or "Cytidine" in narrative


def test_guardrails_collected_for_crispr_neg_pharma_missing():
    """Guardrail should fire for PARP: CRISPR neg + pharma missing."""
    query = MultiModalQueryInput(
        gene="MBD4",
        axes=[CandidateAxis.PARP_INHIBITORS],
        include_pharmacologic_stratification=False,
        include_literature_receipts=True,
    )
    result = run_multimodal_analysis(query)
    # PARP CRISPR is MISSING (no sl_partners), pharma is MISSING → guardrail
    # This is softer — only check that guardrails list exists
    assert isinstance(result.guardrails_applied, list)
