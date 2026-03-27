"""
Unit tests for drug mapper layer (offline — mocks ChEMBL and OncoKB).
"""
import numpy as np
import pandas as pd
import pytest

from sl_agent.core.drug_mapper import (
    compute_rank_score,
    prism_drug_sensitivity,
    _infer_drug_class,
)
from sl_agent.core.models import DrugClass


# ── Fixtures ──────────────────────────────────────────────────────────────────

N_LINES = 60
N_DRUGS = 20
RNG = np.random.default_rng(0)


def _make_prism(n_lines=N_LINES, n_drugs=N_DRUGS) -> pd.DataFrame:
    model_ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    drugs = [f"BRD-K{i:08d}" for i in range(n_drugs)]
    data = RNG.normal(-1.0, 0.5, (n_lines, n_drugs))
    # Inject drug 0 as strongly selective for first 15 lines (mutant)
    data[:15, 0] -= 2.0
    return pd.DataFrame(data, index=model_ids, columns=drugs)


def _make_prism_meta(n_drugs=N_DRUGS) -> pd.DataFrame:
    drugs = [f"BRD-K{i:08d}" for i in range(n_drugs)]
    names = [f"Drug_{i}" for i in range(n_drugs)]
    return pd.DataFrame({"name": names}, index=drugs)


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_prism_sensitivity_returns_evidence():
    prism = _make_prism()
    meta = _make_prism_meta()
    mutant_ids = [f"ACH-{i:06d}" for i in range(15)]
    wt_ids = [f"ACH-{i:06d}" for i in range(15, 50)]

    evidence = prism_drug_sensitivity(
        gene="BRCA1",
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
        prism_df=prism,
        prism_meta=meta,
    )
    assert len(evidence) > 0
    # All evidence should be selective (mutant more sensitive = delta < 0)
    for ev in evidence:
        assert ev.delta_viability < 0


def test_prism_drug0_is_top_hit():
    prism = _make_prism()
    meta = _make_prism_meta()
    mutant_ids = [f"ACH-{i:06d}" for i in range(15)]
    wt_ids = [f"ACH-{i:06d}" for i in range(15, 50)]

    evidence = prism_drug_sensitivity(
        gene="BRCA1",
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
        prism_df=prism,
        prism_meta=meta,
    )
    # Drug_0 should be the top hit by delta
    top_name = evidence[0].drug_name
    assert top_name == "Drug_0", f"Expected Drug_0 as top hit, got {top_name}"


def test_prism_insufficient_lines_returns_empty():
    prism = _make_prism()
    meta = _make_prism_meta()
    evidence = prism_drug_sensitivity(
        gene="X",
        mutant_ids=["ACH-000001"],  # < 3
        wt_ids=["ACH-000002"],
        prism_df=prism,
        prism_meta=meta,
    )
    assert evidence == []


def test_rank_score_range():
    score, components = compute_rank_score(
        sl_delta=-1.2,
        sl_fdr=0.01,
        drug_delta_viability=-2.0,
        drug_fdr=0.05,
        max_phase=4,
        oncokb_level="LEVEL_1",
    )
    assert 0.0 <= score <= 1.0
    assert "sl_signal" in components
    assert "drug_response" in components
    assert "druggability" in components


def test_rank_score_higher_when_stronger_evidence():
    strong, _ = compute_rank_score(
        sl_delta=-1.5, sl_fdr=0.01,
        drug_delta_viability=-2.0, drug_fdr=0.01,
        max_phase=4, oncokb_level="LEVEL_1",
    )
    weak, _ = compute_rank_score(
        sl_delta=-0.2, sl_fdr=0.20,
        drug_delta_viability=None, drug_fdr=None,
        max_phase=0, oncokb_level=None,
    )
    assert strong > weak


def test_drug_class_inference():
    assert _infer_drug_class("PARP trapping inhibitor") == DrugClass.PARP_INHIBITOR
    assert _infer_drug_class("ATR kinase inhibitor") == DrugClass.ATR_INHIBITOR
    assert _infer_drug_class("WEE1 inhibitor") == DrugClass.WEE1_INHIBITOR
    assert _infer_drug_class("HDAC inhibitor") == DrugClass.EPIGENETIC
    assert _infer_drug_class("PD-L1 checkpoint") == DrugClass.IMMUNOTHERAPY
    assert _infer_drug_class("novel compound") == DrugClass.OTHER
    assert _infer_drug_class("") == DrugClass.UNKNOWN
