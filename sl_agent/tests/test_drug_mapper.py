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


# ── MM Sprint 1: drug-to-axis mapping tests ───────────────────────────────────

def test_mm_drug_to_axis_proteasome():
    """Bortezomib and carfilzomib must map to proteasome axis."""
    from sl_agent.multimodal.pharmacologic_analyzer import _match_drug_to_axis
    assert _match_drug_to_axis("bortezomib") == "proteasome"
    assert _match_drug_to_axis("carfilzomib") == "proteasome"
    assert _match_drug_to_axis("ixazomib") == "proteasome"
    assert _match_drug_to_axis("marizomib") == "proteasome"


def test_mm_drug_to_axis_bcl2_mcl1():
    """Venetoclax and navitoclax must map to bcl2_mcl1 axis."""
    from sl_agent.multimodal.pharmacologic_analyzer import _match_drug_to_axis
    assert _match_drug_to_axis("venetoclax") == "bcl2_mcl1"
    assert _match_drug_to_axis("navitoclax") == "bcl2_mcl1"
    assert _match_drug_to_axis("s63845") == "bcl2_mcl1"


def test_mm_drug_to_axis_hdac():
    """Panobinostat and vorinostat must map to hdac axis."""
    from sl_agent.multimodal.pharmacologic_analyzer import _match_drug_to_axis
    assert _match_drug_to_axis("panobinostat") == "hdac"
    assert _match_drug_to_axis("vorinostat") == "hdac"
    assert _match_drug_to_axis("romidepsin") == "hdac"
    assert _match_drug_to_axis("entinostat") == "hdac"


def test_mm_drug_to_axis_imid():
    """Lenalidomide must map to imid axis."""
    from sl_agent.multimodal.pharmacologic_analyzer import _match_drug_to_axis
    assert _match_drug_to_axis("lenalidomide") == "imid"
    assert _match_drug_to_axis("pomalidomide") == "imid"
    assert _match_drug_to_axis("thalidomide") == "imid"


def test_mm_drug_to_axis_no_false_positives():
    """Daratumumab and melphalan must NOT map to any MM axis (not in GDSC)."""
    from sl_agent.multimodal.pharmacologic_analyzer import _match_drug_to_axis
    assert _match_drug_to_axis("daratumumab") is None
    assert _match_drug_to_axis("melphalan") is None
    assert _match_drug_to_axis("dexamethasone") is None


def test_load_mm_panel_returns_expected_genes():
    """load_mm_panel() must return all 29 panel genes (FAM46C included even if absent from CRISPR)."""
    from sl_agent.multimodal.receipt_miner import load_mm_panel
    panel = load_mm_panel()
    assert len(panel) == 29, f"Expected 29 genes, got {len(panel)}: {panel}"
    # Key genes must be present
    for gene in ["NSD2", "CCND2", "MAF", "BCL2", "MCL1", "IKZF1", "IKZF3", "BRD4"]:
        assert gene in panel, f"{gene} missing from MM panel"


def test_mm_candidate_axes_in_enum():
    """All 4 MM Sprint 1 axes must be valid CandidateAxis enum values."""
    from sl_agent.multimodal.models import CandidateAxis
    assert CandidateAxis.PROTEASOME.value == "proteasome"
    assert CandidateAxis.BCL2_MCL1.value == "bcl2_mcl1"
    assert CandidateAxis.HDAC.value == "hdac"
    assert CandidateAxis.IMID.value == "imid"
