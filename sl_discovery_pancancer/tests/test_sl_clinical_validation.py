"""
Unit tests for sl_clinical_validation v2.0 clinical-grade layer.

Covers the engine decision functions with positive + boundary cases:
  - replicate_axis: sign/status on known synthetic input, edge cases
    (insufficient_n, partner_absent, CN not-testable)
  - independent_replication: BH-FDR application + status taxonomy
  - assign_evidence_grade: precedent_class -> E1..E5 mapping
  - assign_trust_tier_v2: (grade x replication) -> tier
  - phase_rank ordering
  - load_sanger_score gene-column cleaning

Run: pytest sl_discovery_pancancer/tests/ -q
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import sl_clinical_validation as scv
from sl_battery import MIN_N


# ---------------------------------------------------------------------------
# Fixtures: a tiny synthetic screen + mutation table with a KNOWN SL signal
# ---------------------------------------------------------------------------
@pytest.fixture
def synthetic_screen():
    """12 cell lines. DRV-mutant lines (5) are strongly dependent on PTNR
    (effect ~ -1.5) vs WT (~0). One CN axis and one partner-absent case tested
    via separate helpers.
    """
    lines = [f"ACH-{i:06d}" for i in range(12)]
    # PTNR: first 5 lines (mutant) very negative, rest near 0
    ptnr = [-1.6, -1.5, -1.7, -1.4, -1.55, -0.05, 0.0, 0.1, -0.1, 0.05, -0.02, 0.03]
    # NEUTRAL gene: no difference
    neutral = [-0.1, 0.0, 0.05, -0.05, 0.02, -0.03, 0.01, 0.0, -0.02, 0.04, -0.01, 0.02]
    df = pd.DataFrame({"PTNR": ptnr, "NEUTRAL": neutral}, index=lines)
    return df


@pytest.fixture
def synthetic_muts():
    """DRV LOF in first 5 lines; ACTdrv missense in lines 0-4 too."""
    lines = [f"ACH-{i:06d}" for i in range(12)]
    recs = []
    for i, ln in enumerate(lines):
        if i < 5:
            recs.append(dict(gene_symbol="DRV", ModelID=ln, effect="frameshift",
                             cancer_driver=True))
    # a driver with too FEW mutants (only 2) -> insufficient_n
    for i, ln in enumerate(lines):
        if i < 2:
            recs.append(dict(gene_symbol="RAREDRV", ModelID=ln, effect="nonsense",
                             cancer_driver=True))
    return pd.DataFrame(recs)


# ---------------------------------------------------------------------------
# replicate_axis
# ---------------------------------------------------------------------------
def test_replicate_axis_positive_signal(synthetic_screen, synthetic_muts):
    """Known SL: DRV-LOF lines depleted for PTNR -> tested, negative delta/d, small p."""
    res = scv.replicate_axis(synthetic_screen, synthetic_muts, "DRV", "PTNR", "LOF")
    assert res["screen_status"] == "tested"
    assert res["screen_n_mut"] == 5
    assert res["screen_n_wt"] == 7
    assert res["screen_delta"] < 0          # partner more depleted in mutant
    assert res["screen_d"] < 0              # negative Cohen's d
    assert res["screen_p"] < 0.05           # significant one-sided


def test_replicate_axis_no_signal(synthetic_screen, synthetic_muts):
    """Neutral gene -> tested but not significant, near-zero effect."""
    res = scv.replicate_axis(synthetic_screen, synthetic_muts, "DRV", "NEUTRAL", "LOF")
    assert res["screen_status"] == "tested"
    assert res["screen_p"] > 0.05


def test_replicate_axis_insufficient_n(synthetic_screen, synthetic_muts):
    """RAREDRV has only 2 mutant lines (< MIN_N) -> insufficient_n."""
    res = scv.replicate_axis(synthetic_screen, synthetic_muts, "RAREDRV", "PTNR", "LOF")
    assert res["screen_status"] == "insufficient_n"
    assert res["screen_n_mut"] < MIN_N


def test_replicate_axis_partner_absent(synthetic_screen, synthetic_muts):
    """Partner gene not in matrix -> partner_absent."""
    res = scv.replicate_axis(synthetic_screen, synthetic_muts, "DRV", "NOT_A_GENE", "LOF")
    assert res["screen_status"] == "partner_absent"


def test_replicate_axis_cn_not_testable(synthetic_screen, synthetic_muts):
    """CN modes are not re-derivable from this matrix -> na_cn."""
    res = scv.replicate_axis(synthetic_screen, synthetic_muts, "DRV", "PTNR", "CN_loss")
    assert res["screen_status"] == "na_cn"


# ---------------------------------------------------------------------------
# independent_replication (BH-FDR + status taxonomy)
# ---------------------------------------------------------------------------
def test_independent_replication_taxonomy(synthetic_screen, synthetic_muts):
    axes = pd.DataFrame([
        dict(driver="DRV", partner="PTNR", mode="LOF", L1_cohens_d=-1.2),   # replicate
        dict(driver="DRV", partner="NEUTRAL", mode="LOF", L1_cohens_d=-0.9),  # concordant/discordant not sig
        dict(driver="RAREDRV", partner="PTNR", mode="LOF", L1_cohens_d=-1.0),  # insufficient_n
        dict(driver="DRV", partner="PTNR", mode="CN_loss", L1_cohens_d=-1.0),  # not_testable_cn
    ])
    out = scv.independent_replication(synthetic_screen, synthetic_muts, axes)
    key = out["driver"] + "|" + out["partner"] + "|" + out["mode"]
    status = dict(zip(key, out["independent_replication_status"]))
    assert status["DRV|PTNR|LOF"] == "replicated_crossplatform"
    assert status["DRV|PTNR|CN_loss"] == "not_testable_cn"
    assert status["RAREDRV|PTNR|LOF"] == "insufficient_n"
    assert status["DRV|NEUTRAL|LOF"] != "replicated_crossplatform"
    # FDR only populated for tested rows
    tested = out["screen_status"] == "tested"
    assert out.loc[tested, "sanger_fdr"].notna().all()
    assert out.loc[~tested, "sanger_fdr"].isna().all()


# ---------------------------------------------------------------------------
# assign_evidence_grade
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pc,expected", [
    ("clinical_target_in_trials", "E2"),
    ("established_SL_drug_dev", "E3"),
    ("preclinical_in_vivo", "E3"),
    ("preclinical_in_vitro", "E4"),
    ("mechanistic_no_drug", "E4"),
    ("none_found", "E5"),
    ("unknown_class", "E5"),   # conservative fallback
])
def test_assign_evidence_grade(pc, expected):
    assert scv.assign_evidence_grade(pc) == expected


# ---------------------------------------------------------------------------
# assign_trust_tier_v2 (2D posture)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("eg,rs,expected", [
    ("E1", "replicated_crossplatform", "T1_anchored_and_replicated"),
    ("E2", "replicated_crossplatform", "T1_anchored_and_replicated"),
    ("E2", "concordant_not_sig", "T2_anchored_precedent"),
    ("E1", "insufficient_n", "T2_anchored_precedent"),
    ("E3", "replicated_crossplatform", "T2_anchored_precedent"),
    ("E4", "replicated_crossplatform", "T2_anchored_precedent"),
    ("E4", "discordant", "T3_precedent_only"),
    ("E5", "replicated_crossplatform", "T3_replicated_no_precedent"),
    ("E5", "discordant", "T4_discovery_only"),
    ("E5", "not_testable_cn", "T4_discovery_only"),
])
def test_assign_trust_tier_v2(eg, rs, expected):
    assert scv.assign_trust_tier_v2(eg, rs) == expected


# ---------------------------------------------------------------------------
# phase_rank ordering
# ---------------------------------------------------------------------------
def test_phase_rank_monotonic():
    assert scv.phase_rank("Launched") > scv.phase_rank("Phase 3")
    assert scv.phase_rank("Phase 3") > scv.phase_rank("Phase 1")
    assert scv.phase_rank("Preclinical") >= 0


# ---------------------------------------------------------------------------
# load_sanger_score gene-column cleaning (uses a tiny temp CSV)
# ---------------------------------------------------------------------------
def test_load_sanger_score_cleans_columns(tmp_path):
    p = tmp_path / "mini.csv"
    pd.DataFrame(
        {"A1BG (1)": [-0.1, -0.2], "TP53 (7157)": [-1.0, -1.1]},
        index=["ACH-000001", "ACH-000002"],
    ).to_csv(p)
    df = scv.load_sanger_score(str(p))
    assert list(df.columns) == ["A1BG", "TP53"]   # Entrez suffix stripped
    assert list(df.index) == ["ACH-000001", "ACH-000002"]
