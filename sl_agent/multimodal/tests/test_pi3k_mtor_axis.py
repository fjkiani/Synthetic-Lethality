"""
Tests for the first-class pi3k_mtor axis promotion.

Covers:
  1. Schema: pi3k_mtor exists; wrn/pkmyt1/pi3k_akt preserved (no axis deleted).
  2. Drug routing: PI3K/AKT/mTOR drugs resolve to pi3k_mtor; other axes' drug
     routing unchanged (no leakage / no suppression). Metformin NOT mapped.
  3. Grounding guard: the pi3k_mtor interpretation carries ONLY validated numbers
     from the VALIDATION_SUMMARY receipt and NONE of the request's fabricated
     strings (HR=1.30, p=0.023, p=0.043, GSE160626, "prognostic killer",
     "toward platinum sensitivity").
  4. Interpretation survives the full fuse pipeline unchanged.
"""
import re
import pytest

from sl_agent.multimodal.models import CandidateAxis, EvidenceRow
from sl_agent.multimodal.matrix_builder import _AXIS_META
from sl_agent.multimodal.modality_fuser import _build_interpretation
from sl_agent.multimodal.pharmacologic_analyzer import _DRUG_TO_AXIS


# ── 1. Schema ────────────────────────────────────────────────────────────────
def test_pi3k_mtor_axis_exists():
    assert CandidateAxis("pi3k_mtor") == CandidateAxis.PI3K_MTOR


@pytest.mark.parametrize("axis_value", ["wrn", "pkmyt1", "pi3k_akt"])
def test_preexisting_axes_not_deleted(axis_value):
    # The request called wrn/pkmyt1 "deprecated" — they are ACTIVE and must remain.
    assert CandidateAxis(axis_value).value == axis_value


# ── 2. Drug routing ──────────────────────────────────────────────────────────
@pytest.mark.parametrize("drug", [
    "everolimus", "temsirolimus", "sirolimus", "rapamycin",
    "alpelisib", "byl719", "byl-719",
    "capivasertib", "azd5363", "inavolisib", "ipatasertib",
    "buparlisib", "bkm120", "pictilisib", "copanlisib",
])
def test_pi3k_mtor_drugs_route_to_axis(drug):
    assert _DRUG_TO_AXIS[drug] == "pi3k_mtor"


def test_metformin_not_mapped():
    # No grounded receipt for metformin -> must NOT be silently invented.
    assert "metformin" not in _DRUG_TO_AXIS


@pytest.mark.parametrize("drug,expected", [
    ("olaparib", "parp_inhibitors"),
    ("rp-6306", "pkmyt1"),
    ("mrtx1719", "wrn"),
    ("bortezomib", "proteasome"),
    ("venetoclax", "bcl2_mcl1"),
])
def test_other_axes_routing_unchanged(drug, expected):
    # No leakage: adding pi3k_mtor must not perturb existing drug->axis routing.
    assert _DRUG_TO_AXIS[drug] == expected


# ── 3. Grounding guard ───────────────────────────────────────────────────────
FABRICATED_STRINGS = [
    "1.30", "HR=1.30", "0.023", "p=0.023",
    "0.043", "p=0.043", "GSE160626",
    "prognostic killer", "toward platinum sensitivity",
]

VALIDATED_TOKENS = [
    "NOT prognostic", "n=299", "n=303", "p=0.88", "p=0.93",
    "RESEARCH USE ONLY",
]


def test_pi3k_mtor_interpretation_has_no_fabricated_claims():
    interp = _AXIS_META[CandidateAxis.PI3K_MTOR]["interpretation"]
    for bad in FABRICATED_STRINGS:
        assert bad not in interp, f"Fabricated string leaked into interpretation: {bad!r}"


def test_pi3k_mtor_interpretation_is_validation_grounded():
    interp = _AXIS_META[CandidateAxis.PI3K_MTOR]["interpretation"]
    for tok in VALIDATED_TOKENS:
        assert tok in interp, f"Expected validated token missing: {tok!r}"


def test_pi3k_mtor_has_evidence_provenance():
    assert "evidence_provenance" in _AXIS_META[CandidateAxis.PI3K_MTOR]


# ── 4. Interpretation survives the fuser ─────────────────────────────────────
def test_validated_interpretation_survives_fuse():
    row = EvidenceRow(
        axis=CandidateAxis.PI3K_MTOR,
        axis_label=_AXIS_META[CandidateAxis.PI3K_MTOR]["label"],
        mechanism=_AXIS_META[CandidateAxis.PI3K_MTOR]["mechanism"],
    )
    out = _build_interpretation(row)
    assert out == _AXIS_META[CandidateAxis.PI3K_MTOR]["interpretation"]
    # and it still contains no fabricated claim after passing through the fuser
    for bad in FABRICATED_STRINGS:
        assert bad not in out
