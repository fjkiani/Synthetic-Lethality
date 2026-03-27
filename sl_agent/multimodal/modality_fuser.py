"""
Modality Fuser — reconciles CRISPR dependency with pharmacologic and literature
evidence to produce overall evidence levels and recommendation tiers.

Doctrine:
  1. CRISPR alone cannot move an axis from "mechanistic candidate" to "negative"
     if pharmacologic and clinical modalities have not been examined.
  2. Compound-level sensitivity in the absence of CRISPR dependency = "compound-
     level vulnerability without genetic dependency" (synthetic sickness, not
     necessarily lethality).
  3. The cytidine-analog axis is the gold-standard pattern; everything is graded
     relative to it.

Recommendation tiers (from the spec):
  - "Validated SL therapeutic lever":  multiple modalities positive
  - "Strong candidate dependency axis": CRISPR + compound screens + mechanism
  - "Mechanistic candidate only":       expression + mechanism, no CRISPR/compound
  - "Not supported / negative":         flat or conflicting, or MSI confounded
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

from .models import (
    CandidateAxis,
    EvidenceMatrix,
    EvidenceRow,
    ModalityAgreementReport,
    ModalityEvidence,
    ModalityStatus,
)
from .replication_stress import RSScore, rs_is_high, rs_is_moderate_or_higher

# RS-sensitive axes: ONLY these two can be modulated by RS score
_RS_SENSITIVE_AXES = {CandidateAxis.ATR_WEE1, CandidateAxis.PKMYT1}

logger = logging.getLogger(__name__)

# ── Evidence-level scoring ────────────────────────────────────────────────────

# Weight map: how much does each modality contribute to the overall score?
_MODALITY_WEIGHTS = {
    "clinical":    4,
    "in_vivo":     3,
    "in_vitro":    2,
    "prism":       2,
    "gdsc":        2,
    "expression":  1,
    "crispr":      1,   # intentionally not dominant
}

_STATUS_SCORES = {
    ModalityStatus.POSITIVE:   1.0,
    ModalityStatus.MIXED:      0.4,
    ModalityStatus.NEGATIVE:  -0.5,
    ModalityStatus.CONFOUNDED: 0.0,
    ModalityStatus.MISSING:    0.0,
}


def _weighted_score(row: EvidenceRow) -> float:
    total = 0.0
    for mod_name, weight in _MODALITY_WEIGHTS.items():
        cell = row.cells().get(mod_name)
        if cell:
            total += weight * _STATUS_SCORES.get(cell.status, 0.0)
    return total


def _assign_evidence_level(score: float, n_positive: int) -> str:
    if n_positive >= 3 and score >= 6:
        return "High"
    elif n_positive >= 2 and score >= 3:
        return "Moderate"
    elif score > 0:
        return "Mechanistic-only"
    else:
        return "Negative"


def _assign_recommendation_tier(
    row: EvidenceRow,
    score: float,
    n_positive: int,
    rs_score: Optional[RSScore] = None,
) -> str:
    """
    Apply the 4-tier recommendation logic from the spec guardrails.

    RS modulation is applied AFTER the base logic for ATR/WEE1 and PKMYT1 only.
    Constraints:
      - RS can NEVER promote directly to Validated.
      - RS can only promote Mechanistic → Strong, and ONLY when CRISPR or pharma
        evidence already exists. RS alone is not sufficient for promotion.
      - MBD4 alone gives RS=low (0.05). Cannot promote.
    """
    cells = row.cells()

    # CRISPR status
    crispr_pos  = cells["crispr"].status == ModalityStatus.POSITIVE
    crispr_neg  = cells["crispr"].status == ModalityStatus.NEGATIVE
    crispr_miss = cells["crispr"].status == ModalityStatus.MISSING

    # Pharmacologic status
    pharma_pos = (
        cells["prism"].status == ModalityStatus.POSITIVE
        or cells["gdsc"].status == ModalityStatus.POSITIVE
    )
    pharma_miss = (
        cells["prism"].status == ModalityStatus.MISSING
        and cells["gdsc"].status == ModalityStatus.MISSING
    )

    # Clinical/in-vivo strong signal
    strong_clinical = (
        cells["clinical"].status == ModalityStatus.POSITIVE
        or cells["in_vivo"].status == ModalityStatus.POSITIVE
    )

    # ── Base tier logic (UNCHANGED from v4.0.0) ───────────────────────────────

    # "Validated SL therapeutic lever" — multiple modalities, including clinical/in-vivo
    if strong_clinical and n_positive >= 3:
        return "Validated SL therapeutic lever"

    # "Strong candidate dependency axis" — CRISPR + pharma + mechanism (no clinical req'd)
    if (crispr_pos or pharma_pos) and n_positive >= 2 and not crispr_neg:
        return "Strong candidate dependency axis"

    # Pharmacologic only (CRISPR neg/missing but pharma positive + mechanism)
    if pharma_pos and (crispr_neg or crispr_miss) and cells["expression"].status in (
        ModalityStatus.POSITIVE, ModalityStatus.MIXED
    ):
        base_tier = "Mechanistic candidate only"
    elif (
        cells["expression"].status == ModalityStatus.POSITIVE
        and not pharma_pos
        and not crispr_pos
    ):
        # Expression / pathway only
        base_tier = "Mechanistic candidate only"
    elif crispr_neg and pharma_miss:
        # Guardrail: CRISPR negative BUT we haven't checked pharma — can't call "negative"
        base_tier = "Mechanistic candidate only"  # NOT "not supported" — pharma unexamined
    elif n_positive == 0 and score <= 0:
        base_tier = "Not supported / negative"
    else:
        base_tier = "Mechanistic candidate only"

    # ── RS modulation (ADDITIVE — runs AFTER base logic) ─────────────────────
    # Applies ONLY to ATR/WEE1 and PKMYT1 axes.
    # Promotes Mechanistic → Strong only when:
    #   (a) RS level is "high"
    #   (b) CRISPR or pharma evidence already exists (RS is co-factor, not standalone)
    #   (c) Base tier was Mechanistic (never promotes from lower tiers or overrides Validated)
    if (
        rs_score is not None
        and row.axis in _RS_SENSITIVE_AXES
        and base_tier == "Mechanistic candidate only"
    ):
        has_independent_evidence = crispr_pos or pharma_pos
        if rs_is_high(rs_score) and has_independent_evidence:
            # High RS + independent evidence → promote to Strong
            return "Strong candidate dependency axis"
        # rs_moderate: no tier change (noted in interpretation only)
        # rs_low: no tier change
        # rs_high but no independent evidence: no tier change (RS alone is insufficient)

    return base_tier


def _crispr_pharma_agreement(row: EvidenceRow) -> str:
    cells = row.cells()
    crispr_status = cells["crispr"].status
    pharma_positive = (
        cells["prism"].status == ModalityStatus.POSITIVE
        or cells["gdsc"].status == ModalityStatus.POSITIVE
    )
    pharma_negative = (
        cells["prism"].status == ModalityStatus.NEGATIVE
        and cells["gdsc"].status == ModalityStatus.NEGATIVE
    )
    pharma_missing = (
        cells["prism"].status == ModalityStatus.MISSING
        and cells["gdsc"].status == ModalityStatus.MISSING
    )

    if pharma_missing or crispr_status == ModalityStatus.MISSING:
        return "insufficient_data"
    if crispr_status == ModalityStatus.POSITIVE and pharma_positive:
        return "agree_positive"
    if crispr_status == ModalityStatus.NEGATIVE and pharma_negative:
        return "agree_negative"
    if crispr_status == ModalityStatus.NEGATIVE and pharma_positive:
        return "disagree_crispr_neg_pharma_pos"
    if crispr_status == ModalityStatus.POSITIVE and pharma_negative:
        return "disagree_crispr_pos_pharma_neg"
    return "insufficient_data"


def fuse_matrix(
    matrix: EvidenceMatrix,
    rs_score: Optional[RSScore] = None,
) -> EvidenceMatrix:
    """
    Compute overall_evidence_level, recommendation_tier, and
    crispr_pharmacologic_agreement for every row in the matrix.

    rs_score: Optional RSScore. When provided, used to modulate tier assignment
    for ATR/WEE1 and PKMYT1 axes only (Mechanistic → Strong with high RS + evidence).
    Returns the updated matrix (in-place mutation + return).
    """
    for row in matrix.rows:
        score = _weighted_score(row)
        n_pos = len(row.positive_modalities())

        row.overall_evidence_level = _assign_evidence_level(score, n_pos)
        row.recommendation_tier = _assign_recommendation_tier(row, score, n_pos, rs_score)
        row.crispr_pharmacologic_agreement = _crispr_pharma_agreement(row)
        row.interpretation = _build_interpretation(row, rs_score)

    return matrix


def _build_interpretation(row: EvidenceRow, rs_score: Optional[RSScore] = None) -> str:
    """
    Generate an interpretation statement for one axis.
    Follows the exact doc templates from Section D of the spec.
    """
    cells = row.cells()
    crispr_cell = cells["crispr"]
    crispr_neg = crispr_cell.status == ModalityStatus.NEGATIVE
    crispr_pos = crispr_cell.status == ModalityStatus.POSITIVE
    pharma_pos = (
        cells["prism"].status == ModalityStatus.POSITIVE
        or cells["gdsc"].status == ModalityStatus.POSITIVE
    )
    pharma_missing = (
        cells["prism"].status == ModalityStatus.MISSING
        and cells["gdsc"].status == ModalityStatus.MISSING
    )
    tier = row.recommendation_tier or "Unknown"

    # Spec Section D — Case 1: CRISPR neg AND pharma flat/missing
    if crispr_neg and not pharma_pos and tier == "Not supported / negative":
        return (
            f"No robust SL signal detected for {row.axis_label}. CRISPR dependency is "
            f"flat (Δdep={crispr_cell.delta_dep:+.2f} "
            if crispr_cell.delta_dep is not None
            else f"No robust SL signal detected for {row.axis_label}. CRISPR dependency is "
            f"negative" +
            " and pharmacologic screens show no consistent sensitivity shift. "
            "This axis is currently unsupported in MBD4 LOF context beyond cytidine analogs."
        )

    # Spec Section D — Case 2: CRISPR neg BUT pharma positive
    if crispr_neg and pharma_pos:
        return (
            f"No CRISPR genetic dependency detected for {row.axis_label}, but a compound-level "
            "sensitivity signal is present in pharmacologic screens. This may represent "
            "non-essential target pharmacology (synthetic sickness, trapping, or partial "
            "inhibition) rather than strict genetic synthetic lethality. Classified as "
            "'Mechanistic candidate only' — not a full SL axis without functional validation."
        )

    # CRISPR missing, pharma also not run
    if crispr_cell.status == ModalityStatus.MISSING and pharma_missing:
        return (
            f"No CRISPR or pharmacologic data available for {row.axis_label}. "
            "Mechanistic rationale exists but no quantitative evidence. "
            "Cannot classify based on current data."
        )

    # CRISPR neg, pharma not examined (guardrail)
    if crispr_neg and pharma_missing:
        return (
            f"CRISPR dependency is negative for {row.axis_label}. However, pharmacologic "
            "screens have NOT been examined. Per doctrine: CRISPR alone cannot rule out "
            "a compound-level vulnerability. Status remains 'Mechanistic candidate only' "
            "pending pharmacologic interrogation."
        )

    # Validated / strong
    if tier == "Validated SL therapeutic lever":
        return (
            f"{row.axis_label}: Validated across multiple modalities including clinical "
            f"and in vivo data. Highest-confidence SL axis. "
            f"Positive modalities: {', '.join(row.positive_modalities())}."
        )

    if tier == "Strong candidate dependency axis":
        return (
            f"{row.axis_label}: CRISPR and/or compound-screen evidence plus mechanism "
            f"support. Positive modalities: {', '.join(row.positive_modalities())}. "
            "Warrants experimental prioritization."
        )

    # Generic
    n_pos = len(row.positive_modalities())
    return (
        f"{row.axis_label}: {n_pos} positive modalities detected. "
        f"Overall evidence: {row.overall_evidence_level}. "
        f"Tier: {row.recommendation_tier}."
    )


def build_agreement_report(matrix: EvidenceMatrix) -> List[ModalityAgreementReport]:
    """
    Build the CRISPR vs pharmacologic agreement report from the fused matrix.
    Surfaces all disagreements and flags those that change the interpretation.
    """
    reports = []
    for row in matrix.rows:
        agreement = row.crispr_pharmacologic_agreement or "insufficient_data"
        cells = row.cells()

        crispr_summ = (
            cells["crispr"].summary
            or f"Δdep={cells['crispr'].delta_dep}, p={cells['crispr'].p_value}"
            if cells["crispr"].delta_dep is not None
            else cells["crispr"].status.value
        )

        prism_summ = cells["prism"].summary or cells["prism"].status.value
        gdsc_summ  = cells["gdsc"].summary  or cells["gdsc"].status.value
        pharma_summ = f"PRISM: {prism_summ} | GDSC: {gdsc_summ}"

        # A disagreement changes interpretation when CRISPR is negative but pharma is positive
        interpretation_change = agreement == "disagree_crispr_neg_pharma_pos"

        if interpretation_change:
            conclusion = (
                "CRISPR and pharmacologic evidence disagree. A compound-level sensitivity "
                "signal is present despite no genetic dependency. This cannot be dismissed as "
                "'no effect' — classified as 'compound-level vulnerability without genetic "
                "dependency.' Functional validation recommended."
            )
        elif agreement == "agree_positive":
            conclusion = (
                "CRISPR and pharmacologic data are concordant (both positive). "
                "Strengthens SL hypothesis."
            )
        elif agreement == "agree_negative":
            conclusion = (
                "Both CRISPR and pharmacologic screens are negative. Low priority axis "
                "unless strong in vitro/in vivo literature receipts exist."
            )
        elif agreement == "disagree_crispr_pos_pharma_neg":
            conclusion = (
                "CRISPR shows dependency but compound screens are flat. May reflect "
                "genetic essentiality without druggable pharmacology, or assay mismatch. "
                "Prioritize in vitro drug testing."
            )
        else:
            conclusion = (
                "Insufficient data to compare modalities. Pharmacologic interrogation needed "
                "before drawing any conclusions about this axis."
            )

        reports.append(ModalityAgreementReport(
            axis_label=row.axis_label,
            agreement_category=agreement,
            crispr_summary=crispr_summ,
            pharmacologic_summary=pharma_summ,
            interpretation_change=interpretation_change,
            conclusion=conclusion,
        ))

    return reports


def build_updated_narrative(
    gene: str,
    cancer_type: Optional[str],
    matrix: EvidenceMatrix,
    guardrails: List[str],
) -> str:
    """
    Build a tumor-board-grade narrative that explicitly references all modalities.
    Never says "CRISPR is negative therefore no vulnerability."
    """
    lines = []
    lines.append(f"=== MULTI-MODAL EVIDENCE NARRATIVE: {gene} ===")
    if cancer_type:
        lines.append(f"Context: {cancer_type}")
    lines.append("")

    # Gold standard first
    cal_row = matrix.get_row(CandidateAxis.CYTIDINE_ANALOGS)
    if cal_row:
        lines.append("── CALIBRATION GOLD STANDARD: Cytidine Analogs ──")
        lines.append(
            f"Evidence level: {cal_row.overall_evidence_level}  |  "
            f"Tier: {cal_row.recommendation_tier}"
        )
        lines.append(
            f"Positive modalities: {', '.join(cal_row.positive_modalities()) or 'none yet evaluated'}"
        )
        if cal_row.interpretation:
            lines.append(cal_row.interpretation)
        lines.append("")

    # All other axes
    lines.append("── CANDIDATE AXES ──")
    for row in matrix.rows:
        if row.axis == CandidateAxis.CYTIDINE_ANALOGS:
            continue
        lines.append(f"\n{row.axis_label}")
        lines.append(f"  Evidence: {row.overall_evidence_level}  |  Tier: {row.recommendation_tier}")
        missing = row.missing_modalities()
        if missing:
            lines.append(f"  Missing modalities (not yet interrogated): {', '.join(missing)}")
        lines.append(f"  CRISPR ↔ Pharmacologic agreement: {row.crispr_pharmacologic_agreement}")
        if row.interpretation:
            lines.append(f"  {row.interpretation}")

    if guardrails:
        lines.append("\n── GUARDRAILS APPLIED ──")
        for g in guardrails:
            lines.append(f"  ⚠ {g}")

    lines.append("\n── DOCTRINE REMINDER ──")
    lines.append(
        "CRISPR dependency is one axis of evidence. A negative CRISPR result does NOT "
        "mean no vulnerability exists, especially when pharmacologic or clinical data "
        "have not been fully interrogated."
    )

    return "\n".join(lines)
