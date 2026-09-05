from __future__ import annotations

from typing import Any, Dict, List, Optional

CYTIDINE_ANALOG_SYNTHETIC_LETHALITY = "CYTIDINE_ANALOG_SYNTHETIC_LETHALITY"
HIGH_CONFIDENCE_TRIAL_CANDIDATE = "HIGH_CONFIDENCE_TRIAL_CANDIDATE"
PRIORITIZE_ATRI_ENROLLMENT = "PRIORITIZE_ATRI_ENROLLMENT"
HIGH_CONFIDENCE_TRIAL_ENROLLMENT = "HIGH_CONFIDENCE_TRIAL_ENROLLMENT"
COLORECTAL_TRIAL_TARGET = "COLORECTAL_TRIAL_TARGET"
CLASS_CONCORDANT_WEE1I_SECONDARY = "CLASS_CONCORDANT_WEE1I_SECONDARY"
HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE = "HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE"
ALLOW_PARPI_TRIAL_EVALUATION = "ALLOW_PARPI_TRIAL_EVALUATION"
ALLOW_PARPI_ROUTING_BYPASS = "ALLOW_PARPI_ROUTING_BYPASS"
SYNERGISTIC_COMBINATION_CANDIDATE = "SYNERGISTIC_COMBINATION_CANDIDATE"

PARP1_Q75 = 7.41
PARP1_GATE_STATUS = "DATASET_DERIVED_EXPLORATORY_THRESHOLD"
ROUTING_SCOPE = "PRECLINICAL_TRIAL_PRIORITIZATION"

_MBD4_LOF_STATES = {
    "loss_of_function",
    "lof",
    "likelylof",
    "biallelic_loss_of_function",
    "homozygous_loss_of_function",
}
_MBD4_HET_LOF_STATES = {
    "heterozygous_loss_of_function",
    "heterozygous_lof",
    "heterozygous_likelylof",
    "monoallelic_loss_of_function",
}
_BOWEL_LINEAGES = {"bowel", "colorectal", "colon", "rectal", "crc"}
_TP53_MUTANT_STATES = {"mutant", "mutation", "loss_of_function", "lof", "pathogenic"}
_MSS_STATES = {"mss", "microsatellite_stable", "stable"}
_PATHOGENIC_STATES = {"pathogenic", "loss_of_function", "lof"}


def _norm(value: Optional[str]) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def ayesha_therapy_fit(
    patient: Dict[str, Any],
    parp1_expression: Optional[float] = None,
) -> Dict[str, Any]:
    mbd4_state = _norm(patient.get("mbd4_status"))
    lineage = _norm(patient.get("lineage"))
    tp53_state = _norm(patient.get("tp53_status"))
    msi_state = _norm(patient.get("msi_status"))
    brca1_state = _norm(patient.get("brca1_status"))
    brca2_state = _norm(patient.get("brca2_status"))
    independent_parpi_indication = bool(patient.get("validated_independent_parpi_indication", False))
    parp1_value = parp1_expression if parp1_expression is not None else patient.get("PARP1_expression")
    parp1_value = float(parp1_value) if parp1_value is not None else None

    confirmed_lof = mbd4_state in _MBD4_LOF_STATES
    heterozygous_lof = mbd4_state in _MBD4_HET_LOF_STATES
    mbd4_framework_eligible = confirmed_lof or heterozygous_lof
    bowel_context = lineage in _BOWEL_LINEAGES
    tp53_mutant = tp53_state in _TP53_MUTANT_STATES
    mss_context = msi_state in _MSS_STATES
    brca_override = brca1_state in _PATHOGENIC_STATES or brca2_state in _PATHOGENIC_STATES
    validated_parpi_override = brca_override or independent_parpi_indication

    routes: List[Dict[str, Any]] = []

    if mbd4_framework_eligible:
        if tp53_mutant:
            atr_action = HIGH_CONFIDENCE_TRIAL_ENROLLMENT
            evidence_context = "MBD4_LOF_TP53_MUTANT"
        elif mss_context:
            atr_action = PRIORITIZE_ATRI_ENROLLMENT
            evidence_context = "MBD4_LOF_MSS"
        else:
            atr_action = HIGH_CONFIDENCE_TRIAL_CANDIDATE
            evidence_context = "MBD4_LOF_PAN_CANCER"

        routes.append({
            "flag": atr_action,
            "action": atr_action,
            "routing_scope": ROUTING_SCOPE,
            "therapy_class": "ATR inhibitor",
            "preferred_agent": "ceralasertib",
            "evidence_context": evidence_context,
            "heterozygous_operational_gate": heterozygous_lof,
            "allele_specific_mechanism_proven": False,
        })

        if bowel_context:
            routes.append({
                "flag": COLORECTAL_TRIAL_TARGET,
                "action": COLORECTAL_TRIAL_TARGET,
                "routing_scope": ROUTING_SCOPE,
                "therapy_class": "ATR inhibitor",
                "preferred_agent": "ceralasertib",
                "canonical_n_biomarker": 5,
                "canonical_n_comparator": 41,
                "delta_mean_ln_ic50": -0.6922974634146344,
                "p_value_one_sided": 0.1263800798684519,
                "cohens_d": -0.4638353720829927,
            })

        if validated_parpi_override:
            parpi_action = ALLOW_PARPI_ROUTING_BYPASS
            basis = "BRCA1_BRCA2_or_other_validated_independent_indication"
        elif parp1_value is None or parp1_value < PARP1_Q75:
            parpi_action = HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE
            basis = "PARP1_expression_missing_or_below_dataset_Q75"
        else:
            parpi_action = ALLOW_PARPI_TRIAL_EVALUATION
            basis = "PARP1_expression_at_or_above_dataset_Q75"

        routes.append({
            "flag": parpi_action,
            "action": parpi_action,
            "routing_scope": ROUTING_SCOPE,
            "therapy_class": "PARP inhibitor",
            "parp1_expression": parp1_value,
            "parp1_q75_threshold": PARP1_Q75,
            "threshold_scale": "DepMap_TPM_log1p",
            "threshold_status": PARP1_GATE_STATUS,
            "prospective_companion_diagnostic_validation_required": True,
            "basis": basis,
        })

        routes.append({
            "flag": SYNERGISTIC_COMBINATION_CANDIDATE,
            "action": SYNERGISTIC_COMBINATION_CANDIDATE,
            "routing_scope": ROUTING_SCOPE,
            "therapy_class": "cytidine analog + ATR inhibitor",
            "example_regimens": [
                ["gemcitabine", "ceralasertib"],
                ["cytarabine", "ceralasertib"],
            ],
            "mechanistic_basis": "BER substrate accumulation plus ATR checkpoint blockade converges on replication-fork failure",
            "combination_synergy_measured": False,
        })

    if confirmed_lof:
        routes.append({
            "flag": CYTIDINE_ANALOG_SYNTHETIC_LETHALITY,
            "action": CYTIDINE_ANALOG_SYNTHETIC_LETHALITY,
            "routing_scope": ROUTING_SCOPE,
            "therapy_class": "cytidine analog",
            "example_agents": ["gemcitabine", "cytarabine"],
            "evidence_grade": "VALIDATED_PRECLINICAL_AXIS",
        })
        routes.append({
            "flag": CLASS_CONCORDANT_WEE1I_SECONDARY,
            "action": CLASS_CONCORDANT_WEE1I_SECONDARY,
            "routing_scope": ROUTING_SCOPE,
            "therapy_class": "WEE1 inhibitor",
            "preferred_agent": "adavosertib",
            "evidence_grade": "DISCOVERY_GRADE_ACTIONABLE",
        })

    return {
        "patient_id": patient.get("patient_id"),
        "routing_scope": ROUTING_SCOPE,
        "mbd4_framework_eligible": mbd4_framework_eligible,
        "routes": routes,
        "flags": [route["flag"] for route in routes],
    }
