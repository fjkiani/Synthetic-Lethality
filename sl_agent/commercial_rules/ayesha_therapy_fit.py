from __future__ import annotations

from typing import Any, Dict, List, Optional

CYTIDINE_ANALOG_SYNTHETIC_LETHALITY = "CYTIDINE_ANALOG_SYNTHETIC_LETHALITY"
HIGH_CONFIDENCE_TRIAL_CANDIDATE = "HIGH_CONFIDENCE_TRIAL_CANDIDATE"
PRIORITIZE_ATRI_ENROLLMENT = "PRIORITIZE_ATRI_ENROLLMENT"
HIGH_CONFIDENCE_TRIAL_ENROLLMENT = "HIGH_CONFIDENCE_TRIAL_ENROLLMENT"
COLORECTAL_TRIAL_TARGET = "COLORECTAL_TRIAL_TARGET"
CLASS_CONCORDANT_WEE1I_SECONDARY = "CLASS_CONCORDANT_WEE1I_SECONDARY"
IMMUNOTHERAPY_CHECKPOINT_BLOCKADE = "IMMUNOTHERAPY_CHECKPOINT_BLOCKADE"
HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE = "HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE"
ALLOW_PARPI_TRIAL_EVALUATION = "ALLOW_PARPI_TRIAL_EVALUATION"
ALLOW_PARPI_ROUTING_BYPASS = "ALLOW_PARPI_ROUTING_BYPASS"
SYNERGISTIC_COMBINATION_CANDIDATE = "SYNERGISTIC_COMBINATION_CANDIDATE"

PARP1_Q75 = 7.41
ROUTING_SCOPE = "PRECLINICAL_AND_CLINICAL_TRIAL_PRIORITIZATION"

_LOF = {"loss_of_function", "lof", "likelylof", "biallelic_loss_of_function", "homozygous_loss_of_function"}
_HET_LOF = {"heterozygous_loss_of_function", "heterozygous_lof", "heterozygous_likelylof", "monoallelic_loss_of_function"}
_BOWEL = {"bowel", "colorectal", "colon", "rectal", "crc"}
_TP53_MUT = {"mutant", "mutation", "loss_of_function", "lof", "pathogenic"}
_MSS = {"mss", "microsatellite_stable", "stable"}
_PATHOGENIC = {"pathogenic", "loss_of_function", "lof"}
_HYPERMUTATOR = {"confirmed", "positive", "hypermutator", "cpg_tpg", "cpg_tpg_hypermutator"}


def _norm(value: Optional[str]) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_").replace(">", "_")


def ayesha_therapy_fit(patient: Dict[str, Any], parp1_expression: Optional[float] = None) -> Dict[str, Any]:
    mbd4 = _norm(patient.get("mbd4_status"))
    lineage = _norm(patient.get("lineage"))
    tp53 = _norm(patient.get("tp53_status"))
    msi = _norm(patient.get("msi_status"))
    brca1 = _norm(patient.get("brca1_status"))
    brca2 = _norm(patient.get("brca2_status"))
    hypermutator = _norm(patient.get("hypermutator_status")) in _HYPERMUTATOR or bool(patient.get("cpg_tpg_signature_confirmed", False))
    parp1 = parp1_expression if parp1_expression is not None else patient.get("PARP1_expression")
    parp1 = float(parp1) if parp1 is not None else None

    confirmed_lof = mbd4 in _LOF
    heterozygous_lof = mbd4 in _HET_LOF
    eligible = confirmed_lof or heterozygous_lof
    routes: List[Dict[str, Any]] = []

    if eligible:
        if tp53 in _TP53_MUT:
            atr_action, context = HIGH_CONFIDENCE_TRIAL_ENROLLMENT, "MBD4_LOF_TP53_MUTANT"
        elif msi in _MSS:
            atr_action, context = PRIORITIZE_ATRI_ENROLLMENT, "MBD4_LOF_MSS"
        else:
            atr_action, context = HIGH_CONFIDENCE_TRIAL_CANDIDATE, "MBD4_LOF_PAN_CANCER"
        routes.append({"flag": atr_action, "action": atr_action, "axis": 2, "therapy_class": "ATR inhibitor", "preferred_agent": "ceralasertib", "evidence_context": context, "routing_scope": ROUTING_SCOPE})

        routes.append({"flag": CYTIDINE_ANALOG_SYNTHETIC_LETHALITY, "action": CYTIDINE_ANALOG_SYNTHETIC_LETHALITY, "axis": 1, "therapy_class": "cytidine analog", "example_agents": ["gemcitabine", "cytarabine"], "routing_scope": ROUTING_SCOPE})

        routes.append({"flag": CLASS_CONCORDANT_WEE1I_SECONDARY, "action": CLASS_CONCORDANT_WEE1I_SECONDARY, "axis": 2, "therapy_class": "WEE1 inhibitor", "preferred_agent": "adavosertib", "routing_scope": ROUTING_SCOPE})

        if lineage in _BOWEL:
            routes.append({"flag": COLORECTAL_TRIAL_TARGET, "action": COLORECTAL_TRIAL_TARGET, "axis": 2, "therapy_class": "ATR inhibitor", "preferred_agent": "ceralasertib", "canonical_n_biomarker": 5, "canonical_n_comparator": 41, "routing_scope": ROUTING_SCOPE})

        if hypermutator:
            routes.append({"flag": IMMUNOTHERAPY_CHECKPOINT_BLOCKADE, "action": IMMUNOTHERAPY_CHECKPOINT_BLOCKADE, "axis": 3, "therapy_class": "immune checkpoint blockade", "evidence_context": "MBD4_LOF_CPG_TPG_HYPERMUTATOR", "human_evidence": "retrospective_metastatic_uveal_melanoma_cohort", "prospective_tissue_agnostic_validation_required": True, "routing_scope": ROUTING_SCOPE})

        independent_override = brca1 in _PATHOGENIC or brca2 in _PATHOGENIC or bool(patient.get("validated_independent_parpi_indication", False))
        if independent_override:
            parpi_action = ALLOW_PARPI_ROUTING_BYPASS
        elif parp1 is None or parp1 < PARP1_Q75:
            parpi_action = HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE
        else:
            parpi_action = ALLOW_PARPI_TRIAL_EVALUATION
        routes.append({"flag": parpi_action, "action": parpi_action, "therapy_class": "PARP inhibitor", "parp1_expression": parp1, "parp1_q75_threshold": PARP1_Q75, "threshold_status": "DATASET_DERIVED_EXPLORATORY_THRESHOLD", "prospective_companion_diagnostic_validation_required": True, "routing_scope": ROUTING_SCOPE})

        routes.append({"flag": SYNERGISTIC_COMBINATION_CANDIDATE, "action": SYNERGISTIC_COMBINATION_CANDIDATE, "axes": [1, 2], "therapy_class": "cytidine analog + ATR/WEE1 inhibitor", "combination_synergy_measured": False, "routing_scope": ROUTING_SCOPE})

    return {"patient_id": patient.get("patient_id"), "mbd4_framework_eligible": eligible, "triaxis_model": True, "routes": routes, "flags": [route["flag"] for route in routes]}
