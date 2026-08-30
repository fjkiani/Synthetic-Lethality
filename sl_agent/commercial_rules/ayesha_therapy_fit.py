from __future__ import annotations

from typing import Any, Dict, List, Optional

DISCOVERY_GRADE_ACTIONABLE = "DISCOVERY_GRADE_ACTIONABLE"
HIGH_CONFIDENCE_TRIAL_CANDIDATE = "HIGH_CONFIDENCE_TRIAL_CANDIDATE"
HARD_ROUTING_RULE = "HARD_ROUTING_RULE"

_MBD4_LOF_STATES = {
    "loss_of_function",
    "lof",
    "biallelic_loss_of_function",
    "homozygous_loss_of_function",
}
_MBD4_HET_LOF_STATES = {
    "heterozygous_loss_of_function",
    "heterozygous_lof",
    "monoallelic_loss_of_function",
}
_BOWEL_LINEAGES = {"bowel", "colorectal", "colon", "rectal", "crc"}


def _norm(value: Optional[str]) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def ayesha_therapy_fit(patient: Dict[str, Any]) -> Dict[str, Any]:
    mbd4_state = _norm(patient.get("mbd4_status"))
    lineage = _norm(patient.get("lineage"))
    tp53_state = _norm(patient.get("tp53_status"))
    brca1_state = _norm(patient.get("brca1_status"))
    brca2_state = _norm(patient.get("brca2_status"))
    independent_parpi_indication = bool(patient.get("validated_independent_parpi_indication", False))

    confirmed_lof = mbd4_state in _MBD4_LOF_STATES
    heterozygous_lof = mbd4_state in _MBD4_HET_LOF_STATES
    mbd4_framework_eligible = confirmed_lof or heterozygous_lof
    bowel_context = lineage in _BOWEL_LINEAGES
    tp53_mutant = tp53_state in {"mutant", "mutation", "loss_of_function", "lof", "pathogenic"}
    brca_override = brca1_state in {"pathogenic", "loss_of_function", "lof"} or brca2_state in {"pathogenic", "loss_of_function", "lof"}
    validated_parpi_override = brca_override or independent_parpi_indication

    routes: List[Dict[str, Any]] = []

    cytidine_qualified = confirmed_lof
    atr_qualified = mbd4_framework_eligible

    if atr_qualified:
        routes.append({
            "flag": "POTENTIAL_ATRI_SENSITIVITY",
            "status": DISCOVERY_GRADE_ACTIONABLE if heterozygous_lof else HIGH_CONFIDENCE_TRIAL_CANDIDATE,
            "therapy_class": "ATR inhibitor",
            "preferred_agent": "ceralasertib",
            "evidence_context": "MBD4_HETEROZYGOUS" if heterozygous_lof else "MBD4_LOF",
            "tp53_enriched_context": tp53_mutant,
            "bowel_discovery_context": bowel_context,
            "trial_action": "PRIORITIZE_ATRI_ENROLLMENT",
            "pending_validation": True,
        })

    if mbd4_framework_eligible:
        if validated_parpi_override:
            routes.append({
                "flag": "PRESERVE_INDEPENDENT_PARPI_INDICATION",
                "status": HARD_ROUTING_RULE,
                "therapy_class": "PARP inhibitor",
                "basis": "BRCA1_BRCA2_or_other_validated_independent_indication",
            })
        else:
            routes.append({
                "flag": "HARD_BLOCK_MBD4_ONLY_PARPI_ROUTE",
                "status": HARD_ROUTING_RULE,
                "therapy_class": "PARP inhibitor",
                "basis": "MBD4_LOF_does_not_upregulate_PARP1_and_is_not_a_sufficient_PARPi_biomarker",
            })

    if confirmed_lof:
        routes.append({
            "flag": "CYTIDINE_ANALOG_SYNTHETIC_LETHALITY",
            "status": "VALIDATED_PRECLINICAL_AXIS",
            "therapy_class": "cytidine analog",
            "example_agents": ["gemcitabine", "cytarabine"],
            "pending_validation": True,
        })
        routes.append({
            "flag": "CLASS_CONCORDANT_WEE1I_SECONDARY",
            "status": DISCOVERY_GRADE_ACTIONABLE,
            "therapy_class": "WEE1 inhibitor",
            "preferred_agent": "adavosertib",
            "trial_action": "PRIORITIZE_WEE1I_COHORT_OR_ARM",
            "pending_validation": True,
        })

    if bowel_context and mbd4_framework_eligible:
        routes.append({
            "flag": "COLORECTAL_PRIMARY_ATRI_TRIAL_LINEAGE",
            "status": DISCOVERY_GRADE_ACTIONABLE,
            "therapy_class": "ATR inhibitor",
            "preferred_agent": "ceralasertib",
            "fold_lower_geometric_mean_ic50": 1.9983012875241446,
            "canonical_n_biomarker": 5,
            "canonical_n_comparator": 41,
            "trial_action": "PRIORITIZE_COLORECTAL_ATRI_COHORT",
            "pending_validation": True,
        })

    if mbd4_framework_eligible and (cytidine_qualified or atr_qualified):
        routes.append({
            "flag": "SYNERGISTIC_COMBINATION_CANDIDATE",
            "status": DISCOVERY_GRADE_ACTIONABLE,
            "therapy_class": "cytidine analog + ATR inhibitor",
            "example_regimens": [
                ["gemcitabine", "ceralasertib"],
                ["cytarabine", "ceralasertib"],
            ],
            "mechanistic_basis": "BER substrate accumulation plus ATR checkpoint blockade converges on replication-fork failure",
            "combination_synergy_measured": False,
            "pending_validation": True,
        })

    return {
        "patient_id": patient.get("patient_id"),
        "mbd4_framework_eligible": mbd4_framework_eligible,
        "routes": routes,
        "flags": [route["flag"] for route in routes],
    }
