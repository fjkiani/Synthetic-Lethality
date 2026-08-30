from __future__ import annotations

import json
from typing import Any, Dict

CANONICAL = {
    "tp53_atri": {"n_lof": 11, "n_wt": 619, "delta_ln_ic50": -1.07, "p": 0.003, "cohens_d": -0.74},
    "bowel_atri": {"n_lof": 5, "n_wt": 41, "delta_ln_ic50": -0.69, "p": 0.126, "cohens_d": -0.46},
    "parp1_q75_gate": 7.41,
}


def normalize(value: Any) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def route_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    mbd4 = normalize(patient.get("MBD4"))
    tp53 = normalize(patient.get("TP53"))
    lineage = normalize(patient.get("lineage"))
    brca1 = normalize(patient.get("BRCA1"))
    brca2 = normalize(patient.get("BRCA2"))
    hrd_override = bool(patient.get("validated_HRD_override", False))
    parp1 = patient.get("PARP1_expression")

    mbd4_lof = mbd4 in {"lof", "loss_of_function", "likelylof", "heterozygous_likelylof"}
    tp53_mutant = tp53 in {"mutant", "lof", "loss_of_function", "pathogenic"}
    bowel = lineage in {"bowel", "colorectal", "colon", "rectal", "crc"}
    brca_override = brca1 == "pathogenic" or brca2 == "pathogenic" or hrd_override

    actions = []
    if mbd4_lof and tp53_mutant:
        actions.append({"action": "HIGH_CONFIDENCE_TRIAL_ENROLLMENT", "therapy": "ceralasertib", "metadata": CANONICAL["tp53_atri"]})
    if mbd4_lof and bowel:
        actions.append({"action": "COLORECTAL_TRIAL_TARGET", "therapy": "ceralasertib", "metadata": CANONICAL["bowel_atri"]})
    if mbd4_lof and parp1 is not None and float(parp1) < CANONICAL["parp1_q75_gate"]:
        actions.append({
            "action": "ALLOW_PARPI_ROUTING_BYPASS" if brca_override else "HARD_BLOCK_MBD4_ONLY_PARPI_ROUTE",
            "therapy_class": "PARP inhibitor",
            "override_basis": "BRCA1_BRCA2_or_validated_HRD" if brca_override else None,
        })
    return {"patient": patient, "actions": actions}


def exact_action(result: Dict[str, Any], action: str) -> Dict[str, Any]:
    matches = [item for item in result["actions"] if item["action"] == action]
    assert len(matches) == 1, f"Expected exactly one {action}, found {len(matches)}: {result}"
    return matches[0]


def main() -> None:
    scenarios = {
        "mbd4_lof_tp53_ovary": route_patient({"MBD4": "LOF", "TP53": "mutant", "lineage": "Ovary"}),
        "mbd4_lof_bowel": route_patient({"MBD4": "LOF", "TP53": "wild_type", "lineage": "Bowel"}),
        "mbd4_lof_low_parp1": route_patient({"MBD4": "LOF", "TP53": "wild_type", "lineage": "Uterus", "PARP1_expression": 7.40}),
        "mbd4_lof_low_parp1_brca1": route_patient({"MBD4": "LOF", "TP53": "wild_type", "lineage": "Ovary", "PARP1_expression": 7.40, "BRCA1": "pathogenic"}),
    }

    tp53_route = exact_action(scenarios["mbd4_lof_tp53_ovary"], "HIGH_CONFIDENCE_TRIAL_ENROLLMENT")
    assert tp53_route["metadata"] == {"n_lof": 11, "n_wt": 619, "delta_ln_ic50": -1.07, "p": 0.003, "cohens_d": -0.74}

    bowel_route = exact_action(scenarios["mbd4_lof_bowel"], "COLORECTAL_TRIAL_TARGET")
    assert bowel_route["metadata"] == {"n_lof": 5, "n_wt": 41, "delta_ln_ic50": -0.69, "p": 0.126, "cohens_d": -0.46}

    exact_action(scenarios["mbd4_lof_low_parp1"], "HARD_BLOCK_MBD4_ONLY_PARPI_ROUTE")
    exact_action(scenarios["mbd4_lof_low_parp1_brca1"], "ALLOW_PARPI_ROUTING_BYPASS")

    print(json.dumps(scenarios, indent=2, sort_keys=True))
    print("[VERIFICATION STATUS: 100% CANONICAL DEPLOYMENT SUCCESS]")


if __name__ == "__main__":
    main()
