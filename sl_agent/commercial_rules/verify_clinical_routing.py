from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sl_agent.commercial_rules.ayesha_therapy_fit import (  # noqa: E402
    ALLOW_PARPI_ROUTING_BYPASS,
    ALLOW_PARPI_TRIAL_EVALUATION,
    CLASS_CONCORDANT_WEE1I_SECONDARY,
    COLORECTAL_TRIAL_TARGET,
    CYTIDINE_ANALOG_SYNTHETIC_LETHALITY,
    HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE,
    HIGH_CONFIDENCE_TRIAL_CANDIDATE,
    HIGH_CONFIDENCE_TRIAL_ENROLLMENT,
    PARP1_GATE_STATUS,
    PRIORITIZE_ATRI_ENROLLMENT,
    SYNERGISTIC_COMBINATION_CANDIDATE,
    ayesha_therapy_fit,
)

RULES_DIR = Path(__file__).resolve().parent
EVIDENCE_PATH = RULES_DIR / "mbd4_discovery_evidence.json"
MANUSCRIPT_PATH = ROOT / "00-mbd4-manuscript/mbd4_parp_response/rxiv/manuscript.md"


def route(patient: Dict[str, Any]) -> Dict[str, Any]:
    return ayesha_therapy_fit(patient)


def exact(result: Dict[str, Any], action: str) -> Dict[str, Any]:
    matches = [item for item in result["routes"] if item["action"] == action]
    assert len(matches) == 1, f"Expected exactly one {action}; got {len(matches)} in {result}"
    return matches[0]


def main() -> None:
    scenarios = {
        "pan_cancer": route({"patient_id": "pan", "mbd4_status": "LOF", "lineage": "Uterus", "PARP1_expression": 7.0}),
        "mss": route({"patient_id": "mss", "mbd4_status": "LOF", "msi_status": "MSS", "PARP1_expression": 7.0}),
        "tp53": route({"patient_id": "tp53", "mbd4_status": "LOF", "tp53_status": "mutant", "lineage": "Ovary", "PARP1_expression": 7.0}),
        "bowel": route({"patient_id": "bowel", "mbd4_status": "LOF", "lineage": "Bowel", "PARP1_expression": 7.0}),
        "heterozygous": route({"patient_id": "het", "mbd4_status": "heterozygous_LikelyLoF", "lineage": "Lung", "PARP1_expression": 7.0}),
        "q75": route({"patient_id": "q75", "mbd4_status": "LOF", "PARP1_expression": 7.41}),
        "brca": route({"patient_id": "brca", "mbd4_status": "LOF", "PARP1_expression": 7.0, "brca1_status": "pathogenic"}),
    }

    exact(scenarios["pan_cancer"], HIGH_CONFIDENCE_TRIAL_CANDIDATE)
    exact(scenarios["mss"], PRIORITIZE_ATRI_ENROLLMENT)
    exact(scenarios["tp53"], HIGH_CONFIDENCE_TRIAL_ENROLLMENT)
    colorectal = exact(scenarios["bowel"], COLORECTAL_TRIAL_TARGET)
    assert (colorectal["canonical_n_biomarker"], colorectal["canonical_n_comparator"]) == (5, 41)
    heterozygous = exact(scenarios["heterozygous"], HIGH_CONFIDENCE_TRIAL_CANDIDATE)
    assert heterozygous["heterozygous_operational_gate"] is True
    assert heterozygous["allele_specific_mechanism_proven"] is False

    parp_block = exact(scenarios["pan_cancer"], HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE)
    assert parp_block["threshold_status"] == PARP1_GATE_STATUS
    assert parp_block["prospective_companion_diagnostic_validation_required"] is True
    exact(scenarios["q75"], ALLOW_PARPI_TRIAL_EVALUATION)
    exact(scenarios["brca"], ALLOW_PARPI_ROUTING_BYPASS)
    exact(scenarios["pan_cancer"], CYTIDINE_ANALOG_SYNTHETIC_LETHALITY)
    exact(scenarios["pan_cancer"], CLASS_CONCORDANT_WEE1I_SECONDARY)
    exact(scenarios["pan_cancer"], SYNERGISTIC_COMBINATION_CANDIDATE)

    evidence = json.loads(EVIDENCE_PATH.read_text())
    serialized = json.dumps(evidence)
    assert "ROUTE_TO_TRIAL_OR_COMPASSIONATE_USE" not in serialized
    actions = {record.get("action") for record in evidence["records"]}
    required_actions = {
        CYTIDINE_ANALOG_SYNTHETIC_LETHALITY,
        HIGH_CONFIDENCE_TRIAL_CANDIDATE,
        PRIORITIZE_ATRI_ENROLLMENT,
        HIGH_CONFIDENCE_TRIAL_ENROLLMENT,
        COLORECTAL_TRIAL_TARGET,
        CLASS_CONCORDANT_WEE1I_SECONDARY,
        SYNERGISTIC_COMBINATION_CANDIDATE,
    }
    assert required_actions <= actions

    manuscript = MANUSCRIPT_PATH.read_text()
    assert "ROUTE_TO_TRIAL_OR_COMPASSIONATE_USE" not in manuscript
    for action in required_actions | {HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE, ALLOW_PARPI_TRIAL_EVALUATION, ALLOW_PARPI_ROUTING_BYPASS}:
        assert action in manuscript, f"Missing manuscript action: {action}"
    for marker in ["n=21", "n=19", "n=14", "n=10", "n=11", "n=5", "n=14 vs 914", "n=11 vs 619", "n=5 vs 41"]:
        assert marker in manuscript, f"Missing denominator marker: {marker}"

    figure_section = manuscript.split("## Figures", 1)[1].split("## Supplementary Material", 1)[0]
    figure_paths = [item.split(")", 1)[0] for item in figure_section.split("](")[1:]]
    assert len(figure_paths) == 4
    for relative in figure_paths:
        assert (MANUSCRIPT_PATH.parent / relative).is_file(), f"Missing figure asset: {relative}"

    print(json.dumps({name: result["flags"] for name, result in scenarios.items()}, indent=2, sort_keys=True))
    print("[VERIFICATION STATUS: CROSS-ARTIFACT RECONCILIATION PASSED]")


if __name__ == "__main__":
    main()
