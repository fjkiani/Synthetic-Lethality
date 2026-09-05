from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]
MD=ROOT/"00-mbd4-manuscript/mbd4_parp_response/rxiv/manuscript.md"
FILES=[MD,ROOT/"sl_agent/commercial_rules/sl_evidence_record_v1.json",ROOT/"sl_agent/commercial_rules/mbd4_discovery_evidence.json"]
text="\n".join(p.read_text() for p in FILES)
for forbidden in ["Rodriguez","Giezland","60%","median_pfs_months","os_hazard_ratio","ici_treated_patients"]: assert forbidden not in text,forbidden
for required in ["Rodrigues","Saint-Ghislain","IMMUNOTHERAPY_CHECKPOINT_BLOCKADE","HIGH_CONFIDENCE_TRIAL_ENROLLMENT","five-patient"]: assert required in text,required
payload=json.loads(FILES[1].read_text())
assert payload["denominators"]["tp53_mutant_pharmacologic_WT_reference"]==619
assert payload["therapeutic_axes"]["axis_3"]["primary_action"]=="IMMUNOTHERAPY_CHECKPOINT_BLOCKADE"
print("[VERIFICATION STATUS: SOURCE-GROUNDED TRIAXIS PASSED]")
