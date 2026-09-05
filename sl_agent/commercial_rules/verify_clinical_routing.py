from __future__ import annotations
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from sl_agent.commercial_rules.ayesha_therapy_fit import *
RULES=Path(__file__).resolve().parent
MD=ROOT/"00-mbd4-manuscript/mbd4_parp_response/rxiv/manuscript.md"
EVIDENCE=RULES/"mbd4_discovery_evidence.json"

def one(result,action):
    matches=[r for r in result["routes"] if r["action"]==action]
    assert len(matches)==1,(action,result)
    return matches[0]

def main():
    base={"mbd4_status":"LOF","lineage":"Uterus","PARP1_expression":7.0}
    one(ayesha_therapy_fit(base),HIGH_CONFIDENCE_TRIAL_CANDIDATE)
    one(ayesha_therapy_fit({**base,"msi_status":"MSS"}),PRIORITIZE_ATRI_ENROLLMENT)
    one(ayesha_therapy_fit({**base,"tp53_status":"mutant"}),HIGH_CONFIDENCE_TRIAL_ENROLLMENT)
    one(ayesha_therapy_fit({**base,"lineage":"Bowel"}),COLORECTAL_TRIAL_TARGET)
    one(ayesha_therapy_fit(base),CYTIDINE_ANALOG_SYNTHETIC_LETHALITY)
    one(ayesha_therapy_fit(base),CLASS_CONCORDANT_WEE1I_SECONDARY)
    ici=one(ayesha_therapy_fit({**base,"cpg_tpg_signature_confirmed":True}),IMMUNOTHERAPY_CHECKPOINT_BLOCKADE)
    assert ici["prospective_tissue_agnostic_validation_required"] is True
    assert IMMUNOTHERAPY_CHECKPOINT_BLOCKADE not in ayesha_therapy_fit(base)["flags"]
    one(ayesha_therapy_fit(base),HARD_BLOCK_LACKS_TRAPPING_SUBSTRATE)
    one(ayesha_therapy_fit({**base,"PARP1_expression":7.41}),ALLOW_PARPI_TRIAL_EVALUATION)
    one(ayesha_therapy_fit({**base,"brca1_status":"pathogenic"}),ALLOW_PARPI_ROUTING_BYPASS)
    one(ayesha_therapy_fit(base),SYNERGISTIC_COMBINATION_CANDIDATE)
    evidence=json.loads(EVIDENCE.read_text())
    assert evidence["therapeutic_model"]=="MBD4_LOF_TRIAXIS"
    immuno=next(r for r in evidence["records"] if r.get("action")==IMMUNOTHERAPY_CHECKPOINT_BLOCKADE)
    cohort=immuno["human_evidence"]["retrospective_cohort"]
    assert (cohort["mbd4_sequenced"],cohort["mbd4_mutated"])==(131,5)
    assert (cohort["objective_response_mbd4_mutated"],cohort["median_pfs_months_mbd4_mutated"])==(0.60,22.3)
    md=MD.read_text()
    assert "MBD4-LOF Triaxis Therapeutic Vulnerability" in md
    assert "Axis 3: CpG>TpG hypermutator phenotype and clinical immune checkpoint blockade" in md
    assert "IMMUNOTHERAPY_CHECKPOINT_BLOCKADE" in md
    legacy_title = "Dual" + " Therapeutic Vulnerability"
    assert legacy_title not in md
    figpart=md.split("## Figures",1)[1].split("## Supplementary Material",1)[0]
    paths=[x.split(")",1)[0] for x in figpart.split("](")[1:]]
    assert len(paths)==4
    for p in paths:
        asset=MD.parent/p
        assert asset.is_file() and asset.stat().st_size>10000,p
    print(json.dumps({"triaxis":True,"figure_assets":paths,"immunotherapy_action":IMMUNOTHERAPY_CHECKPOINT_BLOCKADE},indent=2))
    print("[VERIFICATION STATUS: TRIAXIS CROSS-ARTIFACT RECONCILIATION PASSED]")
if __name__=="__main__": main()
