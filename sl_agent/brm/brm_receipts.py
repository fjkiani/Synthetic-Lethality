"""
Frozen PMID-anchored literature receipts for BrM panel genes.

All claims are anchored to specific PMIDs/DOIs. No unsourced claims.
Receipt keys: gene (uppercase). Returns None for unknown genes.

v2: receipts are Layer 3 context (weight 0.25 in exploit_router),
not overrides. DepMap signal takes precedence when available.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class BrMReceipt:
    gene: str
    primary_exploit_mode: str
    secondary_exploit_modes: List[str]
    cascade_step: str
    cascade_step_label: str
    target_node: str
    target_node_role: str
    drug_class: str
    best_drug_example: Optional[str]
    cns_bbb_feasibility: str
    cns_bbb_rationale: str
    context_specificity: str
    evidence_tier: str
    confidence_score: float
    supporting_pmids: List[str]
    primary_evidence_modality: str
    secondary_evidence_modalities: List[str]
    brm_evidence_class: str
    visibility_axis: str
    primary_vulnerability_hypothesis: str
    data_gaps: List[str]
    # v2 annotation flags
    emt_cluster: bool = False
    broad_state_marker: bool = False


_RECEIPTS: Dict[str, BrMReceipt] = {}


def _register(r: BrMReceipt) -> None:
    _RECEIPTS[r.gene.upper()] = r


# ---------------------------------------------------------------------------
# BACE1 — Gold Standard Reference Row
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="BACE1",
    primary_exploit_mode="direct",
    secondary_exploit_modes=["delivery"],
    cascade_step="colonization",
    cascade_step_label="Brain parenchyma colonization",
    target_node="BACE1",
    target_node_role="self",
    drug_class="BACE1 inhibitor",
    best_drug_example="verubecestat",
    cns_bbb_feasibility="established",
    cns_bbb_rationale=(
        "Verubecestat, lanabecestat, and atabecestat were all designed for CNS "
        "penetration in Alzheimer's disease programs; BBB crossing is established."
    ),
    context_specificity="brm_specific",
    evidence_tier="Strong",
    confidence_score=0.82,
    supporting_pmids=["10.1126/scitranslmed.adu2459"],
    primary_evidence_modality="in_vivo_CRISPRa_screen_NSCLC_PDX",
    secondary_evidence_modalities=[
        "genetic_suppression_BrM_in_vivo",
        "pharmacologic_BACE1i_BrM_suppression",
        "patient_NSCLC_BrM_expression",
    ],
    brm_evidence_class="in_vivo_crispr_screen",
    visibility_axis="absent_from_standard_DDR_and_oncology_panels",
    primary_vulnerability_hypothesis=(
        "BACE1 cleaves EGFR to drive brain parenchyma colonization; "
        "CNS-penetrant BACE1 inhibitors (verubecestat class) block this step."
    ),
    data_gaps=[
        "Phase I/II trial of BACE1i in NSCLC BrM patients",
        "Isogenic BACE1-KO vs WT BrM colonization assay with rescue",
    ],
))

# ---------------------------------------------------------------------------
# MMP9 — Delivery Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="MMP9",
    primary_exploit_mode="delivery",
    secondary_exploit_modes=["pathway"],
    cascade_step="bbb_transit",
    cascade_step_label="BBB transit / extravasation",
    target_node="MMP9",
    target_node_role="delivery_interface",
    drug_class="MMP9-selective inhibitor",
    best_drug_example=None,
    cns_bbb_feasibility="partial",
    cns_bbb_rationale=(
        "MMP9 acts at the BBB endothelium; inhibitors must reach the vascular "
        "compartment but need not cross the BBB themselves."
    ),
    context_specificity="brm_specific",
    evidence_tier="Mechanistic",
    confidence_score=0.58,
    supporting_pmids=[
        "10.1158/0008-5472.can-22-3964",
        "10.1158/0008-5472.can-23-0151",
    ],
    primary_evidence_modality="in_vivo_capillary_remodeling_NSCLC",
    secondary_evidence_modalities=["pharmacologic_pan_MMP_inhibitor_mixed"],
    brm_evidence_class="in_vivo_crispr_screen",
    visibility_axis="absent_from_standard_DDR_panels",
    primary_vulnerability_hypothesis=(
        "Cancer cell-derived MMP9 remodels capillary endothelium for brain "
        "colonization; MMP9-selective inhibition blocks BBB transit."
    ),
    data_gaps=[
        "MMP9-selective inhibitor (not pan-MMP) in NSCLC BrM model",
        "Patient MMP9 expression stratification in BrM vs primary NSCLC",
    ],
))

# ---------------------------------------------------------------------------
# CLDN5 — Delivery Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="CLDN5",
    primary_exploit_mode="delivery",
    secondary_exploit_modes=[],
    cascade_step="bbb_transit",
    cascade_step_label="BBB transit / tight junction modulation",
    target_node="CLDN5",
    target_node_role="delivery_interface",
    drug_class="claudin-targeting antibody",
    best_drug_example=None,
    cns_bbb_feasibility="barrier",
    cns_bbb_rationale=(
        "CLDN5 maintains BBB integrity; targeting it opens the barrier but "
        "creates CNS toxicity risk requiring tumor-cell-specific strategy."
    ),
    context_specificity="cns_relevant",
    evidence_tier="Mechanistic",
    confidence_score=0.38,
    supporting_pmids=["10.1038/s41586-020-2698-6"],
    primary_evidence_modality="genetic_suppression_BBB_permeability",
    secondary_evidence_modalities=[],
    brm_evidence_class="bbb_biology",
    visibility_axis="absent_from_standard_oncology_panels",
    primary_vulnerability_hypothesis=(
        "CLDN5 downregulation in BrM-associated endothelium increases BBB "
        "permeability; tumor-cell-specific CLDN5 modulation may block transit."
    ),
    data_gaps=[
        "Tumor-cell-specific CLDN5 modulation strategy",
        "In vivo BrM model with CLDN5 perturbation",
    ],
))

# ---------------------------------------------------------------------------
# CCL2 — Pathway Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="CCL2",
    primary_exploit_mode="pathway",
    secondary_exploit_modes=["delivery"],
    cascade_step="immune_evasion",
    cascade_step_label="Immune evasion / monocyte recruitment",
    target_node="CCR2",
    target_node_role="pathway_bottleneck",
    drug_class="CCR2 antagonist",
    best_drug_example="BMS-813160",
    cns_bbb_feasibility="partial",
    cns_bbb_rationale=(
        "CCL2/CCR2 axis recruits monocytes across the BBB; "
        "CCR2 antagonists act at the vascular interface."
    ),
    context_specificity="brm_specific",
    evidence_tier="Mechanistic",
    confidence_score=0.44,
    supporting_pmids=["10.1038/s41591-019-0392-0"],
    primary_evidence_modality="patient_expression_BrM_TME",
    secondary_evidence_modalities=["pharmacologic_CCR2_antagonist_preclinical"],
    brm_evidence_class="patient_expression",
    visibility_axis="absent_from_standard_DDR_panels",
    primary_vulnerability_hypothesis=(
        "CCL2 drives monocyte recruitment into BrM TME via CCR2; "
        "CCR2 antagonism (BMS-813160) blocks immunosuppressive niche formation."
    ),
    data_gaps=[
        "CCR2 antagonist in NSCLC BrM model with immune profiling",
        "Patient CCL2 expression stratification vs BrM incidence",
    ],
))

# ---------------------------------------------------------------------------
# ICAM1 — Pathway Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="ICAM1",
    primary_exploit_mode="pathway",
    secondary_exploit_modes=["delivery"],
    cascade_step="bbb_transit",
    cascade_step_label="BBB transit / transendothelial migration",
    target_node="ICAM1",
    target_node_role="microenvironment",
    drug_class="anti-ICAM1 antibody",
    best_drug_example=None,
    cns_bbb_feasibility="partial",
    cns_bbb_rationale=(
        "ICAM1 is expressed on BBB endothelium; antibody-accessible at "
        "vascular surface without requiring CNS penetration."
    ),
    context_specificity="brm_specific",
    evidence_tier="Mechanistic",
    confidence_score=0.41,
    supporting_pmids=["10.1158/0008-5472.can-20-1236"],
    primary_evidence_modality="patient_expression_BrM",
    secondary_evidence_modalities=["genetic_suppression_transendothelial_migration"],
    brm_evidence_class="patient_expression",
    visibility_axis="absent_from_standard_DDR_panels",
    primary_vulnerability_hypothesis=(
        "ICAM1 mediates cancer cell transendothelial migration across the BBB; "
        "anti-ICAM1 antibody blocks this step at the vascular surface."
    ),
    data_gaps=[
        "Anti-ICAM1 antibody in NSCLC BrM in vivo model",
        "ICAM1 expression stratification vs BrM incidence in NSCLC patients",
    ],
    broad_state_marker=True,  # POL-005: zero FDR<0.10 DepMap partners
))

# ---------------------------------------------------------------------------
# MMP2 — Delivery Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="MMP2",
    primary_exploit_mode="delivery",
    secondary_exploit_modes=["pathway"],
    cascade_step="bbb_transit",
    cascade_step_label="BBB transit / basement membrane degradation",
    target_node="MMP2",
    target_node_role="delivery_interface",
    drug_class="MMP2-selective inhibitor",
    best_drug_example=None,
    cns_bbb_feasibility="partial",
    cns_bbb_rationale=(
        "MMP2 degrades BBB basement membrane; inhibitors act at the "
        "vascular compartment and need not cross the BBB."
    ),
    context_specificity="brm_specific",
    evidence_tier="Mechanistic",
    confidence_score=0.40,
    supporting_pmids=["10.1158/0008-5472.can-22-3964"],
    primary_evidence_modality="in_vivo_basement_membrane_remodeling",
    secondary_evidence_modalities=["pharmacologic_pan_MMP_inhibitor_mixed"],
    brm_evidence_class="in_vivo_crispr_screen",
    visibility_axis="absent_from_standard_DDR_panels",
    primary_vulnerability_hypothesis=(
        "MMP2 degrades BBB basement membrane for cancer cell extravasation; "
        "MMP2-selective inhibition blocks this step. Dual MMP2/MMP9 targeting may be required."
    ),
    data_gaps=[
        "MMP2-selective inhibitor in NSCLC BrM model",
        "MMP2 vs MMP9 relative contribution to BBB degradation in NSCLC",
    ],
))

# ---------------------------------------------------------------------------
# TWIST1 — Pathway Exploit
# ---------------------------------------------------------------------------
_register(BrMReceipt(
    gene="TWIST1",
    primary_exploit_mode="pathway",
    secondary_exploit_modes=[],
    cascade_step="multi_step",
    cascade_step_label="EMT / intravasation + colonization",
    target_node="BRD4",
    target_node_role="pathway_bottleneck",
    drug_class="BET bromodomain inhibitor",
    best_drug_example="JQ1 / OTX015",
    cns_bbb_feasibility="partial",
    cns_bbb_rationale=(
        "BETi CNS penetration is variable by compound; "
        "JQ1 has demonstrated CNS exposure in preclinical models."
    ),
    context_specificity="pan_cancer",
    evidence_tier="Mechanistic",
    confidence_score=0.33,
    supporting_pmids=["10.1038/s41467-019-09270-4"],
    primary_evidence_modality="genetic_suppression_EMT_BrM",
    secondary_evidence_modalities=["pharmacologic_BETi_EMT_reversal"],
    brm_evidence_class="genetic_suppression",
    visibility_axis="absent_from_standard_DDR_panels",
    primary_vulnerability_hypothesis=(
        "TWIST1 drives EMT enabling intravasation and colonization; "
        "BRD4 inhibition (JQ1/OTX015) suppresses TWIST1 transcription."
    ),
    data_gaps=[
        "TWIST1-driven BrM model with BETi treatment",
        "HDAC inhibitor as alternative TWIST1 pathway node",
    ],
))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_brm_receipt(gene: str) -> Optional[BrMReceipt]:
    """Return frozen receipt for gene (case-insensitive). None if not found."""
    return _RECEIPTS.get(gene.upper())


def list_receipt_genes() -> List[str]:
    """Return all genes with frozen receipts."""
    return list(_RECEIPTS.keys())


def all_receipts() -> Dict[str, BrMReceipt]:
    """Return all receipts keyed by uppercase gene symbol."""
    return dict(_RECEIPTS)
