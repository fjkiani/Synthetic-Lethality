"""
Multi-modal evidence models.

Core doctrine:
  - CRISPR = one axis of evidence, not the arbiter of truth.
  - Compound-response (PRISM, GDSC) and literature receipts can override or
    nuance CRISPR "no-dependency" when they are mechanistically grounded.
  - Cytidine analogs for MBD4 are the calibration gold standard.
  - Recommendation levels are driven by the matrix, not by any single modality.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


# ── Modality enum ─────────────────────────────────────────────────────────────

class Modality(str, Enum):
    CRISPR_DEPENDENCY   = "crispr_dependency"
    EXPRESSION_ASSOC    = "expression_association"   # RNA/protein co-variation
    PRISM_PHARMACOLOGIC = "prism_pharmacologic"       # PRISM repurposing screen
    GDSC_PHARMACOLOGIC  = "gdsc_pharmacologic"        # GDSC IC50 / AUC
    IN_VITRO_FUNCTIONAL = "in_vitro_functional"       # KO/WT, rescue, isogenic
    IN_VIVO_PDX         = "in_vivo_pdx"               # PDX / xenograft
    CLINICAL            = "clinical"                   # patient outcome data
    PATHWAY_MECHANISM   = "pathway_mechanism"          # known pathway biology


# ── Evidence status for a single modality cell ────────────────────────────────

class ModalityStatus(str, Enum):
    POSITIVE      = "positive"       # clear supporting signal
    NEGATIVE      = "negative"       # tested, no signal
    MIXED         = "mixed"          # contradictory / context-dependent
    MISSING       = "missing"        # not interrogated / data unavailable
    CONFOUNDED    = "confounded"     # signal exists but is MSI/MMR or other confound


# ── Single modality evidence cell ─────────────────────────────────────────────

class ModalityEvidence(BaseModel):
    """One cell in the evidence matrix: (candidate_axis, modality)."""

    modality: Modality
    status: ModalityStatus = ModalityStatus.MISSING

    # Quantitative fields — only those relevant to this modality are filled
    delta_dep: Optional[float] = Field(None, description="Chronos delta-dependency (CRISPR)")
    p_value: Optional[float] = None
    fdr: Optional[float] = None
    effect_size: Optional[float] = Field(None, description="Cohen's d or similar")
    n_mut: Optional[int] = None
    n_wt: Optional[int] = None

    # Pharmacologic fields
    delta_auc: Optional[float] = Field(None, description="ΔAUC or ΔZ-score (drug screen)")
    delta_ic50_log2: Optional[float] = Field(None, description="log2 fold-change in IC50")
    drug_screen_dataset: Optional[str] = None   # "PRISM_repurposing", "GDSC1", etc.
    stratifier: Optional[str] = Field(
        None,
        description="Stratification used: 'MBD4_LOF_vs_WT', 'MSI_H_vs_MSS', etc.",
    )

    # Narrative / receipt
    summary: Optional[str] = Field(None, description="Human-readable evidence summary")
    pmids: List[str] = Field(default_factory=list)
    is_confound_flagged: bool = Field(
        False, description="True if signal may be MSI/MMR confounded not MBD4-specific"
    )
    notes: Optional[str] = None


# ── Candidate axis row ────────────────────────────────────────────────────────

class CandidateAxis(str, Enum):
    """Standard candidate therapy axes for MBD4-related biology."""
    CYTIDINE_ANALOGS = "cytidine_analogs"    # Gemcitabine, cytarabine — calibration gold std
    PARP_INHIBITORS  = "parp_inhibitors"     # Olaparib, niraparib, talazoparib, rucaparib
    ATR_WEE1         = "atr_wee1"            # ATR/WEE1 inhibitors (checkpoint replication)
    WRN              = "wrn"                 # WRN helicase (MSI-context SL)
    IMMUNOTHERAPY    = "immunotherapy"       # PD-1/PD-L1 (hypermutator / MSI-H)
    PKMYT1           = "pkmyt1"              # PKMYT1 checkpoint kinase — RS-specific SL
    CUSTOM           = "custom"              # Any gene/drug pair injected at query time


# ── Full row in the evidence matrix ──────────────────────────────────────────

class EvidenceRow(BaseModel):
    """One candidate axis with all modality cells filled or marked MISSING."""

    axis: CandidateAxis
    axis_label: str   # human-readable (e.g. "PARP Inhibitors (olaparib, niraparib)")
    mechanism: str    # pathway biology / known roles

    # One cell per modality
    crispr: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.CRISPR_DEPENDENCY)
    )
    expression: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.EXPRESSION_ASSOC)
    )
    prism: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.PRISM_PHARMACOLOGIC)
    )
    gdsc: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.GDSC_PHARMACOLOGIC)
    )
    in_vitro: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.IN_VITRO_FUNCTIONAL)
    )
    in_vivo: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.IN_VIVO_PDX)
    )
    clinical: ModalityEvidence = Field(
        default_factory=lambda: ModalityEvidence(modality=Modality.CLINICAL)
    )

    # Derived
    overall_evidence_level: Optional[str] = Field(
        None,
        description="High / Moderate / Mechanistic-only / Negative (computed by fuser)",
    )
    recommendation_tier: Optional[str] = Field(
        None,
        description="Validated SL lever / Strong candidate / Mechanistic candidate / Not supported",
    )
    interpretation: Optional[str] = None
    crispr_pharmacologic_agreement: Optional[str] = Field(
        None,
        description="'agree_positive' | 'agree_negative' | 'disagree_crispr_neg_pharma_pos' | "
                    "'disagree_crispr_pos_pharma_neg' | 'insufficient_data'",
    )

    def cells(self) -> Dict[str, ModalityEvidence]:
        return {
            "crispr": self.crispr,
            "expression": self.expression,
            "prism": self.prism,
            "gdsc": self.gdsc,
            "in_vitro": self.in_vitro,
            "in_vivo": self.in_vivo,
            "clinical": self.clinical,
        }

    def positive_modalities(self) -> List[str]:
        return [k for k, v in self.cells().items() if v.status == ModalityStatus.POSITIVE]

    def missing_modalities(self) -> List[str]:
        return [k for k, v in self.cells().items() if v.status == ModalityStatus.MISSING]


# ── Combinatorial vulnerability ───────────────────────────────────────────────

class CombinationVulnerability(BaseModel):
    """
    Two-pathway convergence creating emergent SL that neither pathway alone produces.

    Examples:
      BER deficiency (MBD4-LOF) × checkpoint inhibition (ATR/WEE1):
        Neither pathway alone reaches Strong, but combined they create mitotic catastrophe.
      Chromatin remodeling loss (ARID1A-LOF) × ATR checkpoint:
        REALIZED — CRISPR data present, combined tier = Strong.
    """
    pathway_a: str                     # e.g. "BER_deficiency"
    pathway_a_gene: str                # e.g. "MBD4"
    pathway_b: str                     # e.g. "checkpoint_inhibition"
    pathway_b_targets: List[str]       # e.g. ["ATR", "CHK1", "WEE1"]
    convergence_mechanism: str         # e.g. "mitotic_catastrophe"

    pathway_a_alone_tier: str          # What pathway A gets alone (e.g. "Mechanistic candidate only")
    pathway_b_alone_tier: str          # What pathway B gets alone (e.g. "Mechanistic candidate only")
    combined_tier: Optional[str] = None  # Projected combined tier (None = insufficient data)

    data_gaps: List[str] = Field(default_factory=list)
    # Exact data needed to unlock combinatorial tier. Empty = tier already realized.

    rs_score: Optional[float] = None   # RS score for this combination
    rs_level: Optional[str] = None     # RS level

    rationale: str                     # One-paragraph mechanistic chain


# ── Full evidence matrix ──────────────────────────────────────────────────────

class EvidenceMatrix(BaseModel):
    """Complete multi-modal evidence matrix for a query gene."""

    query_gene: str
    cancer_type: Optional[str] = None
    depmap_release: str = "unknown"
    rows: List[EvidenceRow] = Field(default_factory=list)
    calibration_axis: str = Field(
        "cytidine_analogs",
        description="The gold-standard axis used to define 'high evidence' threshold",
    )
    gold_standard_summary: Optional[str] = Field(
        None,
        description="One-paragraph description of the calibration axis evidence pattern",
    )
    notes: List[str] = Field(default_factory=list)
    rs_score: Optional[Any] = Field(
        None,
        description="RSScore object if rs_features were provided in the query input",
    )
    combination_vulnerabilities: List["CombinationVulnerability"] = Field(
        default_factory=list,
        description="Combinatorial two-pathway vulnerabilities detected for this gene/cancer context",
    )

    def get_row(self, axis: CandidateAxis) -> Optional[EvidenceRow]:
        for r in self.rows:
            if r.axis == axis:
                return r
        return None

    def recommendation_summary(self) -> Dict[str, str]:
        """Quick lookup: axis_label → recommendation_tier."""
        return {r.axis_label: (r.recommendation_tier or "Unknown") for r in self.rows}


# ── Multi-modal query input ───────────────────────────────────────────────────

class MultiModalQueryInput(BaseModel):
    """Input to the multi-modal analysis endpoint."""

    gene: str = Field(..., description="Hugo gene symbol (e.g. MBD4)")
    mutation_type: str = Field("loss_of_function")
    cancer_type: Optional[str] = None
    axes: List[CandidateAxis] = Field(
        default_factory=lambda: list(CandidateAxis),
        description="Which candidate axes to evaluate. Defaults to all.",
    )
    include_pharmacologic_stratification: bool = True
    include_literature_receipts: bool = True
    include_calibration_narrative: bool = True
    # Stratification sub-filters
    stratify_by_msi: bool = Field(
        True,
        description="Separately stratify pharmacologic screens by MSI-H vs MSS",
    )
    # RS context — optional patient replication stress profile
    rs_features: Optional[Any] = Field(
        None,
        description=(
            "RSFeatureSet: patient genomic features for Replication Stress scoring. "
            "When provided, RS score is computed and used to modulate tier assignment "
            "for ATR/WEE1 and PKMYT1 axes only. "
            "Import RSFeatureSet from sl_agent.multimodal.replication_stress."
        ),
    )


# ── Multi-modal analysis result ───────────────────────────────────────────────

class ModalityAgreementReport(BaseModel):
    """Structured report on where CRISPR and pharmacologic screens agree/disagree."""
    axis_label: str
    agreement_category: str
    crispr_summary: str
    pharmacologic_summary: str
    interpretation_change: bool = Field(
        False,
        description="True if the disagreement materially changes our conclusion",
    )
    conclusion: str


class MultiModalResult(BaseModel):
    """Top-level output of the multi-modal analysis pipeline."""

    query_gene: str
    cancer_type: Optional[str]
    evidence_matrix: EvidenceMatrix
    agreement_report: List[ModalityAgreementReport]
    updated_narrative: str = Field(
        ...,
        description="Updated interpretation referencing all modalities explicitly",
    )
    recommendation_map: Dict[str, str] = Field(
        default_factory=dict,
        description="axis_label → recommendation_tier mapping",
    )
    guardrails_applied: List[str] = Field(
        default_factory=list,
        description="List of guardrail rules that were triggered during analysis",
    )
    warnings: List[str] = Field(default_factory=list)
