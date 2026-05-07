"""
BrM v2 Pydantic models.

Key additions over v1:
- BrMRowClass: CALIBRATION | NOVEL | NEGATIVE (3-section partitioning)
- ScoreBasis: RECEIPT | DEPMAP | HYBRID (score provenance)
- SLPartner: structured DepMap SL discovery result
- emt_cluster: bool flag (POL-001 annotation-only)
- broad_state_marker: bool flag (POL-005 annotation-only)
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel, Field, field_validator


RUO_DISCLAIMER = (
    "Research use only. Not validated for clinical decision-making. "
    "All BrM targetability rows require human expert review before clinical translation."
)


# ---------------------------------------------------------------------------
# v2 enums
# ---------------------------------------------------------------------------

class BrMRowClass(str, Enum):
    CALIBRATION = "calibration"   # receipt-backed, literature-validated
    NOVEL       = "novel"         # DepMap-discovered or receipt-only without DepMap signal
    NEGATIVE    = "negative"      # tested, no signal at FDR<0.25


class ScoreBasis(str, Enum):
    RECEIPT = "receipt"    # score derived from frozen literature receipts
    DEPMAP  = "depmap"     # score derived from DepMap expression-stratified SL
    HYBRID  = "hybrid"     # both sources contributed


class ExploitMode(str, Enum):
    DIRECT   = "direct"
    PATHWAY  = "pathway"
    DELIVERY = "delivery"
    NO_PATH  = "no_current_path"


class BrMEvidenceClass(str, Enum):
    IN_VIVO_CRISPR_SCREEN = "in_vivo_crispr_screen"
    PATIENT_EXPRESSION    = "patient_expression"
    GENETIC_SUPPRESSION   = "genetic_suppression"
    PHARMACOLOGIC         = "pharmacologic"
    PATHWAY_INFERENCE     = "pathway_inference"
    BBB_BIOLOGY           = "bbb_biology"


class CascadeStep(str, Enum):
    INTRAVASATION   = "intravasation"
    CIRCULATION     = "circulation"
    BBB_TRANSIT     = "bbb_transit"
    COLONIZATION    = "colonization"
    IMMUNE_EVASION  = "immune_evasion"
    DORMANCY_ESCAPE = "dormancy_escape"
    PROLIFERATION   = "proliferation"
    MULTI_STEP      = "multi_step"


class TargetNodeRole(str, Enum):
    SELF               = "self"
    PARTNER            = "partner"
    PATHWAY_BOTTLENECK = "pathway_bottleneck"
    MICROENVIRONMENT   = "microenvironment"
    DELIVERY_INTERFACE = "delivery_interface"


class ContextSpecificity(str, Enum):
    BRM_SPECIFIC = "brm_specific"
    CNS_RELEVANT = "cns_relevant"
    PAN_CANCER   = "pan_cancer"


# ---------------------------------------------------------------------------
# SL partner (DepMap discovery)
# ---------------------------------------------------------------------------

class SLPartner(BaseModel):
    """A single SL partner discovered via expression-stratified DepMap analysis."""
    partner_gene: str
    delta_dep: float          # Q25_high - Q25_low dependency; negative = SL
    fdr: float
    n_high: int               # lines in high-expression quartile
    n_low: int                # lines in low-expression quartile
    partner_quality_tier: str  # "strong" | "moderate" | "weak"

    @field_validator("fdr")
    @classmethod
    def fdr_in_range(cls, v: float) -> float:
        if not (0.0 <= v <= 1.0):
            raise ValueError(f"fdr must be in [0, 1], got {v}")
        return v

    @field_validator("partner_quality_tier")
    @classmethod
    def valid_tier(cls, v: str) -> str:
        if v not in ("strong", "moderate", "weak"):
            raise ValueError(f"partner_quality_tier must be strong/moderate/weak, got {v!r}")
        return v


# ---------------------------------------------------------------------------
# Core row model
# ---------------------------------------------------------------------------

class BrMGeneRow(BaseModel):
    # v2 classification
    row_class: BrMRowClass
    score_basis: ScoreBasis

    # Identity
    gene: str
    visibility_axis: str

    # BrM evidence provenance
    brm_evidence_class: BrMEvidenceClass
    primary_evidence_modality: str
    secondary_evidence_modalities: List[str] = Field(default_factory=list)

    # Metastatic cascade
    cascade_step: CascadeStep
    cascade_step_label: str

    # Exploit routing
    primary_exploit_mode: ExploitMode
    secondary_exploit_modes: List[ExploitMode] = Field(default_factory=list)
    primary_vulnerability_hypothesis: str

    # Target specification
    target_node: str
    target_node_role: TargetNodeRole
    drug_class: str
    best_drug_example: Optional[str] = None

    # CNS/BBB feasibility
    cns_bbb_feasibility: str
    cns_bbb_rationale: str

    # Context
    context_specificity: ContextSpecificity

    # Evidence tier and recommendation
    evidence_tier: str
    recommendation_status: str
    confidence_score: float

    # v2 novelty score (DepMap-derived; None for CALIBRATION rows)
    novelty_score: Optional[float] = None

    # v2 SL partners (DepMap-discovered; empty for CALIBRATION/NEGATIVE)
    sl_partners: List[SLPartner] = Field(default_factory=list)

    # v2 annotation flags (POL-001, POL-005)
    emt_cluster: bool = False          # POL-001: EMT co-expression cluster member
    broad_state_marker: bool = False   # POL-005: broad state marker, zero FDR<0.10 partners

    # Provenance
    supporting_pmids: List[str] = Field(default_factory=list)
    depmap_expression_note: Optional[str] = None
    live_literature_note: Optional[str] = None

    # Audit trail
    escalation_path: List[str] = Field(default_factory=list)
    data_gaps: List[str] = Field(default_factory=list)
    frozen_at: str

    @field_validator("confidence_score")
    @classmethod
    def confidence_in_range(cls, v: float) -> float:
        if not (0.0 <= v <= 1.0):
            raise ValueError(f"confidence_score must be in [0, 1], got {v}")
        return v

    @field_validator("evidence_tier")
    @classmethod
    def no_validated_tier(cls, v: str) -> str:
        if v == "Validated":
            raise ValueError(
                "evidence_tier 'Validated' is never auto-assigned. "
                "Use 'Strong', 'Mechanistic', or 'Insufficient'."
            )
        return v

    @field_validator("secondary_exploit_modes")
    @classmethod
    def secondary_differs_from_primary(cls, v: List[ExploitMode], info) -> List[ExploitMode]:
        # Validated after primary_exploit_mode is set
        return v


# ---------------------------------------------------------------------------
# Matrix model
# ---------------------------------------------------------------------------

class BrMTargetabilityMatrix(BaseModel):
    query_context: str
    cancer_type: str
    depmap_release: str = "24Q4"
    gold_standard_gene: str = "BACE1"
    version: str = "v2"

    rows: List[BrMGeneRow]

    # 3-section partitioning
    calibration_rows: List[BrMGeneRow] = Field(default_factory=list)
    novel_rows: List[BrMGeneRow] = Field(default_factory=list)
    negative_rows: List[BrMGeneRow] = Field(default_factory=list)

    exploit_mode_summary: Dict[str, int] = Field(default_factory=dict)
    cascade_step_summary: Dict[str, int] = Field(default_factory=dict)
    row_class_summary: Dict[str, int] = Field(default_factory=dict)

    frozen_at: str
    ruo: bool = True
    disclaimer: str = RUO_DISCLAIMER


# ---------------------------------------------------------------------------
# Query input
# ---------------------------------------------------------------------------

class BrMQueryInput(BaseModel):
    genes: List[str] = Field(
        default_factory=lambda: [
            "BACE1", "MMP9", "CLDN5", "CCL2", "ICAM1", "MMP2", "TWIST1"
        ]
    )
    cancer_type: str = "NSCLC"
    context: str = "brain_metastasis"
    include_depmap_expression: bool = True
    include_live_literature: bool = True
