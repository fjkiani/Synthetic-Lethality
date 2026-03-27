"""
Pydantic v2 domain models — the contract for every layer.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


# ── Enums ─────────────────────────────────────────────────────────────────────

class AnalysisScope(str, Enum):
    CANCER_SPECIFIC = "cancer_specific"
    PAN_CANCER = "pan_cancer"


class DrugClass(str, Enum):
    PARP_INHIBITOR = "PARP inhibitor"
    ATR_INHIBITOR = "ATR inhibitor"
    ATM_INHIBITOR = "ATM inhibitor"
    WEE1_INHIBITOR = "WEE1 inhibitor"
    CHK1_INHIBITOR = "CHK1 inhibitor"
    KINASE_INHIBITOR = "kinase inhibitor"
    EPIGENETIC = "epigenetic modifier"
    PROTEASOME = "proteasome inhibitor"
    TOPOISOMERASE = "topoisomerase inhibitor"
    ANTIMETABOLITE = "antimetabolite"
    ALKYLATING = "alkylating agent"
    IMMUNOTHERAPY = "immunotherapy"
    OTHER = "other"
    UNKNOWN = "unknown"


class MutationType(str, Enum):
    LOF = "loss_of_function"
    GOF = "gain_of_function"
    HOMDEL = "homozygous_deletion"
    AMP = "amplification"
    ANY = "any_mutation"


# ── Input ─────────────────────────────────────────────────────────────────────

class SLQueryInput(BaseModel):
    """Primary query payload sent to the agent."""

    gene: str = Field(..., description="Hugo gene symbol of the driver / query gene")
    mutation_type: MutationType = Field(
        MutationType.ANY,
        description="Which mutation class to use when stratifying cell lines",
    )
    cancer_type: Optional[str] = Field(
        None,
        description="Lineage or subtype (e.g. 'Ovarian Cancer', 'NSCLC'). "
                    "Falls back to pan-cancer if < 30 lines available.",
    )
    top_n_partners: int = Field(
        20, ge=1, le=100, description="Max SL partner genes to return"
    )
    top_n_drugs: int = Field(
        10, ge=1, le=50, description="Max gene–drug pairs per partner"
    )
    fdr_cutoff: float = Field(0.25, ge=0.0, le=1.0)
    delta_dep_cutoff: float = Field(0.15, ge=0.0, le=5.0)
    include_pathway_context: bool = True
    include_codependency: bool = True


# ── SL partner ────────────────────────────────────────────────────────────────

class FrameworkSupport(BaseModel):
    SLIdR: Optional[bool] = None
    SLAYER: Optional[bool] = None
    SL_RFM: Optional[bool] = None
    literature: Optional[bool] = None


class SLPartner(BaseModel):
    gene: str
    delta_dependency: float = Field(..., description="mutant_mean − wt_mean (Chronos)")
    p_value: float
    fdr: float
    n_mut: int
    n_wt: int
    test_type: str = Field(..., description="e.g. 'wilcoxon', 'ttest'")
    effect_size_cohend: Optional[float] = None
    codependency_r: Optional[float] = None     # Pearson r of dependency profiles
    pathway: Optional[List[str]] = None
    supporting_frameworks: FrameworkSupport = Field(default_factory=FrameworkSupport)
    is_essential_in_wt: bool = Field(
        False, description="True if gene is pan-essential (confound flag)"
    )


# ── Drug evidence ─────────────────────────────────────────────────────────────

class DepmapDrugEvidence(BaseModel):
    dataset: str                   # e.g. "PRISM_repurposing", "GDSC2"
    drug_name: str
    delta_viability: Optional[float] = None   # mutant_mean − wt_mean AUC/LFC
    effect_size: Optional[float] = None
    p_value: Optional[float] = None
    fdr: Optional[float] = None
    n_mut: Optional[int] = None
    n_wt: Optional[int] = None


class ExternalDrugEvidence(BaseModel):
    source: str             # "OncoKB", "ChEMBL", "DrugBank", "GDSC_annotation"
    drug_name: str
    drug_chembl_id: Optional[str] = None
    mechanism_of_action: Optional[str] = None
    drug_class: DrugClass = DrugClass.UNKNOWN
    primary_targets: Optional[List[str]] = None
    max_phase: Optional[int] = Field(
        None, description="0=preclinical, 1-3=trials, 4=approved"
    )
    oncokb_level: Optional[str] = None        # e.g. "LEVEL_1", "LEVEL_3A"
    notes: Optional[str] = None


class GeneDrugPair(BaseModel):
    partner_gene: str
    drug_name: str
    drug_class: DrugClass = DrugClass.UNKNOWN
    depmap_evidence: Optional[DepmapDrugEvidence] = None
    external_evidence: Optional[ExternalDrugEvidence] = None
    rank_score: float = Field(
        0.0,
        description="Composite: SL signal strength × drug response differential × druggability",
    )
    rank_components: Dict[str, float] = Field(default_factory=dict)


# ── Output ────────────────────────────────────────────────────────────────────

class InputContext(BaseModel):
    query_gene: str
    mutation_type: MutationType
    cancer_type: Optional[str]
    scope: AnalysisScope
    n_total_lines: int
    n_mut_lines: int
    n_wt_lines: int
    depmap_release: str
    mutation_frequency_depmap: Optional[float] = None
    mutation_frequency_tcga: Optional[float] = None
    notes: Optional[str] = None


class CrossValidation(BaseModel):
    frameworks_checked: List[str]
    confirmed_pairs: List[str]   # partner gene symbols confirmed by ≥1 framework
    speculative_pairs: List[str]
    notes: Optional[str] = None


class SLMapResult(BaseModel):
    """Top-level output — machine-readable JSON / structured tables."""

    input_context: InputContext
    sl_partners: List[SLPartner]
    gene_drug_pairs: List[GeneDrugPair]
    cross_validation: CrossValidation
    notes: Optional[str] = None
    warnings: List[str] = Field(default_factory=list)

    # ── Convenience accessor ──────────────────────────────────────────────
    def top_partners(self, n: int = 5) -> List[SLPartner]:
        return sorted(self.sl_partners, key=lambda x: x.fdr)[:n]

    def top_drug_pairs(self, n: int = 5) -> List[GeneDrugPair]:
        return sorted(self.gene_drug_pairs, key=lambda x: -x.rank_score)[:n]


# ── API transport ─────────────────────────────────────────────────────────────

class SLQueryResponse(BaseModel):
    status: str = "success"
    result: SLMapResult
    processing_time_ms: float


class ErrorResponse(BaseModel):
    status: str = "error"
    error: str
    detail: Optional[Any] = None
