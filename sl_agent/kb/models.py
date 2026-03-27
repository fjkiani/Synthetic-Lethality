"""
Unified KB domain models — the canonical schema all sources normalise into.

Inspired by the VICC meta-knowledgebase harmonization model.
Every source (CIViC, CGI, JAX, PMKB, ClinVar) produces KBEvidence records.
The tiering engine then reduces them to EvidenceTier + drug recommendations.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


# ── Evidence tier (Mars-native, maps from all source levels) ──────────────────
#
#  Tier A  — FDA-approved or guideline-recommended (multi-KB consensus)
#  Tier B  — Clinical evidence (≥1 Phase II/III trial, multiple KBs)
#  Tier C  — Investigational / Phase I, or single-KB only
#  Tier D  — Preclinical / case report / computational only
#  Unknown — Not enough information to tier

class EvidenceTier(str, Enum):
    A = "Tier_A"   # Approved / standard of care
    B = "Tier_B"   # Clinical trial evidence
    C = "Tier_C"   # Investigational / exploratory
    D = "Tier_D"   # Preclinical / in-vitro
    UNKNOWN = "Unknown"


class ResponseType(str, Enum):
    SENSITIVITY = "sensitivity"
    RESISTANCE = "resistance"
    REDUCED_SENSITIVITY = "reduced_sensitivity"
    ADVERSE = "adverse"
    UNKNOWN = "unknown"


class SourceKB(str, Enum):
    CIVIC = "CIViC"
    CGI = "CGI"
    JAX = "JAX-CKB"
    PMKB = "PMKB"
    CLINVAR = "ClinVar"
    CHEMBL = "ChEMBL"
    GDSC = "GDSC"
    INTERNAL = "Internal"


# ── Core evidence record ──────────────────────────────────────────────────────

class KBEvidence(BaseModel):
    """
    Single harmonized evidence record from any source KB.
    All fields use Mars-native names; source_raw holds the original payload.
    """
    # Identity
    source_kb: SourceKB
    source_id: Optional[str] = None       # e.g. CIViC EID12345

    # Genomic context
    gene: str
    variant: Optional[str] = None         # e.g. "V600E", "exon 19 del", "LOF"
    variant_type: Optional[str] = None    # "missense", "deletion", "amplification", etc.

    # Disease context
    disease: Optional[str] = None
    disease_ontology_id: Optional[str] = None  # DOID or NCIt

    # Therapeutic
    drug: str
    drug_aliases: List[str] = Field(default_factory=list)
    drug_class: Optional[str] = None      # e.g. "PARP inhibitor"
    response_type: ResponseType = ResponseType.UNKNOWN

    # Evidence quality
    evidence_level_raw: Optional[str] = None   # source-native level string
    tier: EvidenceTier = EvidenceTier.UNKNOWN
    num_sources: int = 1                         # incremented when multiple KBs agree

    # Human-readable payload
    evidence_statement: Optional[str] = None
    pmids: List[str] = Field(default_factory=list)
    clinical_trials: List[str] = Field(default_factory=list)

    # Provenance
    source_raw: Optional[Dict[str, Any]] = None
    source_url: Optional[str] = None


# ── Aggregated drug recommendation ───────────────────────────────────────────

class DrugRecommendation(BaseModel):
    """
    Final drug recommendation after evidence fusion across KBs.
    This is what /kb/query returns.
    """
    gene: str
    variant: Optional[str] = None
    drug: str
    drug_class: Optional[str] = None
    response_type: ResponseType
    tier: EvidenceTier
    confidence: str   # "high" | "medium" | "low" based on KB consensus

    # Evidence backing
    supporting_kbs: List[SourceKB]
    num_evidence_items: int
    pmids: List[str] = Field(default_factory=list)
    clinical_trials: List[str] = Field(default_factory=list)

    # Explain layer
    evidence_summaries: List[str] = Field(default_factory=list)
    cancer_types: List[str] = Field(default_factory=list)
    notes: Optional[str] = None


# ── Query / Response models ───────────────────────────────────────────────────

class KBQueryInput(BaseModel):
    gene: str = Field(..., description="Hugo gene symbol")
    variant: Optional[str] = Field(None, description="e.g. 'V600E', 'amplification', 'LOF'")
    cancer_type: Optional[str] = Field(None, description="e.g. 'Ovarian Cancer', 'NSCLC'")
    min_tier: EvidenceTier = Field(EvidenceTier.D, description="Minimum tier to include")
    include_resistance: bool = True
    top_n: int = Field(20, ge=1, le=100)


class KBQueryResponse(BaseModel):
    query: KBQueryInput
    recommendations: List[DrugRecommendation]
    sources_queried: List[SourceKB]
    total_evidence_items: int
    notes: Optional[str] = None


class KBExplainResponse(BaseModel):
    gene: str
    drug: str
    raw_evidence: List[KBEvidence]
    recommendation: Optional[DrugRecommendation] = None
