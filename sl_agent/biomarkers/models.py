"""
Biomarker subsystem data models.

Mirrors the pydantic style of sl_agent/multimodal/models.py.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field


# ── Biomarker capability type ─────────────────────────────────────────────────

class BiomarkerType(str, Enum):
    PROGNOSTIC = "prognostic"     # predicts patient outcome (OS/PFS)
    DIAGNOSTIC = "diagnostic"     # classifies disease state / molecular subtype
    PREDICTIVE = "predictive"     # predicts treatment response (kept w/ SL engine)


# ── Validation status (the load-bearing honesty gate) ─────────────────────────

class BiomarkerStatus(str, Enum):
    VALIDATED       = "validated"        # reproduced in an INDEPENDENT, non-overlapping cohort
    DISCOVERY_ONLY  = "discovery_only"   # trained/derived; external replication absent or failed
    FAILED          = "failed"           # did not reproduce; retained for provenance only


# ── A single validation receipt ───────────────────────────────────────────────

class ValidationReceipt(BaseModel):
    """
    One record proving (or refuting) a biomarker on a specific cohort.
    Every metric a biomarker ships MUST trace to a receipt like this.
    """
    kind: str = Field(..., description="'train' | 'internal_cv' | 'external_replication' | 'sensitivity_check'")
    cohort_id: str = Field(..., description="Exact dataset/accession + version")
    n: int
    n_events: Optional[int] = Field(None, description="Events for survival (deaths/progressions)")
    is_independent_of_training: bool = Field(
        ..., description="True only if NO patient/sample overlaps the training cohort"
    )
    # Metrics (fill those relevant to the biomarker type)
    metric_name: Optional[str] = Field(None, description="e.g. 'c_index', 'auroc', 'concordance', 'HR_per_1SD'")
    metric_value: Optional[float] = None
    ci95_low: Optional[float] = None
    ci95_high: Optional[float] = None
    p_value: Optional[float] = None
    ph_test_p: Optional[float] = Field(None, description="Schoenfeld PH test p (Cox only)")
    epv: Optional[float] = Field(None, description="events-per-variable (survival adequacy)")
    seed: Optional[int] = None
    method_detail: Optional[str] = None
    extra: Dict[str, Any] = Field(default_factory=dict)


# ── A registered biomarker model ──────────────────────────────────────────────

class BiomarkerModel(BaseModel):
    """
    A biomarker capability registered in the subsystem.
    Keyed by (cancer, biomarker_type, method_id).
    """
    cancer: str = Field(..., description="e.g. 'ovarian_hgsoc'")
    biomarker_type: BiomarkerType
    method_id: str = Field(..., description="Stable id, e.g. 'riester2014_prognostic'")

    # Provenance of the reused published method (no reinvention)
    source_pmid: Optional[str] = None
    source_doi: Optional[str] = None
    source_note: str = Field("", description="How the published method was reused")

    gene_panel: List[str] = Field(default_factory=list, description="HUGO symbols / feature set receipt")
    artifact_path: Optional[str] = Field(None, description="Path to the frozen trained model, if any")

    status: BiomarkerStatus = BiomarkerStatus.DISCOVERY_ONLY
    receipts: List[ValidationReceipt] = Field(default_factory=list)

    interpretation: str = Field("", description="Honest, receipt-grounded summary string")
    caveats: List[str] = Field(default_factory=list)
    labeling: str = "RESEARCH USE ONLY. Not validated for clinical decision-making."

    def external_replication(self) -> Optional[ValidationReceipt]:
        """Return the first passing external-replication receipt, if any."""
        for r in self.receipts:
            if r.kind == "external_replication" and r.is_independent_of_training:
                return r
        return None

    def is_validated(self) -> bool:
        """
        A model is 'validated' ONLY with a passing external-replication receipt
        on an independent cohort. The status flag alone is never trusted.
        """
        return self.status == BiomarkerStatus.VALIDATED and self.external_replication() is not None


# ── A per-sample / per-query biomarker result ─────────────────────────────────

class BiomarkerResult(BaseModel):
    model_ref: str = Field(..., description="method_id of the BiomarkerModel used")
    biomarker_type: BiomarkerType
    call: Optional[str] = Field(None, description="e.g. risk group / subtype label")
    score: Optional[float] = None
    confidence: Optional[float] = Field(None, description="margin / calibrated prob, if available")
    status: BiomarkerStatus = BiomarkerStatus.DISCOVERY_ONLY
    interpretation: str = ""


# ── Advisory context attached to the multimodal evidence surface ──────────────

class BiomarkerContext(BaseModel):
    """
    Advisory-only block attached to the evidence matrix. It NEVER changes drug
    routing; it exists to surface prognostic/diagnostic context alongside the
    synthetic-lethality recommendation.
    """
    cancer: str
    results: List[BiomarkerResult] = Field(default_factory=list)
    labeling: str = "RESEARCH USE ONLY. Not validated for clinical decision-making."
    advisory_only: bool = True
