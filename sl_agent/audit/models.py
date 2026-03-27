"""
Audit queue data models for the Receipt Miner pipeline.

All models include RUO fields — this layer is the backend of the doctor portal.
"""
from __future__ import annotations
from datetime import datetime
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field


RUO_DISCLAIMER = (
    "Research use only. Not validated for clinical decision-making. "
    "All receipt candidates require human expert review before promotion."
)


class RUOEnvelope(BaseModel):
    """Standard wrapper for all audit API responses."""
    ruo: bool = True
    disclaimer: str = RUO_DISCLAIMER
    data: Any


class ReceiptCandidate(BaseModel):
    """
    Auto-generated evidence candidate awaiting human audit before promotion
    to _FROZEN_RECEIPTS.

    Confidence score components (weights sum to 1.0):
      crispr:     0.35  (delta_dep significance in mutant vs WT)
      pharma:     0.30  (PRISM/GDSC drug screen stratification)
      kb_clinical: 0.25 (CIViC/CGI/ClinVar hits with clinical evidence)
      expression:  0.10 (expression correlation)

    Promotion requires: audit_status == "approved" AND human sets promoted_to_frozen = True.
    Auto-generation NEVER promotes directly. Human approval is mandatory.
    """
    id: Optional[int] = None              # set by DB on insert
    gene: str
    axis: str                             # CandidateAxis.value string
    cancer_type: Optional[str] = None

    # Evidence signals (raw — used to compute confidence_score)
    crispr_delta_dep: Optional[float] = None
    crispr_fdr: Optional[float] = None
    crispr_n_mutant: Optional[int] = None
    crispr_n_wt: Optional[int] = None

    prism_delta_auc: Optional[float] = None
    gdsc_delta_ic50: Optional[float] = None
    kb_clinical_hits: int = 0            # CIViC/CGI/ClinVar hits with clinical evidence
    expression_corr: Optional[float] = None

    # Computed
    confidence_score: float = Field(..., ge=0.0, le=1.0)
    candidate_tier: str                  # projected tier if promoted: Validated/Strong/Mechanistic
    evidence_summary: str                # one-line auto-generated summary

    # Audit lifecycle
    audit_status: str = "pending"        # "pending" | "approved" | "rejected"
    audit_notes: Optional[str] = None
    audited_by: Optional[str] = None
    audited_at: Optional[datetime] = None
    generated_at: datetime = Field(default_factory=datetime.utcnow)
    promoted_to_frozen: bool = False

    # Provenance
    depmap_release: str = "unknown"
    source_pipeline: str = "receipt_miner_v1"


class AuditAction(BaseModel):
    """Request body for approve/reject endpoints."""
    notes: Optional[str] = None
    audited_by: str = Field(..., description="Auditor identifier (email or name)")


class CoverageEntry(BaseModel):
    """
    One cell in the gene × axis coverage grid.
    Feeds the Atlas Panel 1 heat map and doctor portal coverage view.
    """
    gene: str
    axis: str                            # CandidateAxis.value string
    axis_label: str                      # Human-readable label
    tier: str                            # Current tier from frozen receipts / matrix
    frozen_receipt_count: int            # How many frozen receipts exist for this pair
    candidate_count: int                 # Pending candidates in audit queue
    approved_count: int                  # Approved-but-not-yet-frozen candidates
    confidence: Optional[float]          # Best confidence score among pending candidates
    has_clinical_evidence: bool          # Any clinical POSITIVE in frozen receipts


class QueueStats(BaseModel):
    """Summary statistics for the audit queue."""
    total_pending: int
    total_approved: int
    total_rejected: int
    total_promoted: int
    high_confidence_pending: int         # confidence >= 0.70
    coverage_pairs_total: int            # gene × axis pairs tracked
    coverage_pairs_with_frozen: int      # pairs that have at least one frozen receipt
