"""
Audit queue API — backend for the doctor portal.

All endpoints:
  - Tagged "Audit / Doctor Portal"
  - Include RUO wrapper: {"ruo": true, "disclaimer": "...", "data": ...}
  - Designed to feed three doctor portal views:
      View 1: Case view (GET /audit/queue filtered by gene)
      View 2: Audit queue (GET /audit/queue, POST /audit/{id}/approve|reject)
      View 3: Coverage heat map (GET /audit/coverage)

Endpoints:
  GET  /audit/queue              — list pending candidates (sortable, filterable)
  GET  /audit/stats              — queue summary counts
  GET  /audit/coverage           — gene × axis coverage grid (for heat map)
  POST /audit/seed               — inject a test candidate (dev/test only)
  GET  /audit/{id}               — single candidate detail
  POST /audit/{id}/approve       — mark approved (does NOT auto-promote)
  POST /audit/{id}/reject        — mark rejected

NOTE: Static paths (/queue, /stats, /coverage, /seed) MUST be declared BEFORE
      the dynamic /{candidate_id} route to avoid path conflicts.
"""
from __future__ import annotations

import logging
from typing import List, Optional

from fastapi import APIRouter, BackgroundTasks, HTTPException, Query, Response, status

from .models import (
    AuditAction,
    CoverageEntry,
    QueueStats,
    ReceiptCandidate,
    RUO_DISCLAIMER,
    RUOEnvelope,
)
from .queue import AuditQueue, get_coverage_version, init_db, set_last_mine_at

logger = logging.getLogger(__name__)

audit_router = APIRouter(
    prefix="/audit",
    tags=["Audit / Doctor Portal"],
)

# Initialize DB on first import
try:
    init_db()
except Exception as e:
    logger.warning("Audit DB init deferred: %s", e)


def _wrap(data) -> dict:
    """Wrap response in RUO envelope."""
    return {"ruo": True, "disclaimer": RUO_DISCLAIMER, "data": data}


# ── GET /audit/queue ──────────────────────────────────────────────────────────
# MUST come before /{candidate_id}

@audit_router.get(
    "/queue",
    summary="List pending receipt candidates sorted by confidence",
    description=(
        "Returns pending candidates for human review. "
        "Feeds Doctor Portal View 2 (Audit Queue). "
        "Filter by gene for View 1 (Case View)."
    ),
)
async def list_queue(
    gene: Optional[str] = Query(None, description="Filter by gene (e.g. MBD4)"),
    axis: Optional[str] = Query(None, description="Filter by axis (e.g. atr_wee1)"),
    min_confidence: float = Query(0.0, ge=0.0, le=1.0, description="Minimum confidence score"),
    limit: int = Query(100, ge=1, le=500),
):
    candidates = AuditQueue.list_pending(
        min_confidence=min_confidence,
        gene=gene,
        axis=axis,
        limit=limit,
    )
    return _wrap([c.model_dump() for c in candidates])


# ── GET /audit/stats ──────────────────────────────────────────────────────────
# MUST come before /{candidate_id}

@audit_router.get(
    "/stats",
    summary="Audit queue summary statistics",
)
async def queue_stats():
    s = AuditQueue.stats()
    return _wrap(s)


# ── GET /audit/coverage ───────────────────────────────────────────────────────
# MUST come before /{candidate_id}

@audit_router.get(
    "/coverage",
    summary="Gene × axis coverage grid",
    description=(
        "Returns all (gene, axis) pairs with their current tier, frozen receipt count, "
        "pending candidate count, and best confidence score. "
        "Feeds the Atlas Panel 1 heat map and Doctor Portal View 3 (Coverage). "
        "Updates automatically as candidates are approved."
    ),
)
async def coverage_grid(response: Response):
    result = AuditQueue.coverage()  # now returns {coverage, summary, coverage_version}
    # Cache headers — doctor portal polls this endpoint and compares X-Coverage-Version
    # to detect when the heat map needs a refresh (approve/reject/mine all bump the version).
    response.headers["Cache-Control"] = "no-cache"
    response.headers["X-Coverage-Version"] = get_coverage_version()
    # Merge coverage_version into response body too (for clients that can't read headers)
    return {
        "ruo": True,
        "disclaimer": RUO_DISCLAIMER,
        "coverage": result["coverage"],
        "summary": result["summary"],
        "coverage_version": result["coverage_version"],
    }


# ── POST /audit/seed (dev/test only) ─────────────────────────────────────────
# MUST come before /{candidate_id}

@audit_router.post(
    "/seed",
    summary="Seed a test candidate (dev/test only)",
    description="Injects a ReceiptCandidate directly into the queue. Dev/test use only.",
    include_in_schema=True,  # set False in production
)
async def seed_candidate(candidate: ReceiptCandidate):
    saved = AuditQueue.upsert(candidate)
    return _wrap(saved.model_dump())


# ── POST /audit/mine ─────────────────────────────────────────────────────────
# MUST come before /{candidate_id} — route ordering is critical

@audit_router.post(
    "/mine",
    summary="Trigger receipt mining for a gene panel",
    description=(
        "Runs mine_receipts() as a background task across the provided gene panel × axes. "
        "Returns immediately. Check GET /audit/queue for results after completion. "
        "VALIDATED tier is never auto-assigned — all promotions require human approval."
    ),
)
async def trigger_mine(
    background_tasks: BackgroundTasks,
    gene_panel: Optional[List[str]] = Query(None, description="Gene symbols (e.g. BRCA1,BRCA2)"),
    axes: Optional[List[str]] = Query(None, description="Axis values (e.g. parp_inhibitors,atr_wee1)"),
):
    from ..multimodal.receipt_miner import (
        mine_receipts, MinerThresholds, DEFAULT_GENE_PANEL, DEFAULT_AXES,
        load_default_panel,
    )
    from ..multimodal.models import CandidateAxis

    # When no gene_panel provided, use load_default_panel() which reads
    # ddr_panel_v1.json and deduplicates. Falls back to DEFAULT_GENE_PANEL.
    panel = gene_panel if gene_panel else load_default_panel()
    axis_enums = (
        [CandidateAxis(a) for a in axes if a in CandidateAxis._value2member_map_]
        if axes else DEFAULT_AXES
    )

    estimated_pairs = len(panel) * len(axis_enums)

    def _run_mine():
        try:
            summary = mine_receipts(gene_panel=panel, axes=axis_enums, thresholds=MinerThresholds())
            # Pre-warm coverage cache after mine completes + stamp new version
            from .queue import _set_coverage_version
            _set_coverage_version()
            set_last_mine_at()
            AuditQueue.coverage()  # pre-warm — result discarded but data is hot for next GET
            logger.info(
                "Mine complete: %d candidates queued, coverage refreshed",
                summary.candidates_queued,
            )
        except Exception as e:
            logger.error("mine_receipts background task failed: %s", e)

    background_tasks.add_task(_run_mine)

    return _wrap({
        "status": "running",
        "gene_count": len(panel),
        "axis_count": len(axis_enums),
        "estimated_pairs": estimated_pairs,
        "message": "Mining started. Check GET /audit/queue for results.",
    })


# ── GET /audit/{candidate_id} ─────────────────────────────────────────────────
# Dynamic route — MUST come AFTER all static routes

@audit_router.get(
    "/{candidate_id}",
    summary="Get single candidate detail",
)
async def get_candidate(candidate_id: int):
    # INTEGER OVERFLOW GUARD: SQLite INTEGER is 64-bit signed; values > 2^63-1 raise
    # OverflowError which propagates as a 500. Reject oversized IDs with 422.
    if candidate_id > 2**63 - 1 or candidate_id < 1:
        raise HTTPException(status_code=422, detail="candidate_id out of valid range")
    c = AuditQueue.get(candidate_id)
    if not c:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    return _wrap(c.model_dump())


# ── POST /audit/{candidate_id}/approve ───────────────────────────────────────

@audit_router.post(
    "/{candidate_id}/approve",
    summary="Approve a receipt candidate",
    description=(
        "Marks candidate as approved. Does NOT auto-promote to frozen receipts — "
        "human must explicitly paste the returned frozen_receipt_snippet into "
        "literature_receipts.py. Once approved, the candidate will not be "
        "overwritten by future mining runs."
    ),
)
async def approve_candidate(candidate_id: int, action: AuditAction):
    if candidate_id > 2**63 - 1 or candidate_id < 1:
        raise HTTPException(status_code=422, detail="candidate_id out of valid range")
    try:
        result = AuditQueue.approve(candidate_id, action.notes, action.audited_by)
    except ValueError as e:
        # Immutability violation (already rejected) or empty audited_by
        raise HTTPException(status_code=409, detail=str(e))
    if not result:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    candidate = result["candidate"]
    snippet = result["frozen_receipt_snippet"]
    return _wrap({
        "status": "approved",
        "candidate_id": candidate_id,
        "candidate": candidate.model_dump(),
        "frozen_receipt_snippet": snippet,
    })


# ── POST /audit/{candidate_id}/reject ────────────────────────────────────────

@audit_router.post(
    "/{candidate_id}/reject",
    summary="Reject a receipt candidate",
    description="Marks candidate as rejected with auditor notes.",
)
async def reject_candidate(candidate_id: int, action: AuditAction):
    if candidate_id > 2**63 - 1 or candidate_id < 1:
        raise HTTPException(status_code=422, detail="candidate_id out of valid range")
    try:
        updated = AuditQueue.reject(candidate_id, action.notes, action.audited_by)
    except ValueError as e:
        # Immutability violation (already approved) or empty audited_by
        raise HTTPException(status_code=409, detail=str(e))
    if not updated:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    return _wrap(updated.model_dump())
