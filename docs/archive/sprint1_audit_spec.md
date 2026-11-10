# Sprint 1 Spec — Audit Queue + Doctor Portal API Layer

## Context
- FastAPI backend, Python 3.10+, Pydantic v2, SQLite (stdlib only — no new deps)
- Existing routers: `router` (/api/v1), `kb_router` (/api/v1), `mm_router` (/api/v1/analyze)
- New router: `audit_router` mounted at `/api/v1/audit`
- All endpoints tagged RUO. Every response envelope includes `"ruo": true` and `"disclaimer"`.
- Existing pattern: see `multimodal_routes.py` for FastAPI style, `models.py` for Pydantic v2 style.

## Files to Create

```
sl_agent/
├── audit/
│   ├── __init__.py
│   ├── models.py          ← ReceiptCandidate, AuditAction, CoverageEntry, RUOEnvelope
│   ├── queue.py           ← SQLite-backed AuditQueue class
│   └── routes.py          ← audit_router (7 endpoints)
├── multimodal/
│   └── models.py          ← ADD ReceiptCandidate import alias (do not duplicate)
└── api/
    └── app.py             ← ADD audit_router registration
```

## NEW: sl_agent/audit/models.py

```python
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
```

## NEW: sl_agent/audit/queue.py

```python
"""
SQLite-backed audit queue for ReceiptCandidate lifecycle management.

Schema: single table `receipt_candidates` with all ReceiptCandidate fields.
DB file: configured via settings.audit_db_path (default: .cache/audit_queue.db)
Thread-safe: uses connection-per-call pattern (SQLite WAL mode).

No ORM — direct sqlite3 to keep zero new dependencies.
"""
from __future__ import annotations

import json
import logging
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from .models import ReceiptCandidate

logger = logging.getLogger(__name__)

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS receipt_candidates (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    gene              TEXT NOT NULL,
    axis              TEXT NOT NULL,
    cancer_type       TEXT,
    crispr_delta_dep  REAL,
    crispr_fdr        REAL,
    crispr_n_mutant   INTEGER,
    crispr_n_wt       INTEGER,
    prism_delta_auc   REAL,
    gdsc_delta_ic50   REAL,
    kb_clinical_hits  INTEGER DEFAULT 0,
    expression_corr   REAL,
    confidence_score  REAL NOT NULL,
    candidate_tier    TEXT NOT NULL,
    evidence_summary  TEXT NOT NULL,
    audit_status      TEXT NOT NULL DEFAULT 'pending',
    audit_notes       TEXT,
    audited_by        TEXT,
    audited_at        TEXT,
    generated_at      TEXT NOT NULL,
    promoted_to_frozen INTEGER NOT NULL DEFAULT 0,
    depmap_release    TEXT NOT NULL DEFAULT 'unknown',
    source_pipeline   TEXT NOT NULL DEFAULT 'receipt_miner_v1',
    UNIQUE(gene, axis, cancer_type)  -- prevent duplicate candidates
)
"""

_PRAGMA_WAL = "PRAGMA journal_mode=WAL"


def _get_db_path() -> Path:
    from sl_agent.core.config import get_settings
    cfg = get_settings()
    path = getattr(cfg, "audit_db_path", Path(".cache/audit_queue.db"))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _conn() -> sqlite3.Connection:
    db = _get_db_path()
    con = sqlite3.connect(str(db), check_same_thread=False)
    con.row_factory = sqlite3.Row
    con.execute(_PRAGMA_WAL)
    return con


def init_db() -> None:
    """Create table if it doesn't exist. Called at app startup."""
    with _conn() as con:
        con.execute(_CREATE_TABLE)
        con.commit()
    logger.info("Audit queue DB initialized at %s", _get_db_path())


def _row_to_candidate(row: sqlite3.Row) -> ReceiptCandidate:
    d = dict(row)
    d["promoted_to_frozen"] = bool(d["promoted_to_frozen"])
    if d.get("audited_at"):
        d["audited_at"] = datetime.fromisoformat(d["audited_at"])
    if d.get("generated_at"):
        d["generated_at"] = datetime.fromisoformat(d["generated_at"])
    return ReceiptCandidate(**d)


class AuditQueue:
    """Thread-safe audit queue operations (stateless — each call opens its own connection)."""

    @staticmethod
    def upsert(candidate: ReceiptCandidate) -> ReceiptCandidate:
        """
        Insert or update a candidate. Uses (gene, axis, cancer_type) as unique key.
        If a pending candidate already exists for this pair, updates confidence/evidence.
        Approved/rejected candidates are NOT overwritten by new mining runs.
        Returns the candidate with its DB-assigned id.
        """
        with _conn() as con:
            # Check if approved/rejected — don't overwrite human decisions
            existing = con.execute(
                "SELECT id, audit_status FROM receipt_candidates "
                "WHERE gene=? AND axis=? AND (cancer_type IS ? OR cancer_type=?)",
                (candidate.gene.upper(), candidate.axis,
                 candidate.cancer_type, candidate.cancer_type)
            ).fetchone()

            if existing and existing["audit_status"] in ("approved", "rejected"):
                # Return existing — don't stomp human decision
                return _row_to_candidate(
                    con.execute(
                        "SELECT * FROM receipt_candidates WHERE id=?",
                        (existing["id"],)
                    ).fetchone()
                )

            now = datetime.utcnow().isoformat()
            con.execute("""
                INSERT INTO receipt_candidates
                  (gene, axis, cancer_type, crispr_delta_dep, crispr_fdr,
                   crispr_n_mutant, crispr_n_wt, prism_delta_auc, gdsc_delta_ic50,
                   kb_clinical_hits, expression_corr, confidence_score,
                   candidate_tier, evidence_summary, audit_status,
                   generated_at, depmap_release, source_pipeline)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,'pending',?,?,?)
                ON CONFLICT(gene, axis, cancer_type) DO UPDATE SET
                  crispr_delta_dep=excluded.crispr_delta_dep,
                  crispr_fdr=excluded.crispr_fdr,
                  confidence_score=excluded.confidence_score,
                  candidate_tier=excluded.candidate_tier,
                  evidence_summary=excluded.evidence_summary,
                  generated_at=excluded.generated_at
                WHERE receipt_candidates.audit_status = 'pending'
            """, (
                candidate.gene.upper(), candidate.axis, candidate.cancer_type,
                candidate.crispr_delta_dep, candidate.crispr_fdr,
                candidate.crispr_n_mutant, candidate.crispr_n_wt,
                candidate.prism_delta_auc, candidate.gdsc_delta_ic50,
                candidate.kb_clinical_hits, candidate.expression_corr,
                candidate.confidence_score, candidate.candidate_tier,
                candidate.evidence_summary, now,
                candidate.depmap_release, candidate.source_pipeline,
            ))
            con.commit()

            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE gene=? AND axis=? AND (cancer_type IS ? OR cancer_type=?)",
                (candidate.gene.upper(), candidate.axis,
                 candidate.cancer_type, candidate.cancer_type)
            ).fetchone()
            return _row_to_candidate(row)

    @staticmethod
    def list_pending(
        min_confidence: float = 0.0,
        gene: Optional[str] = None,
        axis: Optional[str] = None,
        limit: int = 100,
    ) -> List[ReceiptCandidate]:
        """List pending candidates sorted by confidence descending."""
        clauses = ["audit_status = 'pending'", "confidence_score >= ?"]
        params: list = [min_confidence]
        if gene:
            clauses.append("gene = ?")
            params.append(gene.upper())
        if axis:
            clauses.append("axis = ?")
            params.append(axis)
        params.append(limit)
        with _conn() as con:
            rows = con.execute(
                f"SELECT * FROM receipt_candidates WHERE {' AND '.join(clauses)} "
                f"ORDER BY confidence_score DESC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_candidate(r) for r in rows]

    @staticmethod
    def get(candidate_id: int) -> Optional[ReceiptCandidate]:
        with _conn() as con:
            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
        return _row_to_candidate(row) if row else None

    @staticmethod
    def approve(candidate_id: int, notes: Optional[str], audited_by: str) -> Optional[ReceiptCandidate]:
        """Mark a candidate as approved. Does NOT auto-promote to frozen receipts."""
        now = datetime.utcnow().isoformat()
        with _conn() as con:
            con.execute(
                "UPDATE receipt_candidates SET audit_status='approved', "
                "audit_notes=?, audited_by=?, audited_at=? WHERE id=?",
                (notes, audited_by, now, candidate_id)
            )
            con.commit()
            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
        return _row_to_candidate(row) if row else None

    @staticmethod
    def reject(candidate_id: int, notes: Optional[str], audited_by: str) -> Optional[ReceiptCandidate]:
        """Mark a candidate as rejected."""
        now = datetime.utcnow().isoformat()
        with _conn() as con:
            con.execute(
                "UPDATE receipt_candidates SET audit_status='rejected', "
                "audit_notes=?, audited_by=?, audited_at=? WHERE id=?",
                (notes, audited_by, now, candidate_id)
            )
            con.commit()
            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
        return _row_to_candidate(row) if row else None

    @staticmethod
    def stats() -> dict:
        with _conn() as con:
            total_pending = con.execute(
                "SELECT COUNT(*) FROM receipt_candidates WHERE audit_status='pending'"
            ).fetchone()[0]
            total_approved = con.execute(
                "SELECT COUNT(*) FROM receipt_candidates WHERE audit_status='approved'"
            ).fetchone()[0]
            total_rejected = con.execute(
                "SELECT COUNT(*) FROM receipt_candidates WHERE audit_status='rejected'"
            ).fetchone()[0]
            total_promoted = con.execute(
                "SELECT COUNT(*) FROM receipt_candidates WHERE promoted_to_frozen=1"
            ).fetchone()[0]
            high_conf = con.execute(
                "SELECT COUNT(*) FROM receipt_candidates WHERE audit_status='pending' AND confidence_score>=0.70"
            ).fetchone()[0]
        return {
            "total_pending": total_pending,
            "total_approved": total_approved,
            "total_rejected": total_rejected,
            "total_promoted": total_promoted,
            "high_confidence_pending": high_conf,
        }

    @staticmethod
    def coverage() -> List[dict]:
        """
        Return all (gene, axis) pairs tracked in the queue.
        Merged with frozen receipt counts from literature_receipts.
        Called by GET /audit/coverage.
        """
        from sl_agent.multimodal.literature_receipts import _FROZEN_RECEIPTS
        from sl_agent.multimodal.models import CandidateAxis

        # Axis label map
        axis_labels = {
            "cytidine_analogs": "Cytidine Analogs (gemcitabine, cytarabine)",
            "parp_inhibitors":  "PARP Inhibitors (olaparib, niraparib, talazoparib, rucaparib)",
            "atr_wee1":         "ATR / WEE1 Inhibitors (berzosertib, adavosertib, ceralasertib)",
            "wrn":              "WRN Helicase Inhibitors (VX-803, MRTX1719)",
            "immunotherapy":    "Immunotherapy / Checkpoint Inhibitors (PD-1/PD-L1)",
            "pkmyt1":           "PKMYT1 Kinase Inhibitors (RP-6306 class)",
        }

        # Build frozen receipt index: (gene, axis) → count + clinical flag
        frozen_index: dict = {}
        for (gene, axis_enum), receipts in _FROZEN_RECEIPTS.items():
            key = (gene.upper(), axis_enum.value)
            has_clinical = any(
                r.status.value == "positive"
                for k, r in receipts.items()
                if k == "clinical"
            )
            frozen_index[key] = {
                "count": frozen_index.get(key, {}).get("count", 0) + len(receipts),
                "has_clinical": has_clinical,
            }

        # Pull candidate aggregates from DB
        with _conn() as con:
            rows = con.execute("""
                SELECT gene, axis,
                       SUM(CASE WHEN audit_status='pending' THEN 1 ELSE 0 END) AS pending_count,
                       SUM(CASE WHEN audit_status='approved' THEN 1 ELSE 0 END) AS approved_count,
                       MAX(CASE WHEN audit_status='pending' THEN confidence_score ELSE NULL END) AS best_confidence,
                       MAX(candidate_tier) AS candidate_tier
                FROM receipt_candidates
                GROUP BY gene, axis
            """).fetchall()

        # Merge: frozen pairs + candidate pairs + all known axes for frozen genes
        all_pairs: dict = {}

        # Seed from frozen receipts
        for (gene, axis_val), info in frozen_index.items():
            key = (gene, axis_val)
            all_pairs[key] = {
                "gene": gene,
                "axis": axis_val,
                "axis_label": axis_labels.get(axis_val, axis_val),
                "tier": "Validated SL therapeutic lever" if info["has_clinical"] else "Strong candidate dependency axis",
                "frozen_receipt_count": info["count"],
                "candidate_count": 0,
                "approved_count": 0,
                "confidence": None,
                "has_clinical_evidence": info["has_clinical"],
            }

        # Overlay candidate data
        for row in rows:
            key = (row["gene"].upper(), row["axis"])
            if key not in all_pairs:
                all_pairs[key] = {
                    "gene": row["gene"].upper(),
                    "axis": row["axis"],
                    "axis_label": axis_labels.get(row["axis"], row["axis"]),
                    "tier": row["candidate_tier"] or "Mechanistic candidate only",
                    "frozen_receipt_count": 0,
                    "candidate_count": 0,
                    "approved_count": 0,
                    "confidence": None,
                    "has_clinical_evidence": False,
                }
            all_pairs[key]["candidate_count"] = row["pending_count"] or 0
            all_pairs[key]["approved_count"] = row["approved_count"] or 0
            all_pairs[key]["confidence"] = row["best_confidence"]

        return sorted(all_pairs.values(), key=lambda x: (x["gene"], x["axis"]))
```

## NEW: sl_agent/audit/routes.py

```python
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
  GET  /audit/{id}               — single candidate detail
  POST /audit/{id}/approve       — mark approved (does NOT auto-promote)
  POST /audit/{id}/reject        — mark rejected
  GET  /audit/stats              — queue summary counts
  GET  /audit/coverage           — gene × axis coverage grid (for heat map)
  POST /audit/seed               — inject a test candidate (dev/test only)
"""
from __future__ import annotations

import logging
from typing import List, Optional

from fastapi import APIRouter, HTTPException, Query, status

from .models import (
    AuditAction,
    CoverageEntry,
    QueueStats,
    ReceiptCandidate,
    RUO_DISCLAIMER,
    RUOEnvelope,
)
from .queue import AuditQueue, init_db

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


# ── GET /audit/{id} ───────────────────────────────────────────────────────────

@audit_router.get(
    "/{candidate_id}",
    summary="Get single candidate detail",
)
async def get_candidate(candidate_id: int):
    c = AuditQueue.get(candidate_id)
    if not c:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    return _wrap(c.model_dump())


# ── POST /audit/{id}/approve ──────────────────────────────────────────────────

@audit_router.post(
    "/{candidate_id}/approve",
    summary="Approve a receipt candidate",
    description=(
        "Marks candidate as approved. Does NOT auto-promote to frozen receipts — "
        "human must explicitly run the promotion workflow. "
        "Once approved, the candidate will not be overwritten by future mining runs."
    ),
)
async def approve_candidate(candidate_id: int, action: AuditAction):
    updated = AuditQueue.approve(candidate_id, action.notes, action.audited_by)
    if not updated:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    return _wrap(updated.model_dump())


# ── POST /audit/{id}/reject ───────────────────────────────────────────────────

@audit_router.post(
    "/{candidate_id}/reject",
    summary="Reject a receipt candidate",
    description="Marks candidate as rejected with auditor notes.",
)
async def reject_candidate(candidate_id: int, action: AuditAction):
    updated = AuditQueue.reject(candidate_id, action.notes, action.audited_by)
    if not updated:
        raise HTTPException(status_code=404, detail=f"Candidate {candidate_id} not found")
    return _wrap(updated.model_dump())


# ── GET /audit/stats ──────────────────────────────────────────────────────────

@audit_router.get(
    "/stats",
    summary="Audit queue summary statistics",
)
async def queue_stats():
    s = AuditQueue.stats()
    return _wrap(s)


# ── GET /audit/coverage ───────────────────────────────────────────────────────

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
async def coverage_grid():
    entries = AuditQueue.coverage()
    return _wrap(entries)


# ── POST /audit/seed (dev/test only) ─────────────────────────────────────────

@audit_router.post(
    "/seed",
    summary="Seed a test candidate (dev/test only)",
    description="Injects a ReceiptCandidate directly into the queue. Dev/test use only.",
    include_in_schema=True,  # set False in production
)
async def seed_candidate(candidate: ReceiptCandidate):
    saved = AuditQueue.upsert(candidate)
    return _wrap(saved.model_dump())
```

## MODIFY: sl_agent/api/app.py

Add import:
```python
from ..audit.routes import audit_router
```

Add to `create_app()` after existing `app.include_router(mm_router, ...)`:
```python
    app.include_router(audit_router, prefix="/api/v1")
```

Also add to lifespan startup:
```python
    from ..audit.queue import init_db as init_audit_db
    try:
        init_audit_db()
        log.info("startup", message="Audit queue DB initialized")
    except Exception as exc:
        log.warning("startup_warning", error=str(exc), message="Audit DB init failed")
```

## MODIFY: sl_agent/core/config.py

Add ONE field to the Settings class (after `depmap_cache_dir`):
```python
    audit_db_path: Path = Path(".cache/audit_queue.db")
```

## NEW: sl_agent/audit/__init__.py

```python
"""Audit queue — receipt candidate lifecycle and doctor portal API."""
```

## Tests: sl_agent/audit/tests/test_audit.py

```
sl_agent/audit/tests/
├── __init__.py
└── test_audit.py
```

Tests to write (25 total):

### TestReceiptCandidateModel (5 tests)
- valid candidate with all fields
- confidence_score out of range raises validation error
- audit_status defaults to "pending"
- promoted_to_frozen defaults to False
- generated_at auto-populated

### TestAuditQueue (10 tests)
- upsert new candidate returns with id
- upsert duplicate (pending) updates confidence
- upsert duplicate (approved) does NOT overwrite
- upsert duplicate (rejected) does NOT overwrite
- list_pending returns sorted by confidence desc
- list_pending min_confidence filter works
- list_pending gene filter works
- approve sets status and audited_by
- reject sets status and audited_by
- stats returns correct counts after operations

### TestAuditRoutes (10 tests — use FastAPI TestClient)
- GET /api/v1/audit/queue returns RUO envelope
- GET /api/v1/audit/queue with gene filter
- GET /api/v1/audit/{id} returns candidate
- GET /api/v1/audit/{id} missing returns 404
- POST /api/v1/audit/{id}/approve updates status
- POST /api/v1/audit/{id}/reject updates status
- GET /api/v1/audit/stats returns correct keys
- GET /api/v1/audit/coverage returns list
- GET /api/v1/audit/coverage includes frozen receipt entries
- POST /api/v1/audit/seed injects candidate

IMPORTANT for tests: Use a temporary SQLite DB (tmp_path fixture or monkeypatch settings).
Do NOT use the real .cache/audit_queue.db during tests.

## Run Tests

```bash
cd /home/user/workspace && python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ sl_agent/audit/tests/ -v
```

All 100 existing tests must pass + 25 new audit tests = 125 total.

## INTEGRATION_SPRINT1.md — Wire-In Instructions

Write this file at `sl_patch_sprint1/INTEGRATION_SPRINT1.md`:

```markdown
# Sprint 1 — Audit Queue + Doctor Portal API Layer

## Files Added / Modified

| File | Action |
|------|--------|
| sl_agent/audit/__init__.py | NEW |
| sl_agent/audit/models.py | NEW |
| sl_agent/audit/queue.py | NEW |
| sl_agent/audit/routes.py | NEW |
| sl_agent/audit/tests/__init__.py | NEW |
| sl_agent/audit/tests/test_audit.py | NEW |
| sl_agent/api/app.py | MODIFIED — +audit_router |
| sl_agent/core/config.py | MODIFIED — +audit_db_path |

## Install

1. Drop the `sl_agent/audit/` directory into your project
2. Apply the two-line patch to `sl_agent/api/app.py` (see below)
3. Apply the one-line patch to `sl_agent/core/config.py` (see below)
4. Run: `python -m pytest sl_agent/audit/tests/ -v`

### app.py patch

In `create_app()`, after the last `app.include_router(...)` call:
```python
from ..audit.routes import audit_router
app.include_router(audit_router, prefix="/api/v1")
```

In `lifespan()`, after the DataStore.ensure_loaded block:
```python
from ..audit.queue import init_db as init_audit_db
init_audit_db()
```

### config.py patch

In the Settings class, add:
```python
audit_db_path: Path = Path(".cache/audit_queue.db")
```

## New Endpoints

| Method | Path | Purpose | Doctor Portal View |
|--------|------|---------|-------------------|
| GET | /api/v1/audit/queue | Pending candidates, sorted by confidence | View 2 (Audit Queue) |
| GET | /api/v1/audit/queue?gene=MBD4 | Patient-specific candidates | View 1 (Case View) |
| GET | /api/v1/audit/{id} | Single candidate detail | View 2 detail |
| POST | /api/v1/audit/{id}/approve | Human approval | View 2 action |
| POST | /api/v1/audit/{id}/reject | Human rejection | View 2 action |
| GET | /api/v1/audit/stats | Queue summary | Dashboard widget |
| GET | /api/v1/audit/coverage | Gene × axis grid | View 3 (Heat Map) |
| POST | /api/v1/audit/seed | Inject test candidate | Dev/test only |

## RUO Envelope

Every response is wrapped:
```json
{
  "ruo": true,
  "disclaimer": "Research use only. Not validated for clinical decision-making. All receipt candidates require human expert review before promotion.",
  "data": { ... }
}
```

## Doctor Portal — Future Views (Sprint 2+)

These views consume this API. No backend changes required:

**View 1 — Case View** (per-patient ranked interventions)
- Call: GET /api/v1/audit/queue?gene={gene}&min_confidence=0.5
- Display: ranked table of candidates, approve/flag per axis
- Data already in ReceiptCandidate: gene, axis, candidate_tier, confidence_score, evidence_summary

**View 2 — Audit Queue** (global pending candidates)
- Call: GET /api/v1/audit/queue (no filter)
- Actions: POST /audit/{id}/approve or /reject with AuditAction body
- Batch pattern: loop over selected IDs, POST each

**View 3 — Coverage Heat Map** (Atlas Panel 1)
- Call: GET /api/v1/audit/coverage
- Render as gene × axis grid
- Color by tier: Validated (green) → Strong (blue) → Mechanistic (yellow) → no data (grey)
- Overlay: candidate_count badge when pending > 0
- CoverageEntry fields map directly to heat map cell attributes

## What Sprint 2 Adds

- `receipt_miner.py` — auto-populates the queue from DepMap CRISPR + GDSC ANOVA + KB
- `POST /audit/mine` — triggers mining run for a gene panel
- `gdsc_biomarker_loader.py` — pre-computed GDSC ANOVA gene-drug associations
- The queue you built in Sprint 1 is the destination — no Sprint 1 changes needed
```

## Package

```bash
cd /home/user/workspace
rm -rf sl_patch_sprint1/
mkdir -p sl_patch_sprint1/sl_agent/audit/tests
mkdir -p sl_patch_sprint1/sl_agent/api
mkdir -p sl_patch_sprint1/sl_agent/core

cp sl_agent/audit/__init__.py          sl_patch_sprint1/sl_agent/audit/
cp sl_agent/audit/models.py            sl_patch_sprint1/sl_agent/audit/
cp sl_agent/audit/queue.py             sl_patch_sprint1/sl_agent/audit/
cp sl_agent/audit/routes.py            sl_patch_sprint1/sl_agent/audit/
cp sl_agent/audit/tests/__init__.py    sl_patch_sprint1/sl_agent/audit/tests/
cp sl_agent/audit/tests/test_audit.py  sl_patch_sprint1/sl_agent/audit/tests/
cp sl_agent/api/app.py                 sl_patch_sprint1/sl_agent/api/
cp sl_agent/core/config.py             sl_patch_sprint1/sl_agent/core/

# Write INTEGRATION_SPRINT1.md (see above)

zip -r sl_agent_sprint1_audit.zip sl_patch_sprint1/
rm -rf sl_patch_sprint1/
ls -lh sl_agent_sprint1_audit.zip
```
