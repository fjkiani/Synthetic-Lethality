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
