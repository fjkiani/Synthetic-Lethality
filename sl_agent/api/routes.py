"""
FastAPI route definitions for the SL Mapping Agent.

Endpoints:
  POST /analyze          — primary SL query
  GET  /genes            — list all genes available in current CRISPR matrix
  GET  /cancer_types     — list available cancer lineages
  GET  /health           — readiness check
  GET  /result/{job_id}  — async job polling (for heavy queries)
"""
from __future__ import annotations

import asyncio
import logging
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, Optional

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, status
from fastapi.responses import JSONResponse

from ..core.config import get_settings
from ..core.models import (
    ErrorResponse,
    SLMapResult,
    SLQueryInput,
    SLQueryResponse,
)
from ..core.orchestrator import DataStore, run_sl_analysis

logger = logging.getLogger(__name__)
cfg = get_settings()

router = APIRouter()

# Thread pool for CPU-bound analysis
_executor = ThreadPoolExecutor(max_workers=4)

# In-memory job store (swap for Redis/Celery in production)
_jobs: Dict[str, Any] = {}


# ── Dependency ────────────────────────────────────────────────────────────────

async def verify_data_loaded():
    """Ensure DepMap data is loaded before accepting queries."""
    try:
        DataStore.ensure_loaded(require_prism=True)
    except Exception as exc:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=f"Data not yet loaded: {exc}",
        )


# ── Health ────────────────────────────────────────────────────────────────────

@router.get("/health", tags=["System"])
async def health_check():
    """Readiness probe."""
    data_ready = DataStore._crispr is not None
    return {
        "status": "ok" if data_ready else "loading",
        "data_loaded": data_ready,
        "version": cfg.app_version,
        "depmap_release": cfg.depmap_release,
    }


# ── Gene / lineage lookups ────────────────────────────────────────────────────

@router.get("/genes", tags=["Metadata"], dependencies=[Depends(verify_data_loaded)])
async def list_genes(prefix: Optional[str] = None, limit: int = 100):
    """List genes available in the CRISPR dependency matrix."""
    genes = DataStore.crispr().columns.tolist()
    if prefix:
        genes = [g for g in genes if g.upper().startswith(prefix.upper())]
    return {"count": len(genes), "genes": genes[:limit]}


@router.get("/cancer_types", tags=["Metadata"], dependencies=[Depends(verify_data_loaded)])
async def list_cancer_types():
    """List available cancer lineages / primary disease labels."""
    si = DataStore.sample_info()
    lineage_col = next(
        (c for c in si.columns if c in ("OncotreeLineage", "lineage")), None
    )
    pd_col = next(
        (c for c in si.columns if c in ("OncotreePrimaryDisease", "primary_disease")), None
    )
    result: Dict[str, Any] = {}
    if lineage_col:
        result["lineages"] = sorted(si[lineage_col].dropna().unique().tolist())
    if pd_col:
        result["primary_diseases"] = sorted(si[pd_col].dropna().unique().tolist())
    return result


# ── Synchronous analysis endpoint ─────────────────────────────────────────────

@router.post(
    "/analyze",
    response_model=SLQueryResponse,
    tags=["Analysis"],
    summary="Run full SL analysis for a query gene",
)
async def analyze(
    query: SLQueryInput,
    background_tasks: BackgroundTasks,
    _: None = Depends(verify_data_loaded),
):
    """
    Execute the complete 4-step SL pipeline:
      1. Candidate SL partner identification (stratified Wilcoxon + FDR)
      2. Gene → drug mapping (PRISM + ChEMBL + OncoKB)
      3. Ranking (SL signal × drug differential × druggability)
      4. Structured JSON report

    For large analyses this may take 30–120 s depending on number of partners.
    Use the async `/analyze/async` endpoint for non-blocking operation.
    """
    t0 = time.perf_counter()
    loop = asyncio.get_event_loop()
    try:
        result: SLMapResult = await loop.run_in_executor(
            _executor, run_sl_analysis, query
        )
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc))
    except RuntimeError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(exc)
        )
    except Exception as exc:
        logger.exception("Unhandled error during SL analysis")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Internal error: {exc}",
        )

    return SLQueryResponse(
        status="success",
        result=result,
        processing_time_ms=round((time.perf_counter() - t0) * 1000, 1),
    )


# ── Async job endpoint ────────────────────────────────────────────────────────

@router.post(
    "/analyze/async",
    tags=["Analysis"],
    summary="Submit analysis job (non-blocking)",
    status_code=202,
)
async def analyze_async(
    query: SLQueryInput,
    _: None = Depends(verify_data_loaded),
):
    """
    Submit an SL analysis job. Returns a job_id immediately.
    Poll `GET /result/{job_id}` for results.
    """
    job_id = str(uuid.uuid4())
    _jobs[job_id] = {"status": "queued", "result": None, "error": None}

    async def _run():
        _jobs[job_id]["status"] = "running"
        loop = asyncio.get_event_loop()
        try:
            result = await loop.run_in_executor(_executor, run_sl_analysis, query)
            _jobs[job_id]["status"] = "done"
            _jobs[job_id]["result"] = result.model_dump()
        except Exception as exc:
            _jobs[job_id]["status"] = "error"
            _jobs[job_id]["error"] = str(exc)

    asyncio.create_task(_run())
    return {"job_id": job_id, "status": "queued"}


@router.get(
    "/result/{job_id}",
    tags=["Analysis"],
    summary="Poll async job result",
)
async def get_result(job_id: str):
    if job_id not in _jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    job = _jobs[job_id]
    if job["status"] == "done":
        return {"status": "done", "result": job["result"]}
    return {"status": job["status"], "error": job.get("error")}
