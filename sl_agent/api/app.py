"""
FastAPI application factory.
"""
from __future__ import annotations

import logging
import time
from contextlib import asynccontextmanager

import structlog
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from ..core.config import get_settings
from ..core.orchestrator import DataStore
from .kb_routes import kb_router
from .multimodal_routes import mm_router
from .routes import router
from ..audit.routes import audit_router

cfg = get_settings()

# ── Structured logging ────────────────────────────────────────────────────────
structlog.configure(
    processors=[
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.stdlib.add_log_level,
        structlog.dev.ConsoleRenderer() if cfg.debug else structlog.processors.JSONRenderer(),
    ],
    wrapper_class=structlog.stdlib.BoundLogger,
    logger_factory=structlog.stdlib.LoggerFactory(),
)
logging.basicConfig(level=cfg.log_level.upper())
log = structlog.get_logger()


# ── Lifespan (startup / shutdown) ─────────────────────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    log.info("startup", message="Loading DepMap matrices — this may take a minute…")
    try:
        DataStore.ensure_loaded(require_prism=True)
        log.info("startup", message="Data loaded successfully")
    except Exception as exc:
        log.warning("startup_warning", error=str(exc), message="Data pre-load failed — will retry on first request")

    from ..audit.queue import init_db as init_audit_db
    try:
        init_audit_db()
        log.info("startup", message="Audit queue DB initialized")
    except Exception as exc:
        log.warning("startup_warning", error=str(exc), message="Audit DB init failed")

    # ── GDSC cache warm-up (Sprint 3) ─────────────────────────────────────────
    # Download + cache GDSC1 and GDSC2 parquet files at startup.
    # Timeout: 60 s — if it takes longer the server continues without blocking.
    # Never crashes startup on failure.
    import asyncio
    import concurrent.futures

    def _gdsc_warmup() -> str:
        """Blocking I/O — runs in thread pool so it doesn't block the event loop."""
        from ..data.gdsc_biomarker_loader import load_gdsc1, load_gdsc2
        total = 0
        try:
            df1 = load_gdsc1()
            total += len(df1)
        except Exception as e:
            log.warning("gdsc_warmup", dataset="GDSC1", error=str(e))
        try:
            df2 = load_gdsc2()
            total += len(df2)
        except Exception as e:
            log.warning("gdsc_warmup", dataset="GDSC2", error=str(e))
        return f"GDSC cache warm-up complete: {total} rows"

    loop = asyncio.get_event_loop()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            warmup_msg = await asyncio.wait_for(
                loop.run_in_executor(pool, _gdsc_warmup),
                timeout=60.0,
            )
        log.info("startup", message=warmup_msg)
    except asyncio.TimeoutError:
        log.warning("startup_warning", message="GDSC warm-up timed out after 60 s — skipping; will load on first mine run")
    except Exception as exc:
        log.warning("startup_warning", error=str(exc), message="GDSC warm-up failed — will load on demand")

    yield
    log.info("shutdown", message="SL Agent shutting down")


# ── App factory ───────────────────────────────────────────────────────────────

def create_app() -> FastAPI:
    app = FastAPI(
        title=cfg.app_name,
        version=cfg.app_version,
        description=(
            "Synthetic Lethality Mapping Agent v4. "
            "Integrates DepMap CRISPR screens, PRISM drug data, ChEMBL, and the open KB stack "
            "(CIViC + CGI + JAX) with a multi-modal evidence matrix engine that refuses to "
            "hallucinate SL where receipts are weak, and refuses to discard promising axes "
            "just because CRISPR is negative. Part of the CrisPRO precision oncology platform."
        ),
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=lifespan,
    )

    # CORS — restrict in production
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"] if cfg.debug else [],
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )

    # Request timing middleware
    @app.middleware("http")
    async def add_timing(request: Request, call_next):
        t0 = time.perf_counter()
        response = await call_next(request)
        ms = round((time.perf_counter() - t0) * 1000, 1)
        response.headers["X-Process-Time-ms"] = str(ms)
        return response

    # Global exception handler
    @app.exception_handler(Exception)
    async def global_error_handler(request: Request, exc: Exception):
        log.error("unhandled_exception", path=request.url.path, error=str(exc))
        return JSONResponse(
            status_code=500,
            content={"status": "error", "error": "Internal server error", "detail": str(exc)},
        )

    # Register routes
    app.include_router(router, prefix="/api/v1")
    app.include_router(kb_router, prefix="/api/v1")
    app.include_router(mm_router, prefix="/api/v1")
    app.include_router(audit_router, prefix="/api/v1")

    return app


app = create_app()
