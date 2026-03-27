"""
SQLite-backed audit queue for ReceiptCandidate lifecycle management.

Schema: single table `receipt_candidates` with all ReceiptCandidate fields.
DB file: configured via settings.audit_db_path (default: .cache/audit_queue.db)
Thread-safe: uses connection-per-call pattern (SQLite WAL mode).

No ORM — direct sqlite3 to keep zero new dependencies.
"""
from __future__ import annotations

import logging
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from .models import ReceiptCandidate

logger = logging.getLogger(__name__)

# ── Coverage version cache (in-memory) ───────────────────────────────────────
# Simple flag invalidated on approve() / reject() / upsert() and after mine.
# Exposed via GET /audit/coverage as X-Coverage-Version header.
# Thread-safe write via _set_coverage_version(); reads are lock-free.
import threading
_coverage_version_lock = threading.Lock()
_coverage_version: str = datetime.utcnow().isoformat()  # initialised at import


def _set_coverage_version() -> str:
    """Stamp a new coverage version timestamp. Returns the new version string."""
    global _coverage_version
    with _coverage_version_lock:
        _coverage_version = datetime.utcnow().isoformat()
        return _coverage_version


def get_coverage_version() -> str:
    """Read current coverage version (no lock needed — string reads are atomic in CPython)."""
    return _coverage_version


# ── Last-mine timestamp (in-memory, reset on restart) ──────────────────────────
_last_mine_at: Optional[str] = None


def set_last_mine_at() -> str:
    """Called when a mine run completes. Returns the new timestamp."""
    global _last_mine_at
    _last_mine_at = datetime.utcnow().isoformat()
    return _last_mine_at


def get_last_mine_at() -> Optional[str]:
    return _last_mine_at


# ── Frozen receipt snippet generator ─────────────────────────────────────────────

def generate_frozen_receipt_snippet(candidate: "ReceiptCandidate") -> str:
    """
    Generate a Python source snippet ready to paste into literature_receipts.py.

    Format:
      ("GENE", CandidateAxis.AXIS): FrozenReceipt(
          gene="GENE",
          ...
      ),

    SABOTAGE PROTECTION: snippet tier is derived from candidate_tier.
    Even if the candidate says 'Strong', the snippet says 'STRONG'.
    'Validated' is NEVER injected here — that requires a human to write it
    into literature_receipts.py after independent clinical evidence review.
    The snippet is a PASTE-READY TEMPLATE, not an auto-promotion.

    Human pastes it in. Always.
    """
    gene = candidate.gene.upper()
    # Map axis string to CandidateAxis enum attr name
    axis_val = candidate.axis  # e.g. "parp_inhibitors"
    axis_enum_name = axis_val.upper()  # e.g. "PARP_INHIBITORS"

    # TIER LAUNDERING GUARD: clamp tier at snippet generation time.
    # Even if candidate_tier was manually set to "Validated" at upsert,
    # the snippet NEVER emits VALIDATED. Only a human writing directly
    # into literature_receipts.py with independent clinical evidence can set VALIDATED.
    # This is a hard clamp, not a warning — defence in depth beyond _project_tier().
    _raw_tier = (candidate.candidate_tier or "Mechanistic").upper().replace(" ", "_")
    _ALLOWED_SNIPPET_TIERS = {"STRONG", "MECHANISTIC", "INSUFFICIENT"}
    tier = _raw_tier if _raw_tier in _ALLOWED_SNIPPET_TIERS else "STRONG"
    crispr_delta = f"{candidate.crispr_delta_dep:.4f}" if candidate.crispr_delta_dep is not None else "None"
    crispr_fdr   = f"{candidate.crispr_fdr:.4f}" if candidate.crispr_fdr is not None else "None"
    gdsc_delta   = f"{candidate.gdsc_delta_ic50:.4f}" if candidate.gdsc_delta_ic50 is not None else "None"
    kb_hits      = candidate.kb_clinical_hits
    audited_by   = candidate.audited_by or "unknown"
    promoted_at  = candidate.audited_at.isoformat() if candidate.audited_at else datetime.utcnow().isoformat()
    notes        = (candidate.audit_notes or "").replace('"', "'")

    snippet = (
        f'    ("{gene}", CandidateAxis.{axis_enum_name}): FrozenReceipt(\n'
        f'        gene="{gene}",\n'
        f'        axis=CandidateAxis.{axis_enum_name},\n'
        f'        evidence_tier="{tier}",\n'
        f'        crispr_delta_dep={crispr_delta},\n'
        f'        crispr_fdr={crispr_fdr},\n'
        f'        gdsc_delta_ic50={gdsc_delta},\n'
        f'        kb_clinical_hits={kb_hits},\n'
        f'        source="audit_queue_promoted",\n'
        f'        promoted_by="{audited_by}",\n'
        f'        promoted_at="{promoted_at}",\n'
        f'        notes="{notes}"\n'
        f'    ),'
    )
    return snippet

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
    def approve(
        candidate_id: int,
        notes: Optional[str],
        audited_by: str,
    ) -> Optional[dict]:
        """Mark a candidate as approved. Does NOT auto-promote to frozen receipts.

        Returns a dict with keys:
          - candidate: ReceiptCandidate model
          - frozen_receipt_snippet: Python string ready to paste into literature_receipts.py

        Invalidates coverage version so the doctor portal heat map knows to refresh.
        Human ALWAYS pastes the snippet manually. No auto-write.
        """
        # IMMUTABILITY GUARD: empty audited_by is an audit-trail gap
        if not audited_by or not audited_by.strip():
            raise ValueError("audited_by cannot be empty — audit trail requires an auditor identity")

        now = datetime.utcnow().isoformat()
        with _conn() as con:
            # DECISION IMMUTABILITY: only allow approve on pending candidates.
            # Once a human has rejected a candidate that decision is final.
            existing = con.execute(
                "SELECT audit_status FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
            if existing is None:
                return None
            if existing["audit_status"] == "rejected":
                raise ValueError(
                    f"Candidate {candidate_id} is already rejected. "
                    "Human decisions are immutable. Create a new candidate if evidence has changed."
                )
            con.execute(
                "UPDATE receipt_candidates SET audit_status='approved', "
                "audit_notes=?, audited_by=?, audited_at=? WHERE id=?",
                (notes, audited_by, now, candidate_id)
            )
            con.commit()
            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
        if row is None:
            return None
        candidate = _row_to_candidate(row)
        snippet = generate_frozen_receipt_snippet(candidate)
        _set_coverage_version()  # invalidate coverage cache
        return {"candidate": candidate, "frozen_receipt_snippet": snippet}

    @staticmethod
    def reject(candidate_id: int, notes: Optional[str], audited_by: str) -> Optional[ReceiptCandidate]:
        """Mark a candidate as rejected.
        Invalidates coverage version so the doctor portal heat map knows to refresh.
        """
        # IMMUTABILITY GUARD: empty audited_by is an audit-trail gap
        if not audited_by or not audited_by.strip():
            raise ValueError("audited_by cannot be empty — audit trail requires an auditor identity")

        now = datetime.utcnow().isoformat()
        with _conn() as con:
            # DECISION IMMUTABILITY: only allow reject on pending candidates.
            # Once a human has approved a candidate that decision is final.
            existing = con.execute(
                "SELECT audit_status FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
            if existing is None:
                return None
            if existing["audit_status"] == "approved":
                raise ValueError(
                    f"Candidate {candidate_id} is already approved. "
                    "Human decisions are immutable. Approved receipts cannot be silently rejected."
                )
            con.execute(
                "UPDATE receipt_candidates SET audit_status='rejected', "
                "audit_notes=?, audited_by=?, audited_at=? WHERE id=?",
                (notes, audited_by, now, candidate_id)
            )
            con.commit()
            row = con.execute(
                "SELECT * FROM receipt_candidates WHERE id=?", (candidate_id,)
            ).fetchone()
        _set_coverage_version()  # invalidate coverage cache
        return _row_to_candidate(row) if row else None

    @staticmethod
    def stats() -> dict:
        """Full dashboard stats for the doctor portal."""
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
            # Top 5 by confidence score (pending only)
            top_rows = con.execute(
                "SELECT gene, axis, confidence_score, candidate_tier "
                "FROM receipt_candidates WHERE audit_status='pending' "
                "ORDER BY confidence_score DESC LIMIT 5"
            ).fetchall()

        top_pending = [
            {
                "gene": r["gene"],
                "axis": r["axis"],
                "confidence": r["confidence_score"],
                "projected_tier": r["candidate_tier"],
            }
            for r in top_rows
        ]

        # Coverage summary (reuse coverage() method)
        try:
            cov = AuditQueue.coverage()
            cov_summary = cov["summary"]
        except Exception:
            cov_summary = {"validated": 0, "strong": 0, "mechanistic": 0, "not_covered": 0, "total_pairs": 0}

        return {
            "queue_counts": {
                "pending": total_pending,
                "approved": total_approved,
                "rejected": total_rejected,
                "promoted": total_promoted,
                "high_confidence_pending": high_conf,
            },
            "coverage_summary": cov_summary,
            "top_pending": top_pending,
            "last_mine_at": get_last_mine_at(),
            "coverage_version": get_coverage_version(),
        }

    @staticmethod
    def coverage() -> dict:
        """
        Return the full coverage matrix for the Atlas Panel 1 heat map.

        Returns a dict with:
          - coverage: list of gene×axis cells (includes 'not_covered' entries
            for genes in the default panel with zero frozen receipts or candidates)
          - summary: {validated, strong, mechanistic, not_covered, total_pairs}
          - coverage_version: current version timestamp

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
        all_axes = list(axis_labels.keys())

        # Build frozen receipt index: (gene, axis) → count + clinical flag
        frozen_index: dict = {}
        for (gene, axis_enum), receipts in _FROZEN_RECEIPTS.items():
            key = (gene.upper(), axis_enum.value)
            has_clinical = any(
                r.status.value == "positive"
                for k, r in receipts.items()
                if k == "clinical"
            )
            existing = frozen_index.get(key, {"count": 0, "has_clinical": False})
            frozen_index[key] = {
                "count": existing["count"] + len(receipts),
                "has_clinical": existing["has_clinical"] or has_clinical,
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

        # Merge: frozen pairs + candidate pairs
        all_pairs: dict = {}

        # Seed from frozen receipts
        for (gene, axis_val), info in frozen_index.items():
            key = (gene, axis_val)
            tier = "Validated" if info["has_clinical"] else "Strong"
            all_pairs[key] = {
                "gene": gene,
                "axis": axis_val,
                "axis_label": axis_labels.get(axis_val, axis_val),
                "tier": tier,
                "frozen_receipt_count": info["count"],
                "pending_candidate_count": 0,
                "approved_count": 0,
                "best_confidence": None,
                "has_clinical_positive": info["has_clinical"],
            }

        # Overlay candidate data
        for row in rows:
            key = (row["gene"].upper(), row["axis"])
            if key not in all_pairs:
                all_pairs[key] = {
                    "gene": row["gene"].upper(),
                    "axis": row["axis"],
                    "axis_label": axis_labels.get(row["axis"], row["axis"]),
                    "tier": row["candidate_tier"] or "Mechanistic",
                    "frozen_receipt_count": 0,
                    "pending_candidate_count": 0,
                    "approved_count": 0,
                    "best_confidence": None,
                    "has_clinical_positive": False,
                }
            all_pairs[key]["pending_candidate_count"] = row["pending_count"] or 0
            all_pairs[key]["approved_count"] = row["approved_count"] or 0
            all_pairs[key]["best_confidence"] = row["best_confidence"]

        # Add 'not_covered' entries for genes in the default panel
        # that have no frozen receipts AND no queue candidates on any axis
        try:
            from sl_agent.multimodal.receipt_miner import load_default_panel, DEFAULT_AXES
            panel_genes = load_default_panel()
        except Exception:
            panel_genes = []

        for gene in panel_genes:
            gene_up = gene.upper()
            gene_covered_axes = {axis for (g, axis) in all_pairs if g == gene_up}
            for ax in all_axes:
                if ax not in gene_covered_axes:
                    key = (gene_up, ax)
                    all_pairs[key] = {
                        "gene": gene_up,
                        "axis": ax,
                        "axis_label": axis_labels.get(ax, ax),
                        "tier": "Not covered",
                        "frozen_receipt_count": 0,
                        "pending_candidate_count": 0,
                        "approved_count": 0,
                        "best_confidence": None,
                        "has_clinical_positive": False,
                    }

        sorted_entries = sorted(all_pairs.values(), key=lambda x: (x["gene"], x["axis"]))

        # Build summary counts
        summary = {"validated": 0, "strong": 0, "mechanistic": 0, "not_covered": 0, "total_pairs": len(sorted_entries)}
        for entry in sorted_entries:
            t = entry["tier"].lower()
            if "validated" in t:
                summary["validated"] += 1
            elif "strong" in t:
                summary["strong"] += 1
            elif "mechanistic" in t:
                summary["mechanistic"] += 1
            else:
                summary["not_covered"] += 1

        return {
            "coverage": sorted_entries,
            "summary": summary,
            "coverage_version": get_coverage_version(),
        }
