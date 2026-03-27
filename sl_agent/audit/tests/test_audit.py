"""
25 tests for the audit queue package.

Test groups:
  TestReceiptCandidateModel (5) — Pydantic model validation
  TestAuditQueue            (10) — SQLite CRUD via AuditQueue
  TestAuditRoutes           (10) — FastAPI endpoints via TestClient

IMPORTANT: All tests use a temporary SQLite DB (tmp_path fixture).
           The real .cache/audit_queue.db is NEVER touched.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Generator

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError


# ── Helpers ────────────────────────────────────────────────────────────────────

def _make_candidate(**kwargs):
    """Factory: minimal valid ReceiptCandidate kwargs."""
    defaults = {
        "gene": "BRCA2",
        "axis": "parp_inhibitors",
        "confidence_score": 0.75,
        "candidate_tier": "Strong candidate dependency axis",
        "evidence_summary": "BRCA2 loss → PARP inhibitor sensitivity (CRISPR delta-dep=-0.4)",
    }
    defaults.update(kwargs)
    return defaults


@pytest.fixture()
def tmp_db(tmp_path, monkeypatch) -> Path:
    """
    Redirect _get_db_path() to a fresh temp DB for every test.
    Also clears the lru_cache on get_settings so audit_db_path is re-read.
    """
    db_file = tmp_path / "test_audit.db"

    # Patch the function that AuditQueue calls to locate the DB
    import sl_agent.audit.queue as queue_module
    monkeypatch.setattr(queue_module, "_get_db_path", lambda: db_file)

    # Initialise the temp DB schema
    queue_module.init_db()
    return db_file


@pytest.fixture()
def client(tmp_db, monkeypatch) -> Generator:
    """
    FastAPI TestClient with audit_router wired in and DB redirected to tmp_db.
    We build a minimal app that only includes the audit router so tests are
    self-contained and do not trigger DepMap data loading.
    """
    from fastapi import FastAPI
    from sl_agent.audit.routes import audit_router
    import sl_agent.audit.queue as queue_module

    # tmp_db fixture already patched _get_db_path, so routes will use it too.
    app = FastAPI()
    app.include_router(audit_router, prefix="/api/v1")

    with TestClient(app) as tc:
        yield tc


# ─────────────────────────────────────────────────────────────────────────────
# TestReceiptCandidateModel (5 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestReceiptCandidateModel:

    def test_valid_candidate_all_fields(self):
        """Valid candidate with all optional fields populated."""
        from sl_agent.audit.models import ReceiptCandidate
        c = ReceiptCandidate(
            gene="MBD4",
            axis="cytidine_analogs",
            cancer_type="colorectal",
            crispr_delta_dep=-0.35,
            crispr_fdr=0.05,
            crispr_n_mutant=15,
            crispr_n_wt=120,
            prism_delta_auc=-0.12,
            gdsc_delta_ic50=-1.5,
            kb_clinical_hits=3,
            expression_corr=-0.42,
            confidence_score=0.88,
            candidate_tier="Validated SL therapeutic lever",
            evidence_summary="MBD4 loss confers cytidine analog sensitivity",
            depmap_release="24Q4",
        )
        assert c.gene == "MBD4"
        assert c.confidence_score == 0.88
        assert c.audit_status == "pending"
        assert c.promoted_to_frozen is False

    def test_confidence_score_out_of_range_raises(self):
        """confidence_score > 1.0 must raise ValidationError."""
        from sl_agent.audit.models import ReceiptCandidate
        with pytest.raises(ValidationError):
            ReceiptCandidate(**_make_candidate(confidence_score=1.5))

    def test_confidence_score_negative_raises(self):
        """confidence_score < 0.0 must raise ValidationError."""
        from sl_agent.audit.models import ReceiptCandidate
        with pytest.raises(ValidationError):
            ReceiptCandidate(**_make_candidate(confidence_score=-0.1))

    def test_audit_status_defaults_to_pending(self):
        """audit_status must default to 'pending'."""
        from sl_agent.audit.models import ReceiptCandidate
        c = ReceiptCandidate(**_make_candidate())
        assert c.audit_status == "pending"

    def test_promoted_to_frozen_defaults_false(self):
        """promoted_to_frozen must default to False."""
        from sl_agent.audit.models import ReceiptCandidate
        c = ReceiptCandidate(**_make_candidate())
        assert c.promoted_to_frozen is False

    def test_generated_at_auto_populated(self):
        """generated_at must be auto-populated (not None)."""
        from sl_agent.audit.models import ReceiptCandidate
        from datetime import datetime
        c = ReceiptCandidate(**_make_candidate())
        assert isinstance(c.generated_at, datetime)


# ─────────────────────────────────────────────────────────────────────────────
# TestAuditQueue (10 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestAuditQueue:

    def test_upsert_new_candidate_returns_with_id(self, tmp_db):
        """Upserting a new candidate assigns a DB id."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c = ReceiptCandidate(**_make_candidate(gene="TP53", axis="parp_inhibitors"))
        saved = AuditQueue.upsert(c)
        assert saved.id is not None
        assert saved.id > 0
        assert saved.gene == "TP53"

    def test_upsert_duplicate_pending_updates_confidence(self, tmp_db):
        """Second upsert on a pending candidate updates confidence_score.
        Uses a non-NULL cancer_type so the UNIQUE(gene, axis, cancer_type)
        constraint fires reliably (SQLite NULLs are never equal in UNIQUE).
        """
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c1 = ReceiptCandidate(**_make_candidate(gene="KRAS", cancer_type="pancreatic", confidence_score=0.60))
        AuditQueue.upsert(c1)
        c2 = ReceiptCandidate(**_make_candidate(gene="KRAS", cancer_type="pancreatic", confidence_score=0.85))
        saved = AuditQueue.upsert(c2)
        assert saved.confidence_score == pytest.approx(0.85)

    def test_upsert_duplicate_approved_does_not_overwrite(self, tmp_db):
        """Re-upserting an approved candidate must NOT overwrite human decision."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c = ReceiptCandidate(**_make_candidate(gene="EGFR", confidence_score=0.70))
        saved = AuditQueue.upsert(c)
        result_approve = AuditQueue.approve(saved.id, "looks good", "dr_test")
        assert result_approve is not None  # approve() returns dict now
        # Now try to re-upsert with different confidence
        c2 = ReceiptCandidate(**_make_candidate(gene="EGFR", confidence_score=0.30))
        result = AuditQueue.upsert(c2)
        assert result.audit_status == "approved"
        assert result.confidence_score == pytest.approx(0.70)  # unchanged

    def test_upsert_duplicate_rejected_does_not_overwrite(self, tmp_db):
        """Re-upserting a rejected candidate must NOT overwrite human decision."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c = ReceiptCandidate(**_make_candidate(gene="ALK", confidence_score=0.65))
        saved = AuditQueue.upsert(c)
        AuditQueue.reject(saved.id, "insufficient evidence", "dr_test")
        c2 = ReceiptCandidate(**_make_candidate(gene="ALK", confidence_score=0.90))
        result = AuditQueue.upsert(c2)
        assert result.audit_status == "rejected"
        assert result.confidence_score == pytest.approx(0.65)  # unchanged

    def test_list_pending_sorted_by_confidence_desc(self, tmp_db):
        """list_pending must return candidates sorted by confidence descending."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        for score, gene in [(0.40, "GENE_A"), (0.90, "GENE_B"), (0.65, "GENE_C")]:
            AuditQueue.upsert(ReceiptCandidate(**_make_candidate(gene=gene, confidence_score=score)))
        results = AuditQueue.list_pending()
        scores = [r.confidence_score for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_list_pending_min_confidence_filter(self, tmp_db):
        """list_pending with min_confidence=0.70 must exclude low-confidence candidates."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        for score, gene in [(0.50, "LOW_CONF"), (0.80, "HIGH_CONF")]:
            AuditQueue.upsert(ReceiptCandidate(**_make_candidate(gene=gene, confidence_score=score)))
        results = AuditQueue.list_pending(min_confidence=0.70)
        genes = [r.gene for r in results]
        assert "HIGH_CONF" in genes
        assert "LOW_CONF" not in genes

    def test_list_pending_gene_filter(self, tmp_db):
        """list_pending with gene filter must only return matching gene."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        AuditQueue.upsert(ReceiptCandidate(**_make_candidate(gene="MBD4", axis="cytidine_analogs")))
        AuditQueue.upsert(ReceiptCandidate(**_make_candidate(gene="BRCA1", axis="parp_inhibitors")))
        results = AuditQueue.list_pending(gene="MBD4")
        assert all(r.gene == "MBD4" for r in results)
        assert len(results) >= 1

    def test_approve_sets_status_and_audited_by(self, tmp_db):
        """approve() must set audit_status='approved' and store audited_by."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c = ReceiptCandidate(**_make_candidate(gene="PIK3CA"))
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, "strong mechanistic evidence", "oncologist_1")
        assert result is not None
        updated = result["candidate"]
        assert updated.audit_status == "approved"
        assert updated.audited_by == "oncologist_1"
        assert updated.audit_notes == "strong mechanistic evidence"
        assert updated.audited_at is not None
        assert "frozen_receipt_snippet" in result
        assert isinstance(result["frozen_receipt_snippet"], str)

    def test_reject_sets_status_and_audited_by(self, tmp_db):
        """reject() must set audit_status='rejected' and store audited_by."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        c = ReceiptCandidate(**_make_candidate(gene="PTEN"))
        saved = AuditQueue.upsert(c)
        updated = AuditQueue.reject(saved.id, "artifact in data", "oncologist_2")
        assert updated.audit_status == "rejected"
        assert updated.audited_by == "oncologist_2"
        assert updated.audit_notes == "artifact in data"

    def test_stats_returns_correct_counts(self, tmp_db):
        """stats() must reflect accurate counts after approve/reject operations."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue
        # Insert 3 candidates
        ids = []
        for gene in ["GENE_X", "GENE_Y", "GENE_Z"]:
            saved = AuditQueue.upsert(ReceiptCandidate(**_make_candidate(
                gene=gene, axis="atr_wee1", confidence_score=0.80
            )))
            ids.append(saved.id)
        # Approve one, reject one, leave one pending
        AuditQueue.approve(ids[0], None, "auditor")
        AuditQueue.reject(ids[1], None, "auditor")
        s = AuditQueue.stats()
        qc = s["queue_counts"]
        assert qc["pending"] == 1
        assert qc["approved"] == 1
        assert qc["rejected"] == 1
        assert qc["promoted"] == 0
        assert qc["high_confidence_pending"] == 1  # 0.80 >= 0.70


# ─────────────────────────────────────────────────────────────────────────────
# TestAuditRoutes (10 tests — FastAPI TestClient)
# ─────────────────────────────────────────────────────────────────────────────

def _seed_payload(**kwargs) -> dict:
    """Build a valid seed request body."""
    defaults = {
        "gene": "MBD4",
        "axis": "cytidine_analogs",
        "confidence_score": 0.88,
        "candidate_tier": "Validated SL therapeutic lever",
        "evidence_summary": "MBD4 loss confers cytidine analog sensitivity (CRISPR delta=-0.35)",
    }
    defaults.update(kwargs)
    return defaults


class TestAuditRoutes:

    def test_get_queue_returns_ruo_envelope(self, client):
        """GET /api/v1/audit/queue must return RUO envelope keys."""
        r = client.get("/api/v1/audit/queue")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        assert "disclaimer" in body
        assert "data" in body
        assert isinstance(body["data"], list)

    def test_get_queue_with_gene_filter(self, client):
        """GET /api/v1/audit/queue?gene=MBD4 must filter by gene."""
        # Seed two genes
        client.post("/api/v1/audit/seed", json=_seed_payload(gene="MBD4", axis="cytidine_analogs"))
        client.post("/api/v1/audit/seed", json=_seed_payload(gene="BRCA2", axis="parp_inhibitors"))
        r = client.get("/api/v1/audit/queue?gene=MBD4")
        assert r.status_code == 200
        data = r.json()["data"]
        assert all(c["gene"] == "MBD4" for c in data)

    def test_get_candidate_returns_candidate(self, client):
        """GET /api/v1/audit/{id} must return the seeded candidate."""
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload(gene="TP53", axis="atr_wee1"))
        assert seed_r.status_code == 200
        candidate_id = seed_r.json()["data"]["id"]
        r = client.get(f"/api/v1/audit/{candidate_id}")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        assert body["data"]["id"] == candidate_id

    def test_get_candidate_missing_returns_404(self, client):
        """GET /api/v1/audit/99999 on non-existent id must return 404."""
        r = client.get("/api/v1/audit/99999")
        assert r.status_code == 404

    def test_post_approve_updates_status(self, client):
        """POST /api/v1/audit/{id}/approve must set status to 'approved'."""
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload(gene="RB1", axis="pkmyt1"))
        candidate_id = seed_r.json()["data"]["id"]
        r = client.post(
            f"/api/v1/audit/{candidate_id}/approve",
            json={"audited_by": "dr_smith", "notes": "confirmed in literature"},
        )
        assert r.status_code == 200
        # Sprint 4: approve returns {status, candidate_id, candidate, frozen_receipt_snippet}
        body = r.json()
        data = body["data"]
        assert data["status"] == "approved"
        assert data["candidate"]["audit_status"] == "approved"
        assert data["candidate"]["audited_by"] == "dr_smith"
        assert "frozen_receipt_snippet" in data
        assert isinstance(data["frozen_receipt_snippet"], str)

    def test_post_reject_updates_status(self, client):
        """POST /api/v1/audit/{id}/reject must set status to 'rejected'."""
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload(gene="VHL", axis="immunotherapy"))
        candidate_id = seed_r.json()["data"]["id"]
        r = client.post(
            f"/api/v1/audit/{candidate_id}/reject",
            json={"audited_by": "dr_jones", "notes": "weak evidence"},
        )
        assert r.status_code == 200
        data = r.json()["data"]
        assert data["audit_status"] == "rejected"
        assert data["audited_by"] == "dr_jones"

    def test_get_stats_returns_correct_keys(self, client):
        """GET /api/v1/audit/stats must return all required stat keys."""
        r = client.get("/api/v1/audit/stats")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        data = body["data"]
        # Sprint 4: stats now nested under queue_counts + coverage_summary + top_pending
        assert "queue_counts" in data
        assert "coverage_summary" in data
        assert "top_pending" in data
        assert "last_mine_at" in data
        assert "coverage_version" in data
        qc = data["queue_counts"]
        for key in ["pending", "approved", "rejected", "high_confidence_pending"]:
            assert key in qc, f"Missing queue_counts key: {key}"

    def test_get_coverage_returns_list(self, client):
        """GET /api/v1/audit/coverage must return coverage list + summary (Sprint 4 shape)."""
        r = client.get("/api/v1/audit/coverage")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        # Sprint 4: body["coverage"] is the list; body["summary"] is the summary block
        assert "coverage" in body, f"Missing 'coverage' key. Keys: {list(body.keys())}"
        assert isinstance(body["coverage"], list)
        assert "summary" in body
        assert "coverage_version" in body

    def test_get_coverage_includes_frozen_receipt_entries(self, client):
        """GET /api/v1/audit/coverage must include entries from _FROZEN_RECEIPTS."""
        r = client.get("/api/v1/audit/coverage")
        assert r.status_code == 200
        # Sprint 4: coverage is under body["coverage"] key
        data = r.json()["coverage"]
        # _FROZEN_RECEIPTS has MBD4 + cytidine_analogs as gold standard
        genes = [e["gene"] for e in data]
        assert "MBD4" in genes
        # At least one entry should have frozen_receipt_count > 0
        assert any(e["frozen_receipt_count"] > 0 for e in data)

    def test_post_seed_injects_candidate(self, client):
        """POST /api/v1/audit/seed must inject candidate and return RUO envelope."""
        payload = _seed_payload(
            gene="FANCA",
            axis="parp_inhibitors",
            confidence_score=0.72,
            candidate_tier="Strong candidate dependency axis",
        )
        r = client.post("/api/v1/audit/seed", json=payload)
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        assert body["data"]["gene"] == "FANCA"
        assert body["data"]["id"] is not None
        # Verify it appears in the queue
        q = client.get("/api/v1/audit/queue?gene=FANCA")
        assert q.status_code == 200
        queue_data = q.json()["data"]
        assert any(c["gene"] == "FANCA" for c in queue_data)

    def test_post_mine_returns_running_status(self, client):
        """
        POST /api/v1/audit/mine must return immediately with status='running'
        and the correct RUO envelope.
        mine_receipts is patched so no real data loading occurs.
        """
        from unittest.mock import patch
        from sl_agent.multimodal.receipt_miner import MinerRunSummary

        mock_summary = MinerRunSummary(gene_count=2, axis_count=1)

        with patch("sl_agent.multimodal.receipt_miner.mine_receipts",
                   return_value=mock_summary):
            r = client.post(
                "/api/v1/audit/mine",
                params={"gene_panel": ["BRCA1", "BRCA2"], "axes": ["parp_inhibitors"]},
            )
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        assert "disclaimer" in body
        data = body["data"]
        assert data["status"] == "running"
        assert data["gene_count"] == 2
        assert data["axis_count"] == 1
        assert data["estimated_pairs"] == 2  # 2 genes * 1 axis

    def test_post_mine_uses_default_panel_when_no_params(self, client):
        """
        POST /api/v1/audit/mine with no params uses load_default_panel() (ddr_panel_v1.json)
        and DEFAULT_AXES.  Verifies gene_count and axis_count match the loaded defaults.
        Sprint 4: gene_count reflects ddr_panel_v1.json (36 genes), not the legacy
        DEFAULT_GENE_PANEL constant (29 genes).
        """
        from unittest.mock import patch
        from sl_agent.multimodal.receipt_miner import (
            MinerRunSummary, DEFAULT_AXES, load_default_panel,
        )

        # Use the same loader that routes.py uses so counts always agree
        expected_panel = load_default_panel()

        mock_summary = MinerRunSummary(
            gene_count=len(expected_panel),
            axis_count=len(DEFAULT_AXES),
        )

        with patch("sl_agent.multimodal.receipt_miner.mine_receipts",
                   return_value=mock_summary):
            r = client.post("/api/v1/audit/mine")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        data = body["data"]
        assert data["status"] == "running"
        # Panel now loaded from ddr_panel_v1.json via load_default_panel()
        assert data["gene_count"] == len(expected_panel)
        assert data["axis_count"] == len(DEFAULT_AXES)
        assert data["estimated_pairs"] == len(expected_panel) * len(DEFAULT_AXES)
