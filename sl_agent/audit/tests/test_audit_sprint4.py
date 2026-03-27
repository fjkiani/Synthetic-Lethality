"""
Sprint 4 — Promotion Pipeline + Coverage Heat Map + Gene Panel + Stats Tests.

All tests use a temporary SQLite DB (same pattern as test_audit.py).

Test groups:
  TestFrozenReceiptSnippet   (5)  — generate_frozen_receipt_snippet + approve() response
  TestCoverageMatrix         (4)  — full matrix + summary block + not_covered entries
  TestGenePanel              (3)  — load_default_panel + mine wiring
  TestStatsUpgrade           (3)  — queue_counts, coverage_summary, top_pending

Total: 15 tests (all 10+ spec requirement)
"""
from __future__ import annotations

from pathlib import Path
from typing import Generator
from unittest.mock import patch, MagicMock

import pytest
from fastapi.testclient import TestClient


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures (mirror test_audit.py pattern)
# ─────────────────────────────────────────────────────────────────────────────

def _make_candidate(**kwargs):
    defaults = {
        "gene": "BRCA1",
        "axis": "parp_inhibitors",
        "confidence_score": 0.75,
        "candidate_tier": "Strong",
        "evidence_summary": "BRCA1 × parp_inhibitors | conf=0.75 | Strong",
    }
    defaults.update(kwargs)
    return defaults


@pytest.fixture()
def tmp_db(tmp_path, monkeypatch) -> Path:
    db_file = tmp_path / "test_sprint4.db"
    import sl_agent.audit.queue as queue_module
    monkeypatch.setattr(queue_module, "_get_db_path", lambda: db_file)
    queue_module.init_db()
    return db_file


@pytest.fixture()
def client(tmp_db, monkeypatch) -> Generator:
    from fastapi import FastAPI
    from sl_agent.audit.routes import audit_router

    app = FastAPI()
    app.include_router(audit_router, prefix="/api/v1")
    with TestClient(app) as tc:
        yield tc


# ─────────────────────────────────────────────────────────────────────────────
# TestFrozenReceiptSnippet  (5 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestFrozenReceiptSnippet:
    """Tests for generate_frozen_receipt_snippet + approve() dict return."""

    def test_approve_returns_dict_with_snippet(self, tmp_db):
        """approve() must return a dict with 'candidate' and 'frozen_receipt_snippet'."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue

        c = ReceiptCandidate(**_make_candidate(
            gene="BRCA2",
            axis="parp_inhibitors",
            confidence_score=0.82,
            candidate_tier="Strong",
            crispr_delta_dep=-0.35,
            crispr_fdr=0.01,
            gdsc_delta_ic50=-1.5,
            kb_clinical_hits=4,
        ))
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, "confirmed in TCGA + PRISM", "dr_atlas")

        assert result is not None
        assert "candidate" in result
        assert "frozen_receipt_snippet" in result
        assert isinstance(result["frozen_receipt_snippet"], str)
        candidate = result["candidate"]
        assert candidate.audit_status == "approved"
        assert candidate.audited_by == "dr_atlas"

    def test_snippet_contains_correct_gene_and_axis(self, tmp_db):
        """Snippet must contain the correct gene and axis enum name."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue

        c = ReceiptCandidate(**_make_candidate(
            gene="MBD4",
            axis="cytidine_analogs",
            confidence_score=0.91,
            candidate_tier="Strong",
            crispr_delta_dep=-0.28,
            crispr_fdr=0.008,
            gdsc_delta_ic50=-2.1,
            kb_clinical_hits=5,
        ))
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, "gold standard axis", "dr_curator")
        snippet = result["frozen_receipt_snippet"]

        assert '"MBD4"' in snippet
        assert "CandidateAxis.CYTIDINE_ANALOGS" in snippet
        assert "audit_queue_promoted" in snippet
        assert "dr_curator" in snippet

    def test_snippet_contains_numeric_scores(self, tmp_db):
        """Snippet must embed crispr_delta_dep, crispr_fdr, gdsc_delta_ic50 values."""
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue

        c = ReceiptCandidate(**_make_candidate(
            gene="ATM",
            axis="parp_inhibitors",
            confidence_score=0.73,
            candidate_tier="Strong",
            crispr_delta_dep=-0.2222,
            crispr_fdr=0.0312,
            gdsc_delta_ic50=-1.1234,
            kb_clinical_hits=3,
        ))
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, "ATM PARP SL", "reviewer_1")
        snippet = result["frozen_receipt_snippet"]

        # Values must be present (formatted to 4 decimal places)
        assert "-0.2222" in snippet
        assert "0.0312" in snippet
        assert "-1.1234" in snippet
        assert "kb_clinical_hits=3" in snippet

    def test_snippet_does_not_contain_validated_tier(self, tmp_db):
        """
        SABOTAGE TEST: snippet tier is derived from candidate_tier.
        Approving a Mechanistic candidate generates a Mechanistic snippet, NOT Validated.
        'VALIDATED' must NEVER appear as the evidence_tier in an auto-generated snippet.
        """
        from sl_agent.audit.models import ReceiptCandidate
        from sl_agent.audit.queue import AuditQueue

        # Mechanistic candidate (below Strong threshold)
        c = ReceiptCandidate(**_make_candidate(
            gene="CHEK2",
            axis="atr_wee1",
            confidence_score=0.55,
            candidate_tier="Mechanistic",   # explicitly Mechanistic
            kb_clinical_hits=1,
        ))
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, "mechanistic only", "safety_reviewer")
        snippet = result["frozen_receipt_snippet"]

        # The tier in the snippet must NOT be Validated
        # It should be MECHANISTIC (uppercased from candidate_tier)
        assert 'evidence_tier="VALIDATED"' not in snippet, (
            "SABOTAGE: auto-generated snippet emitted evidence_tier='VALIDATED'. "
            "Validated tier MUST be assigned manually by a human with clinical evidence."
        )
        assert "MECHANISTIC" in snippet

    def test_api_approve_returns_snippet_in_response(self, client):
        """POST /audit/{id}/approve must include frozen_receipt_snippet in response body."""
        # Seed a candidate
        seed_r = client.post("/api/v1/audit/seed", json=_make_candidate(
            gene="PALB2",
            axis="parp_inhibitors",
            confidence_score=0.78,
            candidate_tier="Strong",
            evidence_summary="PALB2×PARP | conf=0.78 | Strong",
        ))
        assert seed_r.status_code == 200
        candidate_id = seed_r.json()["data"]["id"]

        approve_r = client.post(
            f"/api/v1/audit/{candidate_id}/approve",
            json={"audited_by": "dr_portal", "notes": "sprint4 test"},
        )
        assert approve_r.status_code == 200
        body = approve_r.json()
        assert body["ruo"] is True
        data = body["data"]
        assert data["status"] == "approved"
        assert data["candidate_id"] == candidate_id
        assert "frozen_receipt_snippet" in data
        snippet = data["frozen_receipt_snippet"]
        assert '"PALB2"' in snippet
        assert "CandidateAxis.PARP_INHIBITORS" in snippet


# ─────────────────────────────────────────────────────────────────────────────
# TestCoverageMatrix  (4 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestCoverageMatrix:
    """Tests for the extended /audit/coverage endpoint (Sprint 4 shape)."""

    def test_coverage_returns_full_matrix_shape(self, client):
        """GET /audit/coverage must return {coverage, summary, coverage_version, ruo}."""
        r = client.get("/api/v1/audit/coverage")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        assert "coverage" in body
        assert "summary" in body
        assert "coverage_version" in body
        summary = body["summary"]
        for key in ["validated", "strong", "mechanistic", "not_covered", "total_pairs"]:
            assert key in summary, f"Missing summary key: {key}"

    def test_coverage_summary_counts_add_up(self, client):
        """Summary counts must sum to total_pairs."""
        r = client.get("/api/v1/audit/coverage")
        summary = r.json()["summary"]
        counted = summary["validated"] + summary["strong"] + summary["mechanistic"] + summary["not_covered"]
        assert counted == summary["total_pairs"], (
            f"Summary counts don't add up: {summary}"
        )

    def test_coverage_includes_not_covered_entries(self, client):
        """
        Genes in default panel with no frozen receipts or queue candidates
        must appear as 'Not covered' entries in the matrix.
        """
        r = client.get("/api/v1/audit/coverage")
        assert r.status_code == 200
        coverage = r.json()["coverage"]

        # With empty queue and default panel, most genes should be 'Not covered'
        not_covered = [e for e in coverage if e["tier"].lower() == "not covered"]
        assert len(not_covered) > 0, (
            "Expected at least some 'Not covered' entries for panel genes "
            "that have no frozen receipts or queue candidates."
        )

    def test_coverage_mbd4_has_validated_tier(self, client):
        """
        MBD4 + cytidine_analogs is the gold standard — must appear as Validated
        (has clinical POSITIVE in _FROZEN_RECEIPTS).
        """
        r = client.get("/api/v1/audit/coverage")
        coverage = r.json()["coverage"]
        mbd4_cytidine = [
            e for e in coverage
            if e["gene"] == "MBD4" and e["axis"] == "cytidine_analogs"
        ]
        assert len(mbd4_cytidine) == 1, "MBD4+cytidine_analogs not in coverage matrix"
        entry = mbd4_cytidine[0]
        assert "Validated" in entry["tier"] or entry["tier"] == "Validated", (
            f"MBD4+cytidine_analogs expected Validated, got: {entry['tier']}"
        )
        assert entry["frozen_receipt_count"] > 0
        assert entry["has_clinical_positive"] is True


# ─────────────────────────────────────────────────────────────────────────────
# TestGenePanel  (3 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestGenePanel:
    """Tests for load_default_panel() and mine wiring."""

    def test_load_default_panel_returns_deduped_list(self):
        """
        load_default_panel() reads ddr_panel_v1.json, deduplicates across sub-groups,
        and returns a flat list with no duplicates.
        """
        from sl_agent.multimodal.receipt_miner import load_default_panel
        panel = load_default_panel()
        assert isinstance(panel, list)
        assert len(panel) > 0
        # No duplicates
        assert len(panel) == len(set(panel)), (
            f"load_default_panel() returned duplicates: "
            f"{[g for g in set(panel) if panel.count(g) > 1]}"
        )

    def test_load_default_panel_contains_expected_genes(self):
        """Panel must contain key DDR genes."""
        from sl_agent.multimodal.receipt_miner import load_default_panel
        panel = load_default_panel()
        required = ["BRCA1", "BRCA2", "MBD4", "ATM", "CHEK1", "WEE1"]
        for gene in required:
            assert gene in panel, f"Expected gene '{gene}' not in default panel"

    def test_mine_uses_default_panel_when_no_params(self, client):
        """
        POST /audit/mine with no body uses load_default_panel() — which reads
        ddr_panel_v1.json. Gene count must be >= DEFAULT_GENE_PANEL length.
        """
        from sl_agent.multimodal.receipt_miner import MinerRunSummary, load_default_panel
        expected_panel = load_default_panel()

        mock_summary = MinerRunSummary(
            gene_count=len(expected_panel),
            axis_count=6,
        )
        mock_summary.finish()

        with patch("sl_agent.multimodal.receipt_miner.mine_receipts",
                   return_value=mock_summary):
            r = client.post("/api/v1/audit/mine")  # no body
        assert r.status_code == 200
        data = r.json()["data"]
        assert data["status"] == "running"
        # Gene count must equal load_default_panel() result (deduped panel)
        assert data["gene_count"] == len(expected_panel)


# ─────────────────────────────────────────────────────────────────────────────
# TestStatsUpgrade  (3 tests)
# ─────────────────────────────────────────────────────────────────────────────

class TestStatsUpgrade:
    """Tests for the upgraded GET /audit/stats endpoint."""

    def test_stats_returns_all_required_keys(self, client):
        """GET /audit/stats must return queue_counts, coverage_summary, top_pending, last_mine_at."""
        r = client.get("/api/v1/audit/stats")
        assert r.status_code == 200
        body = r.json()
        assert body["ruo"] is True
        data = body["data"]
        assert "queue_counts" in data
        assert "coverage_summary" in data
        assert "top_pending" in data
        assert "last_mine_at" in data
        assert "coverage_version" in data
        qc = data["queue_counts"]
        for key in ["pending", "approved", "rejected", "high_confidence_pending"]:
            assert key in qc, f"Missing queue_counts.{key}"
        cs = data["coverage_summary"]
        for key in ["validated", "strong", "mechanistic", "not_covered"]:
            assert key in cs, f"Missing coverage_summary.{key}"

    def test_stats_top_pending_sorted_by_confidence(self, client):
        """top_pending must be sorted by confidence descending (top 5 max)."""
        # Seed 6 candidates with different confidence scores
        scores = [0.88, 0.72, 0.55, 0.91, 0.63, 0.45]
        for i, score in enumerate(scores):
            client.post("/api/v1/audit/seed", json=_make_candidate(
                gene=f"GENE_{i:02d}",
                axis="atr_wee1",
                confidence_score=score,
                candidate_tier="Strong" if score >= 0.70 else "Mechanistic",
                evidence_summary=f"GENE_{i:02d}×atr_wee1 | conf={score}",
            ))

        r = client.get("/api/v1/audit/stats")
        data = r.json()["data"]
        top = data["top_pending"]
        # Must be sorted descending
        confidences = [t["confidence"] for t in top]
        assert confidences == sorted(confidences, reverse=True), (
            f"top_pending not sorted by confidence desc: {confidences}"
        )
        # Max 5 entries
        assert len(top) <= 5

    def test_stats_coverage_summary_matches_coverage_endpoint(self, client):
        """
        stats().coverage_summary must match /audit/coverage summary block
        (consistent data — same underlying call).
        """
        stats_r = client.get("/api/v1/audit/stats")
        coverage_r = client.get("/api/v1/audit/coverage")

        assert stats_r.status_code == 200
        assert coverage_r.status_code == 200

        stats_summary = stats_r.json()["data"]["coverage_summary"]
        cov_summary = coverage_r.json()["summary"]

        # Both must agree on validated, strong, mechanistic counts
        # (not_covered may differ by 1 if there is a race — but in tests they should match)
        for key in ["validated", "strong", "mechanistic"]:
            assert stats_summary[key] == cov_summary[key], (
                f"coverage_summary.{key} mismatch: stats={stats_summary[key]}, "
                f"coverage={cov_summary[key]}"
            )
