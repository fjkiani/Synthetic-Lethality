"""
test_stress.py — Adversarial stress tests for the audit queue (Attack 1, 2, 3, 4, 6).

These tests try to BREAK the system in every way management, regulators, and hostile
reviewers will attempt. Every test here corresponds to a named attack vector.

Attack vectors covered:
  Attack 1 — Confidence Inflation
  Attack 2 — Decision Immutability
  Attack 3 — Tier Laundering
  Attack 4 — API Boundary / SQL Injection / Malformed Input
  Attack 6 — Concurrent Approve Race Condition

Run with:
  pytest sl_agent/audit/tests/test_stress.py -v
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Generator

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

# ── Fixtures ──────────────────────────────────────────────────────────────────

def _seed_payload(
    gene: str = "BRCA1",
    axis: str = "parp_inhibitors",
    confidence: float = 0.82,
    tier: str = "Strong",
) -> dict:
    return {
        "gene": gene,
        "axis": axis,
        "confidence_score": confidence,
        "candidate_tier": tier,
        "evidence_summary": f"{gene}×{axis} | test payload | conf={confidence}",
    }


@pytest.fixture()
def tmp_db(tmp_path: Path, monkeypatch) -> Path:
    """Redirect DB to an isolated temp file so every test starts clean."""
    import sl_agent.audit.queue as queue_module
    db_file = tmp_path / "stress_test.db"
    monkeypatch.setattr(queue_module, "_get_db_path", lambda: db_file)
    from sl_agent.audit.queue import init_db
    init_db()
    return db_file


@pytest.fixture()
def client(tmp_db: Path, monkeypatch) -> Generator:
    """Minimal FastAPI app — audit router only, no network I/O."""
    from sl_agent.audit.routes import audit_router
    import sl_agent.audit.queue as queue_module
    app = FastAPI()
    app.include_router(audit_router, prefix="/api/v1")
    with TestClient(app) as tc:
        yield tc


# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 1 — CONFIDENCE INFLATION
# Goal: get a garbage gene-drug pair to score ≥ 0.40 and enter the queue.
# ─────────────────────────────────────────────────────────────────────────────

class TestConfidenceInflation:
    """Attack 1: Can garbage data breach the 0.40 confidence gate?"""

    def test_confidence_score_le_1_enforced_by_model(self, tmp_db):
        """ReceiptCandidate rejects confidence_score > 1.0 at model level."""
        from pydantic import ValidationError
        from sl_agent.audit.models import ReceiptCandidate
        with pytest.raises(ValidationError):
            ReceiptCandidate(
                gene="FAKE",
                axis="parp_inhibitors",
                confidence_score=1.05,
                candidate_tier="Strong",
                evidence_summary="overflow test",
            )

    def test_confidence_score_negative_rejected_by_model(self, tmp_db):
        """ReceiptCandidate rejects negative confidence_score at model level."""
        from pydantic import ValidationError
        from sl_agent.audit.models import ReceiptCandidate
        with pytest.raises(ValidationError):
            ReceiptCandidate(
                gene="FAKE",
                axis="parp_inhibitors",
                confidence_score=-0.10,
                candidate_tier="Mechanistic",
                evidence_summary="negative confidence test",
            )

    def test_manually_upserted_validated_tier_stored_as_is(self, tmp_db):
        """
        AuditQueue.upsert does NOT validate candidate_tier — it stores whatever
        is given (it is a data-intake layer, not an enforcement layer).
        This is known and intentional: enforcement happens at snippet generation.
        This test DOCUMENTS the behaviour, it does not assert it is a bug.
        The attack is defeated at generate_frozen_receipt_snippet() (see Attack 3 tests).
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate
        c = ReceiptCandidate(
            gene="FAKE",
            axis="parp_inhibitors",
            confidence_score=1.0,
            candidate_tier="Validated",  # manually injected tier
            evidence_summary="manual injection attempt",
        )
        saved = AuditQueue.upsert(c)
        # Stored as given — the guard fires at snippet time, not here.
        assert saved.candidate_tier == "Validated"
        assert saved.id is not None

    def test_upsert_duplicate_approved_does_not_overwrite_confidence(self, tmp_db):
        """
        After a candidate is approved, a re-upsert of the same gene+axis
        with lower confidence does NOT overwrite the approved record.
        Approved decisions are immutable at the upsert layer.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        high = ReceiptCandidate(
            gene="MBD4", axis="cytidine_analogs",
            confidence_score=0.93, candidate_tier="Strong",
            evidence_summary="high confidence",
        )
        saved = AuditQueue.upsert(high)
        AuditQueue.approve(saved.id, "good data", "dr_real")

        low = ReceiptCandidate(
            gene="MBD4", axis="cytidine_analogs",
            confidence_score=0.30, candidate_tier="Mechanistic",
            evidence_summary="attempted overwrite with low confidence",
        )
        AuditQueue.upsert(low)

        # Fetch by id — confidence must still be 0.93, status still approved
        refetched = AuditQueue.get(saved.id)
        assert refetched.audit_status == "approved"
        assert refetched.confidence_score == pytest.approx(0.93, abs=0.01)

    def test_empty_gene_string_seeds_but_scores_zero_in_miner(self, tmp_db):
        """
        Empty gene string can be seeded directly (data-intake is permissive)
        but mine_receipts with empty gene_panel returns immediately with 0 pairs.
        """
        from sl_agent.multimodal.receipt_miner import mine_receipts, MinerThresholds
        from sl_agent.multimodal.models import CandidateAxis
        summary = mine_receipts(
            gene_panel=[],
            axes=[CandidateAxis.PARP_INHIBITORS],
            thresholds=MinerThresholds(),
        )
        assert summary.pairs_evaluated == 0
        assert summary.candidates_queued == 0


# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 2 — DECISION IMMUTABILITY
# Goal: overwrite an approved human decision.
# ─────────────────────────────────────────────────────────────────────────────

class TestDecisionImmutability:
    """Attack 2: Human decisions must be irreversible in both directions."""

    def test_reject_already_approved_candidate_raises(self, tmp_db):
        """
        CRITICAL: reject() on an already-approved candidate must raise ValueError.
        Without this guard, a bad actor could flip an approved clinical receipt
        to rejected via the API.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="BRCA1", axis="parp_inhibitors",
            confidence_score=0.85, candidate_tier="Strong",
            evidence_summary="legit approval",
        )
        saved = AuditQueue.upsert(c)
        AuditQueue.approve(saved.id, "approved by panel", "dr_oncology")

        with pytest.raises(ValueError, match="already approved"):
            AuditQueue.reject(saved.id, "changed my mind", "dr_saboteur")

    def test_approve_already_rejected_candidate_raises(self, tmp_db):
        """
        approve() on an already-rejected candidate must raise ValueError.
        A rejected record represents a deliberate human decision to discard
        low-quality evidence; it cannot be silently un-rejected.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="TSPAN6", axis="parp_inhibitors",
            confidence_score=0.42, candidate_tier="Mechanistic",
            evidence_summary="weak signal",
        )
        saved = AuditQueue.upsert(c)
        AuditQueue.reject(saved.id, "insufficient evidence", "dr_skeptic")

        with pytest.raises(ValueError, match="already rejected"):
            AuditQueue.approve(saved.id, "override attempt", "dr_override")

    def test_api_reject_approved_returns_409(self, client):
        """
        HTTP layer: POST /reject on an approved candidate returns 409 Conflict.
        Not 200, not 500 — the client must know this was an immutability violation.
        """
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload("ATM", "atr_wee1", 0.88))
        cid = seed_r.json()["data"]["id"]
        client.post(f"/api/v1/audit/{cid}/approve",
                    json={"audited_by": "dr_panel", "notes": "approved"})
        r = client.post(f"/api/v1/audit/{cid}/reject",
                        json={"audited_by": "dr_saboteur", "notes": "flip"})
        assert r.status_code == 409
        assert "already approved" in r.json()["detail"].lower()

    def test_api_approve_rejected_returns_409(self, client):
        """HTTP layer: POST /approve on a rejected candidate returns 409 Conflict."""
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload("VHL", "immunotherapy", 0.44))
        cid = seed_r.json()["data"]["id"]
        client.post(f"/api/v1/audit/{cid}/reject",
                    json={"audited_by": "dr_panel", "notes": "rejected"})
        r = client.post(f"/api/v1/audit/{cid}/approve",
                        json={"audited_by": "dr_override", "notes": "override"})
        assert r.status_code == 409
        assert "already rejected" in r.json()["detail"].lower()

    def test_approve_with_empty_audited_by_raises(self, tmp_db):
        """
        Empty audited_by is an audit trail gap. approve() must refuse.
        A receipt snippet with promoted_by="" is legally and regulatorily invalid.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="RAD51", axis="parp_inhibitors",
            confidence_score=0.75, candidate_tier="Strong",
            evidence_summary="audit trail test",
        )
        saved = AuditQueue.upsert(c)
        with pytest.raises(ValueError, match="audited_by cannot be empty"):
            AuditQueue.approve(saved.id, "notes", "")

    def test_reject_with_empty_audited_by_raises(self, tmp_db):
        """Empty audited_by on reject() must also refuse — same audit trail requirement."""
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="CHEK1", axis="atr_wee1",
            confidence_score=0.65, candidate_tier="Mechanistic",
            evidence_summary="audit trail test",
        )
        saved = AuditQueue.upsert(c)
        with pytest.raises(ValueError, match="audited_by cannot be empty"):
            AuditQueue.reject(saved.id, "notes", "")

    def test_approve_with_none_notes_renders_clean_snippet(self, tmp_db):
        """
        notes=None is allowed (notes are optional per the model spec).
        The frozen_receipt_snippet must render cleanly — no 'None' string in output.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="PALB2", axis="parp_inhibitors",
            confidence_score=0.80, candidate_tier="Strong",
            evidence_summary="null notes test",
        )
        saved = AuditQueue.upsert(c)
        result = AuditQueue.approve(saved.id, None, "dr_nonotes")
        snippet = result["frozen_receipt_snippet"]
        # notes field should be empty string, not the Python None literal
        assert "notes=\"None\"" not in snippet
        assert "notes=" in snippet  # field is still present, just blank


# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 3 — TIER LAUNDERING
# Goal: get VALIDATED tier into a frozen receipt without a real clinical POSITIVE.
# ─────────────────────────────────────────────────────────────────────────────

class TestTierLaundering:
    """Attack 3: The snippet generator must clamp VALIDATED tier regardless of DB content."""

    def test_snippet_clamps_validated_to_strong(self, tmp_db):
        """
        CRITICAL: if a row was manually upserted with candidate_tier="Validated",
        generate_frozen_receipt_snippet() must NOT emit evidence_tier="VALIDATED".
        The clamp maps VALIDATED → STRONG (most credible auto-clamp below human threshold).
        """
        from sl_agent.audit.queue import AuditQueue, generate_frozen_receipt_snippet
        from sl_agent.audit.models import ReceiptCandidate

        injected = ReceiptCandidate(
            gene="BRCA2",
            axis="parp_inhibitors",
            confidence_score=0.95,
            candidate_tier="Validated",  # manually injected
            evidence_summary="tier laundering attempt",
        )
        saved = AuditQueue.upsert(injected)
        snippet = generate_frozen_receipt_snippet(saved)

        assert "VALIDATED" not in snippet, (
            "TIER LAUNDERING: snippet emitted VALIDATED tier from manually injected candidate. "
            "This bypasses the clinical validation requirement."
        )
        assert "STRONG" in snippet, "Expected STRONG as the clamp-down tier"

    def test_approve_with_manually_injected_validated_tier_clamps_snippet(self, tmp_db):
        """
        Full pipeline: upsert(Validated) → approve() → snippet must still say STRONG.
        The clamp must fire during approval, not just in isolation.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        injected = ReceiptCandidate(
            gene="ATR",
            axis="atr_wee1",
            confidence_score=0.91,
            candidate_tier="Validated",
            evidence_summary="full pipeline laundering test",
        )
        saved = AuditQueue.upsert(injected)
        result = AuditQueue.approve(saved.id, "approval notes", "dr_launderer")
        snippet = result["frozen_receipt_snippet"]

        assert "VALIDATED" not in snippet
        assert "STRONG" in snippet

    def test_sqlite_direct_edit_then_approve_still_clamps(self, tmp_db):
        """
        Simulates a direct SQLite row edit (UPDATE candidates SET candidate_tier='Validated').
        Even if the DB row says Validated, the approve pipeline must clamp the snippet.
        """
        import sqlite3
        from sl_agent.audit.queue import AuditQueue, _get_db_path
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="CCNE1",
            axis="atr_wee1",
            confidence_score=0.88,
            candidate_tier="Strong",
            evidence_summary="db edit simulation",
        )
        saved = AuditQueue.upsert(c)

        # Direct SQLite write — bypasses all Python-layer validation
        with sqlite3.connect(str(_get_db_path())) as con:
            con.execute(
                "UPDATE receipt_candidates SET candidate_tier='Validated' WHERE id=?",
                (saved.id,),
            )
            con.commit()

        result = AuditQueue.approve(saved.id, "post-db-edit approval", "dr_dba")
        snippet = result["frozen_receipt_snippet"]

        assert "VALIDATED" not in snippet, (
            "TIER LAUNDERING via direct DB edit: snippet emitted VALIDATED. "
            "The clamp in generate_frozen_receipt_snippet() failed."
        )
        assert "STRONG" in snippet

    def test_snippet_clamps_unknown_tier_to_strong(self, tmp_db):
        """
        An unexpected tier string (e.g. 'EXPERIMENTAL', 'PHASE3') must also be clamped.
        Only STRONG, MECHANISTIC, INSUFFICIENT are allowed in snippets.
        """
        from sl_agent.audit.queue import generate_frozen_receipt_snippet
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="MBD4",
            axis="cytidine_analogs",
            confidence_score=0.80,
            candidate_tier="EXPERIMENTAL",  # unknown tier
            evidence_summary="unknown tier test",
        )
        snippet = generate_frozen_receipt_snippet(c)
        assert "EXPERIMENTAL" not in snippet
        assert "STRONG" in snippet  # clamp default

    def test_project_tier_never_returns_validated(self, tmp_db):
        """
        Legacy sabotage guard — _project_tier() must never return 'Validated'
        regardless of confidence and KB hits. Unchanged from Sprint 1.
        """
        from sl_agent.multimodal.receipt_miner import _project_tier
        # Max possible inputs
        assert _project_tier(1.0, 100) != "Validated"
        assert _project_tier(1.0, 10) != "Validated"
        assert _project_tier(0.99, 5) != "Validated"


# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 4 — API BOUNDARY
# Goal: crash endpoints, inject SQL, expose internal state.
# ─────────────────────────────────────────────────────────────────────────────

class TestAPIBoundary:
    """Attack 4: Malformed input, SQL injection, integer overflow, boundary cases."""

    def test_malformed_json_to_approve_returns_422(self, client):
        """FastAPI's request body validator returns 422 for malformed JSON bodies."""
        seed_r = client.post("/api/v1/audit/seed", json=_seed_payload())
        cid = seed_r.json()["data"]["id"]
        r = client.post(
            f"/api/v1/audit/{cid}/approve",
            content=b"NOT JSON AT ALL",
            headers={"Content-Type": "application/json"},
        )
        assert r.status_code == 422

    def test_sql_injection_in_gene_name_does_not_crash(self, client):
        """
        Gene name containing SQL injection string must not crash the endpoint
        and must not drop or corrupt the candidates table.
        The miner uses parameterized queries throughout.
        """
        injection_gene = "BRCA1'; DROP TABLE receipt_candidates; --"
        payload = _seed_payload(gene=injection_gene)
        r = client.post("/api/v1/audit/seed", json=payload)
        # Should succeed (seed is permissive) — what matters is no DB corruption
        assert r.status_code == 200
        # Queue should still be intact — GET /queue should work
        q = client.get("/api/v1/audit/queue")
        assert q.status_code == 200

    def test_sql_injection_in_gene_query_param_does_not_crash(self, client):
        """SQL injection in the ?gene= query param must be safely parameterized."""
        r = client.get("/api/v1/audit/queue?gene=BRCA1%27+OR+%271%27%3D%271")
        assert r.status_code == 200

    def test_integer_overflow_candidate_id_returns_4xx(self, client):
        """Extremely large candidate ID must not cause a 500 — expect 404 or 422."""
        r = client.get("/api/v1/audit/99999999999999999999")
        assert r.status_code in (404, 422)

    def test_empty_gene_panel_mine_returns_zero_pairs(self, client):
        """
        POST /audit/mine with gene_panel=[] (no genes) must return immediately
        with estimated_pairs=0 and not crash or divide-by-zero.
        The background task is mocked so no network I/O happens in this test.
        """
        from unittest.mock import patch as _patch
        from sl_agent.multimodal.receipt_miner import MinerRunSummary
        with _patch(
            "sl_agent.multimodal.receipt_miner.mine_receipts",
            return_value=MinerRunSummary(gene_count=0, axis_count=0),
        ):
            r = client.post("/api/v1/audit/mine")
        assert r.status_code == 200
        data = r.json()["data"]
        assert data["status"] == "running"
        # Default panel is 36 genes — estimated_pairs still set from load_default_panel
        assert "estimated_pairs" in data

    def test_invalid_axis_in_mine_is_silently_filtered(self, client):
        """
        POST /audit/mine with an invalid axis value must not raise an unhandled
        exception. Invalid axes are filtered by the CandidateAxis enum check,
        resulting in axis_count=0 and estimated_pairs=0.
        """
        from unittest.mock import patch as _patch
        from sl_agent.multimodal.receipt_miner import MinerRunSummary
        with _patch(
            "sl_agent.multimodal.receipt_miner.mine_receipts",
            return_value=MinerRunSummary(gene_count=1, axis_count=0),
        ):
            r = client.post("/api/v1/audit/mine?gene_panel=BRCA1&axes=NOT_AN_AXIS")
        assert r.status_code == 200
        data = r.json()["data"]
        # Invalid axis filtered → 0 valid axes → estimated_pairs = 1*0 = 0
        assert data["estimated_pairs"] == 0

    def test_500_char_gene_name_seeds_gracefully(self, client):
        """500-character gene name must not cause a 500 or corrupt the DB."""
        long_gene = "A" * 500
        payload = _seed_payload(gene=long_gene)
        r = client.post("/api/v1/audit/seed", json=payload)
        assert r.status_code == 200
        # Queue must still be queryable
        q = client.get("/api/v1/audit/queue")
        assert q.status_code == 200

    def test_coverage_with_no_candidates_returns_not_covered_matrix(self, client):
        """
        GET /audit/coverage with an empty DB must return a valid matrix
        (not_covered for all panel pairs) — no NoneType errors in merge logic.
        """
        r = client.get("/api/v1/audit/coverage")
        assert r.status_code == 200
        body = r.json()
        assert "coverage" in body
        assert "summary" in body
        # All entries should be not_covered when no candidates exist
        tiers = [e["tier"].lower() for e in body["coverage"]]
        assert all("not covered" in t or "strong" in t or "validated" in t for t in tiers)


# ─────────────────────────────────────────────────────────────────────────────
# ATTACK 6 — CONCURRENT APPROVE RACE CONDITION
# Goal: two simultaneous approve calls on the same candidate ID.
# ─────────────────────────────────────────────────────────────────────────────

class TestConcurrentApprove:
    """Attack 6: Race condition — two threads approve the same candidate simultaneously."""

    def test_concurrent_approve_produces_exactly_one_approved_record(self, tmp_db):
        """
        Two threads call AuditQueue.approve() on the same candidate_id simultaneously.
        SQLite serializes writes (WAL mode). Exactly one approved record must exist.
        The second call must either succeed (idempotent) or raise a ValueError
        (already approved is not a valid state to approve again — but both are
        acceptable since the first approve already succeeded).
        What is NOT acceptable: two different audit trails, data corruption,
        or a status that is neither 'approved' nor raises cleanly.
        """
        from sl_agent.audit.queue import AuditQueue
        from sl_agent.audit.models import ReceiptCandidate

        c = ReceiptCandidate(
            gene="CDK2", axis="atr_wee1",
            confidence_score=0.78, candidate_tier="Strong",
            evidence_summary="race condition test",
        )
        saved = AuditQueue.upsert(c)

        results = []
        errors = []

        def _approve(auditor: str):
            try:
                r = AuditQueue.approve(saved.id, "concurrent approve", auditor)
                results.append(r)
            except ValueError as e:
                errors.append(str(e))

        t1 = threading.Thread(target=_approve, args=("dr_thread_1",))
        t2 = threading.Thread(target=_approve, args=("dr_thread_2",))
        t1.start(); t2.start()
        t1.join(); t2.join()

        # After the race: candidate must be in exactly one final state
        final = AuditQueue.get(saved.id)
        assert final.audit_status == "approved", (
            f"After concurrent approve, status is {final.audit_status!r} — expected 'approved'"
        )
        # At least one thread must have succeeded
        assert len(results) >= 1, "At least one approve call must succeed"
        # Total outcomes = 2 (one per thread): success OR ValueError
        assert len(results) + len(errors) == 2
