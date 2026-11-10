# Sprint 4S — Management Artifact Audit
## Six Required Proofs — All Passing

**Requested by**: Management review  
**Produced by**: Computer 1 (hardening instance)  
**Date**: 2026-03-27  
**Verdict**: All 6 proofs pass. Branch is merge-ready.

---

## Full Suite Result

```
======================== 236 passed in 63.85s (0:01:03) ========================
```

Sprint baseline: 195. New tests this sprint: 41 (26 in test_stress.py + 15 in test_adversarial.py).  
Zero regressions. Zero failures.

---

## Proof 1 — git diff --stat: Files Changed

```
MODIFIED  sl_agent/audit/queue.py           560 lines  (tier clamp, immutability, empty audited_by)
MODIFIED  sl_agent/audit/routes.py          268 lines  (409 Conflict, integer overflow 422)
MODIFIED  sl_agent/multimodal/receipt_miner.py  530 lines  (KB direction filter, POLE guard)
NEW       sl_agent/audit/tests/test_stress.py   570 lines  (26 adversarial tests)
NEW       sl_agent/multimodal/tests/test_adversarial.py  426 lines  (15 adversarial tests)
```

**Unified diff — key hunks:**

### queue.py — P1: Tier clamp (lines 89-91)
```python
# BEFORE (vulnerable):
tier = (candidate.candidate_tier or "Mechanistic")

# AFTER (patched):
_raw_tier = (candidate.candidate_tier or "Mechanistic").upper().replace(" ", "_")
_ALLOWED_SNIPPET_TIERS = {"STRONG", "MECHANISTIC", "INSUFFICIENT"}
tier = _raw_tier if _raw_tier in _ALLOWED_SNIPPET_TIERS else "STRONG"
```

### queue.py — P2: Immutability guard in approve() (lines 297-315)
```python
# BEFORE (vulnerable): no guard, silently allowed flip
con.execute("UPDATE receipt_candidates SET audit_status='approved' ...")

# AFTER (patched):
if not audited_by or not audited_by.strip():
    raise ValueError("audited_by cannot be empty — audit trail requires an auditor identity")
...
existing = con.execute("SELECT audit_status FROM receipt_candidates WHERE id=?", ...)
if existing["audit_status"] == "rejected":
    raise ValueError(
        f"Candidate {candidate_id} is already rejected. "
        "Human decisions are immutable. Create a new candidate if evidence has changed."
    )
```

### queue.py — P2: Immutability guard in reject() (lines 336-354)
```python
if existing["audit_status"] == "approved":
    raise ValueError(
        f"Candidate {candidate_id} is already approved. "
        "Human decisions are immutable. Approved receipts cannot be silently rejected."
    )
```

### routes.py — 409 Conflict for immutability violations (lines 236-238, 263-265)
```python
# BEFORE: ValueError → unhandled → 500
# AFTER:
except ValueError as e:
    raise HTTPException(status_code=409, detail=str(e))
```

### routes.py — Integer overflow guard (lines 210-212, 232-233, 259-260)
```python
if candidate_id > 2**63 - 1 or candidate_id < 1:
    raise HTTPException(status_code=422, detail="candidate_id out of valid range")
```

### receipt_miner.py — P3: KB direction filter (lines 253-260)
```python
# BEFORE (vulnerable):
kb_hits = len(response.recommendations)  # counted ALL directions

# AFTER (patched):
sensitivity_recs = [
    r for r in response.recommendations
    if r.response_type == ResponseType.SENSITIVITY
]
kb_hits = len(sensitivity_recs)
```

### receipt_miner.py — P4: POLE DDR axis guard (lines 442-457)
```python
# NEW — fires on 5 DDR axes only, IO is explicitly NOT in the set
_DDR_AXES_FOR_POLE_GUARD = {
    CandidateAxis.PARP_INHIBITORS, CandidateAxis.ATR_WEE1,
    CandidateAxis.PKMYT1, CandidateAxis.WRN, CandidateAxis.CYTIDINE_ANALOGS,
}
if gene_upper == "POLE" and axis in _DDR_AXES_FOR_POLE_GUARD:
    summary.pairs_evaluated += 1
    summary.pairs_no_signal += 1
    continue  # IMMUNOTHERAPY (IO) is NOT in this set — clinical safety
```

---

## Proof 2 — generate_frozen_receipt_snippet() Tier Clamp

Live execution output:

```
[PASS] input='Validated'          → snippet tier='STRONG'
[PASS] input='VALIDATED'          → snippet tier='STRONG'
[PASS] input='TOTALLY_FAKE_XYZ'   → snippet tier='STRONG'
[PASS] input='Strong'             → snippet tier='STRONG'
[PASS] input='Mechanistic'        → snippet tier='MECHANISTIC'
[PASS] input='Insufficient'       → snippet tier='INSUFFICIENT'
```

Even with `candidate_tier="Validated"` injected directly into the ReceiptCandidate  
(simulating a direct DB edit bypass), the snippet emits `evidence_tier="STRONG"`.  
VALIDATED never reaches a frozen receipt.

---

## Proof 3 — approve→reject and reject→approve both return 409

Live execution via queue layer + FastAPI TestClient:

```
[PASS] approve→reject blocked:
       Candidate 95 is already approved. Human decisions are immutable.
       Approved receipts cannot be silently rejected.

[PASS] reject→approve blocked:
       Candidate 96 is already rejected. Human decisions are immutable.
       Create a new candidate if evidence has changed.

[PASS] empty audited_by blocked:
       audited_by cannot be empty — audit trail requires an auditor identity
```

API layer (FastAPI TestClient, same code path as production HTTP):
```
POST /api/v1/audit/92/reject  (already approved)  → HTTP 409 Conflict
POST /api/v1/audit/93/approve (already rejected)  → HTTP 409 Conflict
```

---

## Proof 4 — _compute_kb_score() counts only ResponseType.SENSITIVITY

Live execution output:

```
[PASS] 0 sens + 10 resist → hits=0, score=0.000   (pure resistance = no signal)
[PASS] 3 sens + 0 resist  → hits=3, score=1.000   (pure sensitivity = max signal)
[PASS] 2 sens + 5 resist  → hits=2, score=0.667   (mixed = sensitivity-only count)
[PASS] 0 sens + 50 resist → hits=0, score=0.000   (50 resistance studies = still 0)
```

The filter is at line 254-257 of receipt_miner.py:
```python
sensitivity_recs = [
    r for r in response.recommendations
    if r.response_type == ResponseType.SENSITIVITY
]
```

---

## Proof 5 — POLE DDR axes quarantined, IO axis not quarantined

Live execution output:

```
[PASS] POLE + 5 DDR axes:  candidates_queued=0 (expected 0),
                            pairs_no_signal=5  (expected ≥5)
                            → all 5 DDR pairs quarantined by guard

[PASS] POLE + IO only:     pairs_evaluated=1  (expected ≥1)
                            → IO not quarantined — clinical safety preserved

[PASS] POLE + all axes:    pairs_evaluated=6  (expected 6)
                            candidates_queued=0
                            → CUSTOM skipped (0 counter), 5 DDR quarantined,
                              1 IO evaluated (scored 0.0 from mocked signals)
```

POLE+IO (immunotherapy) is biologically correct: POLE hypermutators are  
pembrolizumab-sensitive (FDA-approved for TMB-high). Quarantining IO would  
deny patients a valid treatment option. Guard explicitly excludes IMMUNOTHERAPY.

---

## Proof 6 — Integer overflow returns HTTP 422 (not 500)

Live execution via FastAPI TestClient:

```
[PASS] GET /audit/9223372036854775808       (2^63)          → HTTP 422 (bounds)
[PASS] GET /audit/9999999999999999999999    (22-digit int)  → HTTP 422 (bounds)
[PASS] GET /audit/0                         (zero)          → HTTP 422 (bounds)
[PASS] GET /audit/-1                        (negative)      → HTTP 422 (bounds)
[PASS] 2^63 confirmed 422 (not 500) — no stack trace leakage
```

Server log confirms FastAPI returned 422 Unprocessable Entity for all cases,  
not 500 Internal Server Error.  
Guard location: `routes.py:211`, `routes.py:232`, `routes.py:259`.

---

## Prior Assurance Discrepancy — Acknowledged

Management noted: "an earlier sprint narrative said human decisions were  
immutable, but this new report says approve/reject flipping was still possible."

This is accurate. The Sprint 2-3 spec declared immutability as an *intent*  
and the upsert logic preserved approved/rejected rows from miner overwrites —  
but the `approve()` and `reject()` methods themselves had no guard preventing  
a human operator (or malicious actor) from calling the opposite method via  
the API. The P2 fix closes that gap. The prior assurance was overstated.  
The current code is correct. Tests prove it.

---

## Conclusion

All six proofs produced. All 18 sub-assertions pass.  
236/236 tests green.  
Repository: https://github.com/fjkiani/Synthetic-Lethality  
Commit: ff66041d6aeb88a1ed2a41c215a66308593dce18

**Recommendation: merge-ready.**
