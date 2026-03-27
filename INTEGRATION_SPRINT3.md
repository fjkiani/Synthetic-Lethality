# INTEGRATION_SPRINT3.md — Sprint 3: PRISM + Cache Warm-up + Coverage Auto-refresh

**Version**: v4.2.0  
**Baseline**: Sprint 2 (158/158 passing)  
**Test count after**: 180/180 passing  
**New tests added**: 22 (15 sprint-specific + 7 PRISM loader)

---

## What changed

### Addition 1 — PRISM viability signal (`data/prism_loader.py`)

New file. Loads PRISM Repurposing Secondary Screen from DepMap Figshare
and computes drug sensitivity associations per gene × axis.

**Key properties:**
- `prism_stratify(gene, axis)` → `PRISMResult | None` — never raises
- Stratifies by mutant vs WT cell lines (reuses DepMap mutation table)
- Mann-Whitney U on AUC values (`alternative="less"` — mutant more sensitive)
- Parquet cache at `.cache/depmap/prism_secondary.parquet`
- Returns `None` on any failure: no data, insufficient lines, network error

**Bonus wiring in `receipt_miner.py`:**
- PRISM is a BONUS signal only — original weights unchanged
  (`CRISPR 0.35 / GDSC 0.30 / KB 0.25 / Expr 0.10 = 1.0`)
- If `prism_fdr < 0.05`: `confidence += 0.05` (capped at 1.0)
- `prism_delta_auc` stored in `ReceiptCandidate` for transparency
- All 158 Sprint 2 tests pass unchanged

### Addition 2 — GDSC cache warm-up (`api/app.py` lifespan)

On startup, the lifespan handler triggers GDSC1 + GDSC2 downloads
in a thread pool (non-blocking to the event loop).

**Key properties:**
- Runs in `ThreadPoolExecutor` via `asyncio.wait_for` — 60 s timeout
- On success: logs `"GDSC cache warm-up complete: N rows"`
- On timeout (> 60 s): logs warning, continues startup — no crash
- On download failure: logs warning, continues startup — no crash
- First `mine_receipts()` call gets instant cache hit instead of 29 MB download

### Addition 3 — Coverage heat map auto-refresh (`audit/queue.py` + `audit/routes.py`)

`GET /audit/coverage` now returns two new headers:

```
Cache-Control: no-cache
X-Coverage-Version: <ISO-8601 UTC timestamp>
```

The timestamp changes whenever:
- `AuditQueue.approve()` is called
- `AuditQueue.reject()` is called
- A `POST /audit/mine` background task completes

**Doctor portal polling pattern:**
```js
// Poll every 30 s; refresh heat map only when version changes
const prev = localStorage.getItem('coverageVersion');
const res = await fetch('/api/v1/audit/coverage');
const version = res.headers.get('x-coverage-version');
if (version !== prev) {
  localStorage.setItem('coverageVersion', version);
  renderHeatMap(await res.json());
}
```

---

## Files changed

| File | Change |
|------|--------|
| `sl_agent/data/prism_loader.py` | **NEW** — `PRISMResult`, `prism_stratify()`, parquet cache |
| `sl_agent/multimodal/receipt_miner.py` | PRISM bonus wiring: `_compute_prism_score()`, `PRISM_BONUS`, `PRISM_FDR_GATE` |
| `sl_agent/api/app.py` | GDSC warm-up in lifespan handler (thread pool, 60 s timeout) |
| `sl_agent/audit/queue.py` | `_coverage_version` state, `_set_coverage_version()`, `get_coverage_version()`, invalidation in `approve()` / `reject()` |
| `sl_agent/audit/routes.py` | `Response` import, `get_coverage_version` import, headers in `coverage_grid()`, coverage pre-warm in `_run_mine()` |
| `sl_agent/multimodal/tests/test_prism_sprint3.py` | **NEW** — 15 Sprint 3 tests |
| `sl_agent/audit/tests/test_audit_sprint3.py` | **NEW** — 7 coverage/warmup tests |

---

## Drop-in instructions

```bash
# 1. Unzip into your project root (same level as sl_agent/)
unzip sl_agent_sprint3_prism.zip -d /path/to/your/project

# 2. No new pip deps required:
#    scipy, httpx, pandas already installed from Sprint 2

# 3. Verify
cd /path/to/your/project
python -m pytest sl_agent/tests/ sl_agent/kb/tests/ \
                 sl_agent/multimodal/tests/ sl_agent/audit/tests/ -v
# Expected: 180 passed

# 4. Start the server
uvicorn sl_agent.api.app:app --reload --port 8000
# On startup you will see:
#   startup | message=GDSC cache warm-up complete: N rows
#   (or a warning if GDSC unreachable — server still starts)

# 5. Trigger a mine run
curl -X POST http://localhost:8000/api/v1/audit/mine \
  -H "Content-Type: application/json" \
  --get --data-urlencode "gene_panel=BRCA1" \
  --data-urlencode "gene_panel=MBD4" \
  --data-urlencode "axes=parp_inhibitors" \
  --data-urlencode "axes=cytidine_analogs"

# 6. Check coverage version header
curl -I http://localhost:8000/api/v1/audit/coverage
# Look for: X-Coverage-Version: 2026-03-26T...
```

---

## Invariants preserved (do not relax)

| Invariant | Status |
|-----------|--------|
| `_project_tier()` NEVER returns "Validated" | ✅ sabotage test passes (180/180) |
| Base weights CRISPR/GDSC/KB/Expr sum to 1.0 | ✅ unchanged |
| PRISM is additive bonus only (not rebalanced) | ✅ `PRISM_BONUS = 0.05`, capped at 1.0 |
| Startup never crashes on GDSC/PRISM failure | ✅ all exceptions swallowed with log warning |
| All 158 Sprint 2 tests pass | ✅ |
| RS constraints unchanged | ✅ |

---

## Sprint 4 preview (planned, not built)

- PRISM compound metadata integration (pull compound targets for axis refinement)
- Background coverage auto-push via WebSocket (replace polling with push)
- GDSC + PRISM cold-start benchmark: track cache hit rate, surface in `/audit/stats`
