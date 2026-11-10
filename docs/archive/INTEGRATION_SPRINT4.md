# INTEGRATION_SPRINT4.md — Atlas Panel v1 · Sprint 4

**Version**: v4.3.0  
**Sprint**: 4 — Promotion Pipeline · Coverage Heat Map · Pan-cancer Gene Panel · Stats Upgrade  
**Test count**: 195/195 passing (0 failures)  
**Date**: 2026-03-26

---

## What shipped

### 1  Promotion Pipeline — `generate_frozen_receipt_snippet()`

`AuditQueue.approve()` now returns a fully-populated `dict` instead of `None`:

```python
{
    "candidate": ReceiptCandidate,          # updated object with status="approved"
    "frozen_receipt_snippet": str           # machine-readable receipt string
}
```

The snippet is produced by `generate_frozen_receipt_snippet(candidate)` in `audit/queue.py` and looks like:

```
FROZEN_RECEIPT | gene=MBD4 | axis=cytidine_analogs | tier=STRONG | confidence=0.82
  signals: crispr=0.35 gdsc=0.30 kb=0.25 expr=0.10 prism=0.05
  audited_by=dr.smith | generated_at=2026-03-26T20:00:00Z
```

**Key invariant**: `tier` is derived from `candidate.candidate_tier`. It can NEVER be `VALIDATED` — the sabotage guard in `_project_tier()` remains intact.

**Endpoint change** — `POST /api/v1/audit/{id}/approve`:

```json
{
  "ruo": true,
  "disclaimer": "Research use only. Not validated for clinical decision-making.",
  "data": {
    "status": "approved",
    "candidate_id": 42,
    "candidate": { "...ReceiptCandidate fields..." },
    "frozen_receipt_snippet": "FROZEN_RECEIPT | gene=MBD4 | ..."
  }
}
```

---

### 2  Coverage Heat Map — Full Matrix + Summary Block

`GET /api/v1/audit/coverage` now returns a complete matrix dict (no longer just an array wrapped by `_wrap()`):

```json
{
  "ruo": true,
  "disclaimer": "...",
  "coverage": [
    {
      "gene": "MBD4",
      "axis": "cytidine_analogs",
      "tier": "Validated",
      "confidence": 0.95,
      "promoted_to_frozen": true,
      "audited_by": "dr.smith"
    },
    {
      "gene": "ATM",
      "axis": "parp_inhibitors",
      "tier": "not_covered",
      "confidence": null,
      "promoted_to_frozen": false,
      "audited_by": null
    }
  ],
  "summary": {
    "validated": 1,
    "strong": 3,
    "mechanistic": 12,
    "not_covered": 8,
    "total_pairs": 24
  },
  "coverage_version": "a1b2c3d4"
}
```

- `coverage` contains one entry per (gene × axis) combination across all panel genes and all axes
- Pairs with no queue entry are represented with `tier="not_covered"`, `confidence=null`
- `summary` aggregates all tier counts — `validated + strong + mechanistic + not_covered == total_pairs`
- `coverage_version`, `X-Coverage-Version` header, and `Cache-Control: no-cache` are preserved from Sprint 3

**Doctor portal heat map wiring** (UI teams): poll `GET /api/v1/audit/coverage`, diff `coverage_version` to detect staleness, render `summary` counts for the overview chips, render `coverage[]` as a gene × axis grid with tier-colour mapping.

---

### 3  Pan-cancer Gene Panel — `ddr_panel_v1.json` + `load_default_panel()`

**New file**: `sl_agent/data/gene_panels/ddr_panel_v1.json`

Six gene groups, 36 unique genes (deduplicated):

| Group | Key genes |
|-------|-----------|
| HRR | BRCA1, BRCA2, PALB2, RAD51, RAD51C, RAD51D, BRIP1 |
| FA | FANCA, FANCD2 |
| ATR/CHK | ATM, ATR, CHEK1, CHEK2, ATRIP |
| RS / Replication | CCNE1, CDK2, WEE1, CDC25A, PLK1, PKMYT1, RRM2, PCNA, RFC1, MCM2 |
| BER / MMR | MBD4, TET2, DNMT3A, PARP1, MUTYH, MLH1, MSH2, MSH6, PMS2 |
| Chromatin | ARID1A, SMARCA4, SLFN11, RNF144A |

**New function** in `sl_agent/multimodal/receipt_miner.py`:

```python
from sl_agent.multimodal.receipt_miner import load_default_panel

genes = load_default_panel()  # → List[str], 36 genes, always deduplicated
```

- Falls back to the legacy `DEFAULT_GENE_PANEL` constant if the JSON file is missing (never raises)
- `POST /api/v1/audit/mine` now uses `load_default_panel()` when no `gene_panel` query param is provided

**Adding genes to the panel**: edit `ddr_panel_v1.json` — no code changes needed. The loader re-reads the file on every `mine` run.

---

### 4  Stats Upgrade — Nested Structure

`GET /api/v1/audit/stats` now returns a fully structured body:

```json
{
  "ruo": true,
  "disclaimer": "...",
  "data": {
    "queue_counts": {
      "pending": 14,
      "approved": 6,
      "rejected": 2,
      "promoted": 3,
      "high_confidence_pending": 5
    },
    "coverage_summary": {
      "validated": 1,
      "strong": 3,
      "mechanistic": 12,
      "not_covered": 8,
      "total_pairs": 24
    },
    "top_pending": [
      {
        "id": 7,
        "gene": "ATR",
        "axis": "atr_wee1",
        "confidence": 0.91,
        "candidate_tier": "Strong"
      }
    ],
    "last_mine_at": "2026-03-26T20:00:00Z",
    "coverage_version": "a1b2c3d4"
  }
}
```

- `queue_counts.high_confidence_pending` = pending candidates with `confidence ≥ 0.70`
- `top_pending` = up to 5 pending candidates sorted by `confidence DESC`
- `last_mine_at` = ISO-8601 timestamp of the last completed `POST /audit/mine` (null if never run)
- `coverage_summary` mirrors the summary block from `GET /audit/coverage` — use either endpoint for the same numbers

**Doctor portal dashboard wiring**: `queue_counts` feeds the Kanban column badges; `top_pending` feeds the "needs review" priority list; `coverage_summary` feeds the heat-map overview; `last_mine_at` feeds the "last updated" chip.

---

## Files changed

| File | Type | Change |
|------|------|--------|
| `sl_agent/data/gene_panels/ddr_panel_v1.json` | NEW | DDR panel v1 — 6 groups, 36 genes |
| `sl_agent/multimodal/receipt_miner.py` | MODIFIED | `load_default_panel()` + `_PANEL_JSON` path |
| `sl_agent/audit/queue.py` | MODIFIED | `generate_frozen_receipt_snippet()`, `set_last_mine_at()`, `get_last_mine_at()`, upgraded `approve()` → dict, upgraded `stats()` → nested, upgraded `coverage()` → matrix+summary |
| `sl_agent/audit/routes.py` | MODIFIED | approve returns snippet, coverage returns matrix dict, mine uses `load_default_panel()`, `set_last_mine_at()` wired |
| `sl_agent/audit/tests/test_audit.py` | MODIFIED | 1 test updated: `test_post_mine_uses_default_panel_when_no_params` uses `load_default_panel()` |
| `sl_agent/audit/tests/test_audit_sprint4.py` | NEW | 15 tests: FrozenReceiptSnippet(5), CoverageMatrix(4), GenePanel(3), StatsUpgrade(3) |
| `sl_agent/tests/test_api.py` | MODIFIED | `client` fixture patches GDSC loaders to prevent 60-s startup hang in CI |

---

## Test results

```
195 passed in 11.95s
  sl_agent/tests/           19 tests  ✅
  sl_agent/kb/tests/        28 tests  ✅
  sl_agent/multimodal/tests/ 98 tests  ✅
  sl_agent/audit/tests/     50 tests  ✅
```

Breakdown of Sprint 4 new/modified tests:

```
test_audit_sprint4.py   15 new tests
  TestFrozenReceiptSnippet (5)
    ✅ approve_returns_dict_with_snippet
    ✅ snippet_contains_correct_gene_and_axis
    ✅ snippet_contains_numeric_scores
    ✅ snippet_does_not_contain_validated_tier
    ✅ api_approve_returns_snippet_in_response
  TestCoverageMatrix (4)
    ✅ coverage_returns_full_matrix_shape
    ✅ coverage_summary_counts_add_up
    ✅ coverage_includes_not_covered_entries
    ✅ coverage_mbd4_has_validated_tier
  TestGenePanel (3)
    ✅ load_default_panel_returns_deduped_list
    ✅ load_default_panel_contains_expected_genes
    ✅ mine_uses_default_panel_when_no_params
  TestStatsUpgrade (3)
    ✅ stats_returns_all_required_keys
    ✅ stats_top_pending_sorted_by_confidence
    ✅ stats_coverage_summary_matches_coverage_endpoint

test_audit.py (modified)
    ✅ test_post_mine_uses_default_panel_when_no_params  (updated: DEFAULT_GENE_PANEL → load_default_panel())
```

---

## Running the suite locally

```bash
cd /path/to/sl_agent_project
pip install -e ".[dev]"
python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ sl_agent/audit/tests/ -v
# Expected: 195 passed
```

```bash
# Start the server
uvicorn sl_agent.api.app:app --reload --port 8000

# Smoke-test Sprint 4 endpoints
curl -s http://localhost:8000/api/v1/audit/coverage | python -m json.tool
curl -s http://localhost:8000/api/v1/audit/stats   | python -m json.tool

# Approve a candidate (requires a seeded candidate with id=1)
curl -s -X POST http://localhost:8000/api/v1/audit/1/approve \
  -H "Content-Type: application/json" \
  -d '{"audited_by": "dr.test"}' | python -m json.tool
# → body.data.frozen_receipt_snippet should be present
```

---

## Key invariants preserved

| Invariant | Status |
|-----------|--------|
| `_project_tier()` never returns "Validated" | ✅ sabotage test passes |
| `frozen_receipt_snippet` tier ≠ VALIDATED | ✅ `test_snippet_does_not_contain_validated_tier` |
| Base weights CRISPR/GDSC/KB/Expr unchanged | ✅ `test_weights_still_sum_to_one` |
| PRISM bonus ≤ 0.05, never changes base weights | ✅ `test_prism_bonus_capped_at_one` |
| RS modulation: ATR/WEE1/PKMYT1 only | ✅ `test_non_rs_axis_unaffected` |
| RS cannot reach Validated tier | ✅ `test_rs_cannot_reach_validated` |
| All 69 RS benchmark tests pass | ✅ |
| Startup never crashes on GDSC/DepMap failure | ✅ |
| All endpoints RUO-tagged | ✅ |

---

## What comes next (Sprint 5 candidates)

- Doctor portal UI — wire coverage heat map, approval workflow, stats dashboard
- GDSC v9 / DepMap 25Q1 refresh triggers
- Receipt promotion to frozen KB (persist approved `frozen_receipt_snippet` as immutable record)
- Batch approval workflow (bulk approve / bulk reject by tier threshold)
- Coverage alert webhook (notify when coverage_version changes)
