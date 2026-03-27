# RS Patch — Remaining Work Spec

## Critical Bug (Fix First)

`_build_interpretation` in `/home/user/workspace/sl_agent/multimodal/modality_fuser.py` is called as
`_build_interpretation(row, rs_score)` at line 221 but the function signature at line 226 is
`def _build_interpretation(row: EvidenceRow) -> str:` — it does NOT accept `rs_score`.

Fix: change signature to `def _build_interpretation(row: EvidenceRow, rs_score: Optional[RSScore] = None) -> str:`

The function body does NOT need to use rs_score — the tier has already been set before _build_interpretation
is called. The rs_score param is accepted but not used in the interpretation text (tier is already correct).

## File 1: literature_receipts.py
Path: `/home/user/workspace/sl_agent/multimodal/literature_receipts.py`

Add PKMYT1 frozen receipts. Read the existing file first to understand the pattern.
The existing file uses `_RECEIPT_CACHE: Dict[Tuple[str, CandidateAxis], Dict[str, ModalityStatus]]`
pattern for frozen receipts (look at how ATR_WEE1 entries are structured).

Add this entry:
```python
("CCNE1", CandidateAxis.PKMYT1): {
    "crispr":      ModalityStatus.POSITIVE,   # Sanger/DepMap CCNE1-amp → PKMYT1 dependency
    "expression":  ModalityStatus.POSITIVE,
    "in_vitro":    ModalityStatus.POSITIVE,
    "clinical":    ModalityStatus.MISSING,    # RP-6306 Phase I ongoing
},
```

Also add a MBD4+PKMYT1 entry (MBD4 has no direct PKMYT1 evidence — all MISSING):
```python
("MBD4", CandidateAxis.PKMYT1): {
    "crispr":      ModalityStatus.MISSING,
    "expression":  ModalityStatus.MISSING,
    "in_vitro":    ModalityStatus.MISSING,
    "clinical":    ModalityStatus.MISSING,
},
```

## File 2: matrix_builder.py
Path: `/home/user/workspace/sl_agent/multimodal/matrix_builder.py`

### 2a. Add PKMYT1 to `_AXIS_META`
Find the `_AXIS_META` dict (maps CandidateAxis → metadata). Add:
```python
CandidateAxis.PKMYT1: AxisMeta(
    label="PKMYT1",
    description="PKMYT1 kinase inhibition (RP-6306 class); replication stress axis",
    target_genes=["PKMYT1"],
),
```
(Use whatever AxisMeta fields already exist in the file — match the pattern.)

### 2b. Add PKMYT1 to `_axis_to_target_genes()`
In the function `_axis_to_target_genes()` (or equivalent), add:
```python
CandidateAxis.PKMYT1: ["PKMYT1"],
```

### 2c. Add PKMYT1 inhibitors to `pharmacologic_analyzer._DRUG_TO_AXIS`
In `pharmacologic_analyzer.py` at `/home/user/workspace/sl_agent/multimodal/pharmacologic_analyzer.py`,
find `_DRUG_TO_AXIS` dict. Add:
```python
"rp-6306":    CandidateAxis.PKMYT1,
"rp6306":     CandidateAxis.PKMYT1,
"pkmyt1i":    CandidateAxis.PKMYT1,
"pkmyt1 inhibitor": CandidateAxis.PKMYT1,
```

### 2d. Wire RS scoring in `matrix_builder.py`

Find the main `build_matrix()` or equivalent entry-point function that constructs the EvidenceMatrix
and calls `fuse_matrix()`. Add:

1. Import at top of file (if not already there):
```python
from sl_agent.multimodal.replication_stress import compute_rs_score
```

2. Before calling `fuse_matrix(matrix)`, compute rs_score:
```python
rs_score = compute_rs_score(query.rs_features) if query.rs_features else None
```

3. Pass rs_score to fuse_matrix:
```python
fuse_matrix(matrix, rs_score=rs_score)
```

4. Store rs_score on matrix (if EvidenceMatrix has rs_score field — it was added to models.py):
```python
matrix.rs_score = rs_score
```

IMPORTANT: Check if `query` has `rs_features` attribute (it was added to `MultiModalQueryInput`
in models.py). If the function signature uses a different variable name for the query object,
use that.

## File 3: test_rs_benchmark.py (NEW)
Path: `/home/user/workspace/sl_agent/multimodal/tests/test_rs_benchmark.py`

Create a complete pytest test file with these test classes/functions:

### Class TestRSScoring — RS unit tests
```python
def test_mbd4_alone():
    """MBD4 alone: rs=0.05, level='low'"""
    features = RSFeatureSet(mbd4_lof=True)
    result = compute_rs_score(features)
    assert abs(result.score - 0.05) < 1e-9
    assert result.level == "low"

def test_mbd4_tp53():
    """MBD4+TP53: 0.05+0.10=0.15, level='low'"""
    features = RSFeatureSet(mbd4_lof=True, tp53_lof=True)
    result = compute_rs_score(features)
    assert abs(result.score - 0.15) < 1e-9
    assert result.level == "low"

def test_mbd4_tp53_arid1a():
    """MBD4+TP53+ARID1A: 0.05+0.10+0.20=0.35, level='moderate'"""
    features = RSFeatureSet(mbd4_lof=True, tp53_lof=True, arid1a_lof=True)
    result = compute_rs_score(features)
    assert abs(result.score - 0.35) < 1e-9
    assert result.level == "moderate"

def test_ccne1_tp53():
    """CCNE1+TP53: 0.30+0.10=0.40, level='high'"""
    features = RSFeatureSet(ccne1_amplified=True, tp53_lof=True)
    result = compute_rs_score(features)
    assert abs(result.score - 0.40) < 1e-9
    assert result.level == "high"

def test_ccne1_myc_arid1a():
    """CCNE1+MYC+ARID1A: 0.30+0.20+0.20=0.70, level='high'"""
    features = RSFeatureSet(ccne1_amplified=True, myc_amplified=True, arid1a_lof=True)
    result = compute_rs_score(features)
    assert result.score <= 1.0
    assert result.level == "high"

def test_no_features():
    """No features: score=0.0, level='none'"""
    features = RSFeatureSet()
    result = compute_rs_score(features)
    assert result.score == 0.0
    assert result.level == "none"
```

### Class TestRSTierPromotion
Use the fuse_matrix function directly with mock EvidenceRow objects.

```python
def test_ccne1_pkmyt1_promotes_to_strong():
    """CCNE1 with PKMYT1 axis, high RS + CRISPR receipt → Strong"""
    # Build a minimal EvidenceMatrix with one PKMYT1 row for CCNE1
    # CRISPR=POSITIVE, PRISM/GDSC=MISSING → n_pos includes CRISPR
    # RS level = high (CCNE1+TP53)
    # Expected: recommendation_tier == "Strong"
    ...

def test_mbd4_alone_stays_mechanistic():
    """MBD4 alone: rs=low → ATR/WEE1 stays Mechanistic (no promotion)"""
    # ATR/WEE1 row with typical Mechanistic evidence
    # RS = low (mbd4_lof only)
    # Expected: recommendation_tier == "Mechanistic"
    ...

def test_moderate_rs_no_pharma_stays_mechanistic():
    """MBD4+TP53+ARID1A: rs=moderate but no pharma/CRISPR → stays Mechanistic"""
    # ATR/WEE1 row with Mechanistic evidence
    # RS = moderate
    # No CRISPR, no pharma (PRISM/GDSC MISSING)
    # Expected: tier stays "Mechanistic" (RS moderate ≠ sufficient)
    ...

def test_high_rs_no_evidence_cannot_promote():
    """Sabotage: high RS (CCNE1+TP53) but all modalities MISSING → no promotion"""
    # ATR/WEE1 row with ALL modalities MISSING
    # RS = high
    # Expected: NOT "Strong" (RS alone insufficient — need crispr_pos or pharma_pos)
    ...
```

### Class TestRSRegressionGuard
```python
def test_non_rs_axis_unaffected():
    """PARP axis should not be affected by RS score"""
    # Build EvidenceMatrix with PARP axis
    # Provide high RS score
    # Expected: tier for PARP is unchanged vs. without RS
    ...

def test_rs_cannot_reach_validated():
    """RS modulation must never promote directly to Validated"""
    # Even with high RS + CRISPR, tier must be Strong, not Validated
    ...
```

### How to build test rows
Look at how existing `test_multimodal.py` builds EvidenceMatrix objects.
Use whatever fixtures/helpers already exist in the test file. Import from:
```python
from sl_agent.multimodal.replication_stress import RSFeatureSet, RSScore, compute_rs_score
from sl_agent.multimodal.models import (
    EvidenceMatrix, EvidenceRow, ModalityEvidence, ModalityStatus, CandidateAxis
)
from sl_agent.multimodal.modality_fuser import fuse_matrix
```

### IMPORTANT constraints
- RS can NEVER promote to "Validated" — check the promotion code in modality_fuser.py
- Promotion only works: Mechanistic → Strong, when rs_score.level=="high" AND (crispr_pos OR pharma_pos)
- Non-RS-sensitive axes (PARP, POLQ, BASE_EXCISION_REPAIR, etc.) must be unaffected

## Run Tests

After all edits:
```bash
cd /home/user/workspace && python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ -v 2>&1
```

Must be: all existing tests pass + new RS tests pass. Zero failures.

## Package

After tests pass, create zip:
```bash
cd /home/user/workspace
mkdir -p sl_patch_rs/sl_agent/multimodal/tests

cp sl_agent/multimodal/replication_stress.py   sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/models.py               sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/modality_fuser.py       sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/matrix_builder.py       sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/literature_receipts.py  sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/pharmacologic_analyzer.py sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/tests/test_rs_benchmark.py sl_patch_rs/sl_agent/multimodal/tests/

# Write INTEGRATION_RS.md (see below)

zip -r sl_agent_rs_patch.zip sl_patch_rs/
rm -rf sl_patch_rs/
```

## INTEGRATION_RS.md content

```markdown
# SL Agent — RS Axis Patch (v4.0.0 → v4.1.0)

## What's New
- Replication Stress (RS) scoring axis (ATR/WEE1 and PKMYT1)
- RS-aware tier modulation: Mechanistic → Strong when rs_score.level=="high" AND (CRISPR or pharma evidence present)
- PKMYT1 frozen receipts: CCNE1-amp → PKMYT1 dependency (Sanger/DepMap); RP-6306 Phase I ongoing
- RS score is fully transparent — all weights documented in replication_stress.py

## Files Changed
| File | Status | Description |
|------|--------|-------------|
| sl_agent/multimodal/replication_stress.py | NEW | RSFeatureSet, RSScore, compute_rs_score() |
| sl_agent/multimodal/models.py | MODIFIED | +PKMYT1 axis, +rs_features field, +rs_score field |
| sl_agent/multimodal/modality_fuser.py | MODIFIED | RS-aware tier modulation |
| sl_agent/multimodal/matrix_builder.py | MODIFIED | RS wiring, PKMYT1 axis metadata |
| sl_agent/multimodal/literature_receipts.py | MODIFIED | PKMYT1 frozen receipts |
| sl_agent/multimodal/pharmacologic_analyzer.py | MODIFIED | PKMYT1 inhibitor entries |
| sl_agent/multimodal/tests/test_rs_benchmark.py | NEW | RS unit + benchmark + regression tests |

## Installation
1. Overwrite the 7 files listed above (keep all other files unchanged)
2. Run tests: `python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ -v`
3. All tests must pass (existing 84 + new RS tests)

## RS Score Weights
| Feature | Weight |
|---------|--------|
| CCNE1 amplified | 0.30 |
| MYC amplified | 0.20 |
| ARID1A LOF | 0.20 |
| MSI-H | 0.15 |
| TP53 LOF | 0.10 |
| High ploidy | 0.10 |
| MBD4 LOF | 0.05 |

Thresholds: high ≥ 0.40 | moderate ≥ 0.20 | low > 0.00 | none = 0.0

## RS Tier Promotion Rules (Strict)
- RS modulation is ATR/WEE1 and PKMYT1 ONLY
- Promotion: Mechanistic → Strong ONLY when (rs_score.level == "high") AND (CRISPR+ OR pharma+)
- RS CANNOT promote to Validated
- MBD4 alone → rs=0.05 (low) → no promotion
- All other axes: unaffected by RS score
```
