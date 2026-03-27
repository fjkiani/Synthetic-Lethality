# RS Phase 2 Spec — Combinatorial Vulnerability Framework

## What's Already Done (DO NOT TOUCH)
- `replication_stress.py` — RSFeatureSet, RSScore, compute_rs_score(), rs_is_high(), etc. COMPLETE.
- `models.py` — PKMYT1, rs_features, rs_score all already added. COMPLETE.
- `modality_fuser.py` — RS-aware tier modulation COMPLETE and tested.
- `matrix_builder.py` — RS wiring COMPLETE (Steps 6-7 already present).
- `literature_receipts.py` — PKMYT1 receipts COMPLETE.
- `tests/test_rs_benchmark.py` — 12 tests, all passing. DO NOT BREAK THEM.

Current test count: 96 tests passing. All must continue to pass after Phase 2.

---

## Step 1 — Add CombinationVulnerability to models.py

File: `/home/user/workspace/sl_agent/multimodal/models.py`

Add AFTER the `EvidenceRow` class and BEFORE the `EvidenceMatrix` class:

```python
# ── Combinatorial vulnerability ───────────────────────────────────────────────

class CombinationVulnerability(BaseModel):
    """
    Two-pathway convergence creating emergent SL that neither pathway alone produces.
    
    Examples:
      BER deficiency (MBD4-LOF) × checkpoint inhibition (ATR/WEE1):
        Neither pathway alone reaches Strong, but combined they create mitotic catastrophe.
      Chromatin remodeling loss (ARID1A-LOF) × ATR checkpoint:
        REALIZED — CRISPR data present, combined tier = Strong.
    """
    pathway_a: str                     # e.g. "BER_deficiency"
    pathway_a_gene: str                # e.g. "MBD4"
    pathway_b: str                     # e.g. "checkpoint_inhibition"
    pathway_b_targets: List[str]       # e.g. ["ATR", "CHK1", "WEE1"]
    convergence_mechanism: str         # e.g. "mitotic_catastrophe"

    pathway_a_alone_tier: str          # What pathway A gets alone (e.g. "Mechanistic candidate only")
    pathway_b_alone_tier: str          # What pathway B gets alone (e.g. "Mechanistic candidate only")
    combined_tier: Optional[str] = None  # Projected combined tier (None = insufficient data)

    data_gaps: List[str] = Field(default_factory=list)
    # Exact data needed to unlock combinatorial tier. Empty = tier already realized.

    rs_score: Optional[float] = None   # RS score for this combination
    rs_level: Optional[str] = None     # RS level

    rationale: str                     # One-paragraph mechanistic chain
```

Also add `combination_vulnerabilities` field to `EvidenceMatrix`:

Find the `EvidenceMatrix` class and add to its fields (after `rs_score`):
```python
    combination_vulnerabilities: List["CombinationVulnerability"] = Field(
        default_factory=list,
        description="Combinatorial two-pathway vulnerabilities detected for this gene/cancer context",
    )
```

---

## Step 2 — Add assess_combinatorial_vulnerability() to replication_stress.py

File: `/home/user/workspace/sl_agent/multimodal/replication_stress.py`

Add AFTER the existing `rs_is_actionable()` function.

IMPORTANT: CombinationVulnerability is defined in models.py. Use TYPE_CHECKING to avoid circular imports:

```python
from __future__ import annotations
# (already at top of file)
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .models import EvidenceMatrix, CombinationVulnerability
```

Then add:

```python
# ── Combinatorial vulnerability assessment ────────────────────────────────────

# Genes known to drive BER/DDR deficiency relevant to RS
_BER_DDR_GENES = {"MBD4", "MLH1", "MSH2", "MSH6", "BRCA1", "BRCA2"}

# Checkpoint axes that can form combinatorial vulnerability
_CHECKPOINT_AXES_FOR_COMBINATION = {
    "atr_wee1": ["ATR", "CHK1", "WEE1"],
    "pkmyt1":   ["PKMYT1"],
}

def assess_combinatorial_vulnerability(
    gene: str,
    rs_score: "RSScore",
    evidence_matrix: "EvidenceMatrix",
) -> Optional["CombinationVulnerability"]:
    """
    Check if a gene's BER/DDR deficiency creates a combinatorial vulnerability
    with checkpoint kinases, where neither pathway alone is SL.

    Logic:
      1. Gene must be in known BER/DDR deficiency set (or rs_score > 0)
      2. ATR/WEE1 axis must be "Mechanistic candidate only" individually (not already Strong/Validated)
      3. There must be a mechanistic chain: gene_LOF → fork_stress → checkpoint_dependency
      4. Identify exact data gaps required to unlock promotion

    Returns CombinationVulnerability or None if not applicable.

    CONSTRAINT: This function ONLY identifies the vulnerability — it does NOT modify tiers.
    Tier promotion still requires independent CRISPR or pharma evidence via fuse_matrix().
    """
    # Late import to avoid circular dependency
    from .models import CombinationVulnerability

    gene_upper = gene.upper()

    # Must have some RS contribution
    if rs_score.score == 0.0:
        return None

    # Find ATR/WEE1 row in matrix
    atr_row = None
    for row in evidence_matrix.rows:
        if row.axis.value == "atr_wee1":
            atr_row = row
            break

    if atr_row is None:
        return None

    # Both pathways must be below Strong for a "combinatorial" story to make sense
    atr_tier = atr_row.recommendation_tier or "Unknown"
    if atr_tier in ("Validated SL therapeutic lever", "Strong candidate dependency axis"):
        # Already realized — still worth reporting, but as REALIZED combinatorial
        combined_tier = atr_tier
        data_gaps: List[str] = []
        convergence = "replication_fork_collapse"
    else:
        # Not yet realized — identify gaps
        combined_tier = None
        data_gaps = _identify_data_gaps(gene_upper, atr_row)
        convergence = "mitotic_catastrophe"

    # Build mechanistic rationale
    if gene_upper == "MBD4":
        rationale = (
            f"MBD4 LOF → BER collapse → 5mC>T lesions → replication fork stalling → "
            f"ATR/CHK1 checkpoint dependency → if checkpoint inhibited → mitotic catastrophe. "
            f"RS score = {rs_score.score:.2f} ({rs_score.level}). "
            f"{'Data gaps remain before combinatorial tier can be promoted.' if data_gaps else 'Combinatorial tier realized by available evidence.'}"
        )
    elif gene_upper == "ARID1A":
        rationale = (
            f"ARID1A LOF → chromatin access failure → fork instability → "
            f"ATR dependency (Williamson 2016). RS score = {rs_score.score:.2f} ({rs_score.level}). "
            f"{'Combinatorial tier unlocked by available CRISPR + functional data.' if not data_gaps else 'Additional evidence needed.'}"
        )
    else:
        rationale = (
            f"{gene_upper} LOF contributes to replication stress (RS={rs_score.score:.2f}, {rs_score.level}), "
            f"creating potential checkpoint dependency. "
            f"{'Data gaps identified for combinatorial tier promotion.' if data_gaps else 'Evidence sufficient for combinatorial assessment.'}"
        )

    return CombinationVulnerability(
        pathway_a=_gene_to_pathway(gene_upper),
        pathway_a_gene=gene_upper,
        pathway_b="checkpoint_inhibition",
        pathway_b_targets=["ATR", "CHK1", "WEE1"],
        convergence_mechanism=convergence,
        pathway_a_alone_tier=_get_gene_alone_tier(gene_upper),
        pathway_b_alone_tier=atr_tier,
        combined_tier=combined_tier,
        data_gaps=data_gaps,
        rs_score=rs_score.score,
        rs_level=rs_score.level,
        rationale=rationale,
    )


def _gene_to_pathway(gene: str) -> str:
    """Map gene to its deficiency pathway label."""
    pathway_map = {
        "MBD4":  "BER_deficiency",
        "ARID1A": "chromatin_remodeling_loss",
        "BRCA1": "HR_deficiency",
        "BRCA2": "HR_deficiency",
        "MLH1":  "MMR_deficiency",
        "MSH2":  "MMR_deficiency",
        "MSH6":  "MMR_deficiency",
        "CCNE1": "replication_stress_amplification",
    }
    return pathway_map.get(gene, f"{gene}_deficiency")


def _get_gene_alone_tier(gene: str) -> str:
    """Conservative standalone tier for pathway A gene (without combination)."""
    # These are the tier the gene's OWN axis gets individually
    # For BER/chromatin genes, standalone checkpoint dependency is typically Mechanistic
    return "Mechanistic candidate only"


def _identify_data_gaps(gene: str, atr_row) -> List[str]:
    """Identify what data is missing to unlock the combinatorial tier."""
    from .models import ModalityStatus
    gaps = []
    cells = atr_row.cells()

    if cells["crispr"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"ATR/WEE1 CRISPR dependency in {gene}-LOF lines "
            f"(DepMap stratification: Chronos Δdep in {gene}-LOF vs WT)"
        )
    if cells["prism"].status != ModalityStatus.POSITIVE and cells["gdsc"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"ATR/WEE1 drug sensitivity stratified by {gene} status "
            f"(GDSC/PRISM: ceralasertib / adavosertib / berzosertib AUC in {gene}-LOF vs WT)"
        )
    if cells["in_vitro"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"Isogenic {gene}-KO + ATRi synergy experiment "
            f"(Bliss/Loewe combination index in {gene}-KO vs WT)"
        )
    if not gaps:
        # All present — tier should have been promoted already
        gaps = []  # REALIZED

    return gaps
```

---

## Step 3 — Add get_checkpoint_dependency_in_context() to depmap_loader.py

File: `/home/user/workspace/sl_agent/data/depmap_loader.py`

Add AFTER the existing `get_mutant_wt_lines()` function.

First, add these imports at the top of the file (after existing imports):
```python
from typing import Dict  # already imported via other typing imports — check first
# Also add if not present:
from scipy import stats as _scipy_stats  # for t-test; use try/except if scipy not installed
```

Then add:

```python
# ── Stratified Result model ───────────────────────────────────────────────────

class StratifiedResult:
    """
    Result of stratified CRISPR dependency analysis for one checkpoint gene.
    Plain dataclass (not Pydantic) to avoid circular import with multimodal models.
    """
    __slots__ = [
        "checkpoint_gene", "context_gene",
        "mean_dep_mutant", "mean_dep_wt",
        "delta_dep", "p_value", "effect_size_cohend",
        "n_mutant", "n_wt",
        "is_more_essential_in_mutant",
        "notes",
    ]

    def __init__(
        self,
        checkpoint_gene: str,
        context_gene: str,
        mean_dep_mutant: float,
        mean_dep_wt: float,
        delta_dep: float,
        p_value: float,
        effect_size_cohend: float,
        n_mutant: int,
        n_wt: int,
        is_more_essential_in_mutant: bool,
        notes: str = "",
    ):
        self.checkpoint_gene = checkpoint_gene
        self.context_gene = context_gene
        self.mean_dep_mutant = mean_dep_mutant
        self.mean_dep_wt = mean_dep_wt
        self.delta_dep = delta_dep
        self.p_value = p_value
        self.effect_size_cohend = effect_size_cohend
        self.n_mutant = n_mutant
        self.n_wt = n_wt
        self.is_more_essential_in_mutant = is_more_essential_in_mutant
        self.notes = notes

    def __repr__(self) -> str:
        return (
            f"StratifiedResult({self.checkpoint_gene} in {self.context_gene}-LOF: "
            f"Δdep={self.delta_dep:+.3f}, p={self.p_value:.4f}, "
            f"d={self.effect_size_cohend:.3f}, n_mut={self.n_mutant})"
        )


def get_checkpoint_dependency_in_context(
    context_gene: str,
    checkpoint_genes: List[str],
    mutation_type: str = "loss_of_function",
    cancer_type: Optional[str] = None,
) -> Dict[str, "StratifiedResult"]:
    """
    For each checkpoint gene, compute Chronos dependency stratified by context_gene status.

    This answers: "Are {context_gene}-LOF cells more dependent on each checkpoint gene?"

    Args:
        context_gene:     Stratification gene (e.g. "MBD4")
        checkpoint_genes: List of checkpoint/target genes (e.g. ["ATR", "CHEK1", "WEE1", "PKMYT1"])
        mutation_type:    How to define context_gene-LOF (see get_mutant_wt_lines)
        cancer_type:      Optional cancer type filter

    Returns:
        Dict mapping checkpoint_gene → StratifiedResult

    Key threshold for combinatorial tier unlock:
        delta_dep < -0.15 in context_gene-LOF vs WT → combinatorial tier unlockable
        (Chronos: more negative = more essential)

    Note:
        This function DOWNLOADS DepMap data (CRISPR + mutations + sample_info).
        Use the cached versions — first call may be slow, subsequent calls are fast.
        If DepMap data is unavailable, returns empty dict with logged warning.
    """
    results: Dict[str, StratifiedResult] = {}

    try:
        crispr_df = load_crispr_gene_effect()
        mut_df = load_mutations()
        sample_df = load_sample_info()
    except Exception as e:
        logger.warning("get_checkpoint_dependency_in_context: data load failed: %s", e)
        return results

    try:
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=context_gene,
            mutation_df=mut_df,
            sample_info=sample_df,
            cancer_type=cancer_type,
            mutation_type=mutation_type,
        )
    except Exception as e:
        logger.warning("get_checkpoint_dependency_in_context: stratification failed: %s", e)
        return results

    if len(mutant_ids) < 3 or len(wt_ids) < 3:
        logger.warning(
            "get_checkpoint_dependency_in_context: insufficient lines "
            "(n_mut=%d, n_wt=%d) for %s",
            len(mutant_ids), len(wt_ids), context_gene,
        )
        return results

    # Filter CRISPR df to lines we have
    mut_lines_in_crispr = [m for m in mutant_ids if m in crispr_df.index]
    wt_lines_in_crispr  = [w for w in wt_ids     if w in crispr_df.index]

    if len(mut_lines_in_crispr) < 3 or len(wt_lines_in_crispr) < 3:
        logger.warning(
            "get_checkpoint_dependency_in_context: insufficient CRISPR-matched lines "
            "(n_mut=%d, n_wt=%d)",
            len(mut_lines_in_crispr), len(wt_lines_in_crispr),
        )
        return results

    for chk_gene in checkpoint_genes:
        chk_gene_upper = chk_gene.upper()
        if chk_gene_upper not in crispr_df.columns:
            logger.debug("Checkpoint gene %s not in CRISPR matrix", chk_gene_upper)
            continue

        mut_deps = crispr_df.loc[mut_lines_in_crispr, chk_gene_upper].dropna().values
        wt_deps  = crispr_df.loc[wt_lines_in_crispr,  chk_gene_upper].dropna().values

        if len(mut_deps) < 3 or len(wt_deps) < 3:
            continue

        mean_mut = float(mut_deps.mean())
        mean_wt  = float(wt_deps.mean())
        delta    = mean_mut - mean_wt  # negative = more essential in mutant

        # Two-sample t-test (Welch's)
        try:
            import scipy.stats as stats
            t_stat, p_val = stats.ttest_ind(mut_deps, wt_deps, equal_var=False)
        except ImportError:
            # Fallback: no scipy — p_value undefined
            p_val = float("nan")

        # Cohen's d
        pooled_std = float(
            (
                ((len(mut_deps) - 1) * mut_deps.std() ** 2 +
                 (len(wt_deps)  - 1) * wt_deps.std()  ** 2)
                / (len(mut_deps) + len(wt_deps) - 2)
            ) ** 0.5
        )
        cohend = (mean_mut - mean_wt) / pooled_std if pooled_std > 0 else 0.0

        results[chk_gene_upper] = StratifiedResult(
            checkpoint_gene=chk_gene_upper,
            context_gene=context_gene.upper(),
            mean_dep_mutant=round(mean_mut, 4),
            mean_dep_wt=round(mean_wt, 4),
            delta_dep=round(delta, 4),
            p_value=round(float(p_val), 6) if not _is_nan(p_val) else float("nan"),
            effect_size_cohend=round(cohend, 3),
            n_mutant=len(mut_deps),
            n_wt=len(wt_deps),
            is_more_essential_in_mutant=delta < -0.15,
            notes=(
                f"Δdep threshold for combinatorial unlock: < -0.15 "
                f"(current: {delta:+.3f})"
            ),
        )

    return results


def _is_nan(v) -> bool:
    """Safe NaN check that works for both float and non-float."""
    try:
        import math
        return math.isnan(v)
    except (TypeError, ValueError):
        return False
```

IMPORTANT: After adding the function, check if `List` and `Dict` are already imported from `typing` in `depmap_loader.py`. They are — `from typing import Optional, Tuple` is there. Add `List, Dict` to that import if not present.

---

## Step 4 — Wire combinatorial assessment into matrix_builder.py

File: `/home/user/workspace/sl_agent/multimodal/matrix_builder.py`

### 4a. Add import at top
After the existing `from .replication_stress import compute_rs_score` line, add:
```python
from .replication_stress import compute_rs_score, assess_combinatorial_vulnerability
```

### 4b. Add combinatorial assessment after RS wiring (Step 6/7 which already exist)

Find in `build_evidence_matrix()`:
```python
    # ── Step 7: Fuse — compute levels, tiers, agreements ─────────────────────
    matrix = fuse_matrix(matrix, rs_score=rs_score)
    matrix.rs_score = rs_score
```

Add AFTER these two lines:
```python
    # ── Step 8: Assess combinatorial vulnerabilities ──────────────────────────
    if rs_score is not None and rs_score.score > 0.0:
        combo = assess_combinatorial_vulnerability(
            gene=gene,
            rs_score=rs_score,
            evidence_matrix=matrix,
        )
        if combo is not None:
            matrix.combination_vulnerabilities = [combo]
```

---

## Step 5 — Add combinatorial fixture tests to test_rs_benchmark.py

File: `/home/user/workspace/sl_agent/multimodal/tests/test_rs_benchmark.py`

Add a new class at the END of the file (after `TestRSRegressionGuard`):

```python
# ─────────────────────────────────────────────────────────────────────────────
# Class 4: Combinatorial Vulnerability Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestCombinatorialVulnerability:

    def test_mbd4_combinatorial_identified_no_crispr(self):
        """
        Fixture 12: MBD4 BER × checkpoint — no ATR CRISPR data.
        assess_combinatorial_vulnerability returns a CombinationVulnerability
        with combined_tier=None and non-empty data_gaps.
        """
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability
        from sl_agent.multimodal.models import CombinationVulnerability

        # Build matrix with ATR/WEE1 row, CRISPR MISSING
        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.MISSING,
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="MBD4", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet(mbd4_lof=True, tp53_lof=True))
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("MBD4", rs, fused)

        assert combo is not None, "Should identify combinatorial vulnerability for MBD4"
        assert isinstance(combo, CombinationVulnerability)
        assert combo.pathway_a_gene == "MBD4"
        assert combo.pathway_b == "checkpoint_inhibition"
        assert "ATR" in combo.pathway_b_targets
        assert combo.combined_tier is None, "No combined tier without CRISPR data"
        assert len(combo.data_gaps) > 0, "Should have data gaps"
        assert combo.rs_score is not None
        assert combo.convergence_mechanism == "mitotic_catastrophe"

    def test_arid1a_combinatorial_realized(self):
        """
        Fixture 13: ARID1A combinatorial — has CRISPR POSITIVE → combined_tier realized.
        assess_combinatorial_vulnerability returns non-None combined_tier.
        """
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability
        from sl_agent.multimodal.models import CombinationVulnerability

        # ARID1A has ATR CRISPR POSITIVE
        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            crispr_status=ModalityStatus.POSITIVE,
            expression_status=ModalityStatus.POSITIVE,
            in_vitro_status=ModalityStatus.POSITIVE,
            in_vivo_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="ARID1A", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet(arid1a_lof=True, tp53_lof=True))
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("ARID1A", rs, fused)

        assert combo is not None
        assert combo.pathway_a_gene == "ARID1A"
        assert combo.combined_tier is not None, "CRISPR data present → tier should be realized"
        assert len(combo.data_gaps) == 0, "No gaps when all evidence present"

    def test_combinatorial_stored_on_matrix(self):
        """
        build_evidence_matrix wires combinatorial assessment onto matrix.combination_vulnerabilities
        when rs_features are provided and ATR/WEE1 axis is present.
        """
        from sl_agent.multimodal.matrix_builder import build_evidence_matrix
        from sl_agent.multimodal.models import MultiModalQueryInput, CandidateAxis
        from sl_agent.multimodal.replication_stress import RSFeatureSet

        query = MultiModalQueryInput(
            gene="MBD4",
            mutation_type="loss_of_function",
            axes=[CandidateAxis.ATR_WEE1, CandidateAxis.CYTIDINE_ANALOGS],
            include_pharmacologic_stratification=False,
            rs_features=RSFeatureSet(mbd4_lof=True, tp53_lof=True),
        )
        matrix = build_evidence_matrix(query)

        # RS score should be stored
        assert matrix.rs_score is not None
        assert matrix.rs_score.level in ("low", "moderate", "high")

        # Combinatorial vulnerabilities should be assessed
        assert isinstance(matrix.combination_vulnerabilities, list)
        # With MBD4+TP53 RS (low), should still identify the combinatorial structure
        assert len(matrix.combination_vulnerabilities) >= 1

    def test_zero_rs_no_combinatorial(self):
        """Gene with zero RS features should not produce combinatorial assessment."""
        from sl_agent.multimodal.replication_stress import assess_combinatorial_vulnerability

        atr_row = _make_row(
            axis=CandidateAxis.ATR_WEE1,
            label="ATR/WEE1",
            expression_status=ModalityStatus.POSITIVE,
        )
        matrix = EvidenceMatrix(query_gene="UNKNOWNGENE", rows=[atr_row])
        rs = compute_rs_score(RSFeatureSet())  # No features → score = 0.0
        fused = fuse_matrix(matrix, rs_score=rs)

        combo = assess_combinatorial_vulnerability("UNKNOWNGENE", rs, fused)
        assert combo is None, "Zero RS score → no combinatorial vulnerability"
```

---

## Step 6 — Run tests

```bash
cd /home/user/workspace && python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ -v 2>&1
```

All 96 existing tests + new combinatorial tests must pass. Zero failures.

If there are circular import issues:
- `replication_stress.py` must use `TYPE_CHECKING` guard for any models.py imports
- `matrix_builder.py` already imports from `replication_stress.py` — check no new cycles

If `assess_combinatorial_vulnerability` causes circular import:
- Use `from __future__ import annotations` (already present)
- Use `if TYPE_CHECKING:` for the type hints
- Use late `from .models import CombinationVulnerability` inside the function body

---

## Step 7 — Update INTEGRATION_RS.md and package zip

After tests pass, regenerate the patch zip:

```bash
cd /home/user/workspace
rm -rf sl_patch_rs/
mkdir -p sl_patch_rs/sl_agent/multimodal/tests
mkdir -p sl_patch_rs/sl_agent/data

cp sl_agent/multimodal/replication_stress.py       sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/models.py                   sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/modality_fuser.py           sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/matrix_builder.py           sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/literature_receipts.py      sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/pharmacologic_analyzer.py   sl_patch_rs/sl_agent/multimodal/
cp sl_agent/multimodal/tests/test_rs_benchmark.py  sl_patch_rs/sl_agent/multimodal/tests/
cp sl_agent/data/depmap_loader.py                  sl_patch_rs/sl_agent/data/
```

Write `sl_patch_rs/INTEGRATION_RS.md` with this content:

```markdown
# SL Agent — RS Axis Patch v4.1.0 (Phase 1 + Phase 2)

## What's New in Phase 2 (Combinatorial Vulnerability Framework)
- `CombinationVulnerability` model — two-pathway emergent SL that neither pathway alone produces
- `assess_combinatorial_vulnerability()` — detects BER×checkpoint combinatorial vulnerabilities
- `get_checkpoint_dependency_in_context()` in depmap_loader — DepMap stratification engine
- Combinatorial assessment wired into `build_evidence_matrix()` — auto-runs when rs_features provided
- `EvidenceMatrix.combination_vulnerabilities` — list of detected combinatorial vulnerabilities
- 4 new combinatorial tests (fixtures 12–15)

## What Was in Phase 1 (RS Scoring + PKMYT1)
- RS scoring engine (RSFeatureSet, RSScore, compute_rs_score)
- RS-aware tier modulation (Mechanistic → Strong for ATR/WEE1 + PKMYT1 when high RS + evidence)
- PKMYT1 frozen receipts (CCNE1-amp, Sanger/DepMap)
- 12 RS benchmark tests

## Files Changed
| File | Status | Description |
|------|--------|-------------|
| sl_agent/multimodal/replication_stress.py | MODIFIED | +assess_combinatorial_vulnerability() |
| sl_agent/multimodal/models.py | MODIFIED | +CombinationVulnerability, +combination_vulnerabilities on EvidenceMatrix |
| sl_agent/multimodal/modality_fuser.py | PHASE 1 ONLY | No Phase 2 changes |
| sl_agent/multimodal/matrix_builder.py | MODIFIED | +Step 8: combinatorial assessment wiring |
| sl_agent/multimodal/literature_receipts.py | PHASE 1 ONLY | No Phase 2 changes |
| sl_agent/multimodal/pharmacologic_analyzer.py | PHASE 1 ONLY | No Phase 2 changes |
| sl_agent/data/depmap_loader.py | MODIFIED | +StratifiedResult, +get_checkpoint_dependency_in_context() |
| sl_agent/multimodal/tests/test_rs_benchmark.py | MODIFIED | +TestCombinatorialVulnerability (4 tests) |

## Installation
1. Overwrite the 8 files listed above (keep all other files unchanged)
2. Run tests: `python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ -v`
3. All tests must pass (existing 96 + new combinatorial tests)

## RS Score Weights (Phase 1 — unchanged)
| Feature | Weight |
|---------|--------|
| CCNE1 amplified | 0.30 |
| MYC amplified | 0.20 |
| ARID1A LOF | 0.20 |
| MSI-H | 0.15 |
| TP53 LOF | 0.10 |
| High ploidy | 0.10 |
| MBD4 LOF | 0.05 |

## Combinatorial Tier Logic
- `assess_combinatorial_vulnerability(gene, rs_score, matrix)` → Optional[CombinationVulnerability]
- Returns None when: rs_score == 0, no ATR/WEE1 axis in matrix
- combined_tier = None → data gaps exist, tier not yet realized
- combined_tier = tier string → combinatorial tier realized (CRISPR/pharma evidence present)
- Data gaps list specifies exactly what DepMap/isogenic experiments would unlock promotion

## DepMap Stratification (NEW)
```python
from sl_agent.data.depmap_loader import get_checkpoint_dependency_in_context

results = get_checkpoint_dependency_in_context(
    context_gene="MBD4",
    checkpoint_genes=["ATR", "CHEK1", "WEE1", "PKMYT1"],
    mutation_type="loss_of_function",
)
# results["ATR"].delta_dep < -0.15 → combinatorial tier unlockable
```

## The Paper Sentence This Enables
"We define MBD4 deficiency not as a single-axis targetable sensitivity, but as a systemic
replication-stress amplifier. Our multi-modal engine identifies an emergent combinatorial
vulnerability: BER-pathway collapse (MBD4-LOF) necessitates checkpoint-mediated replication
fork stabilization (ATR/WEE1). The engine explicitly models this convergence through replication
stress scoring and requires independent experimental evidence before promoting the combinatorial
axis, maintaining the same multi-source triangulation discipline applied to all benchmark genotypes."
```

Then:
```bash
cd /home/user/workspace
zip -r sl_agent_rs_patch.zip sl_patch_rs/
rm -rf sl_patch_rs/
ls -lh sl_agent_rs_patch.zip
```

---

## Final Report Required
1. Each file modified and key changes
2. Full pytest output (last 30 lines minimum)
3. Total test count (96 existing + new combinatorial tests)
4. Zip file size and confirmation it was created
5. Any import issues encountered and how they were resolved
