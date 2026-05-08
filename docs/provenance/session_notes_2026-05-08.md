# BrM SL Discovery Engine — DepMap 24Q4 Live Run

## Objective
Wire real DepMap 24Q4 data into the BrM targetability matrix builder (fjkiani/Synthetic-Lethality) to produce genuine NOVEL rows. Freeze as brm_targetability_matrix_v3.json. Close the provenance chain back to crispro-backend-v3 SL handoff metadata.

## Scientific Context & Assumptions
- Disease: brain_metastasis (BrM), NSCLC context
- DepMap release: 24Q4 (figshare article 27993248)
- NSCLC filter: OncotreePrimaryDisease == "Non-Small Cell Lung Cancer" (164 model lines; 95 with CRISPR+expression overlap)
- Stratification: Q75/Q25 expression split, Mann-Whitney U, BH correction (FDR_STRONG=0.10, FDR_MODERATE=0.25, DELTA_MIN=0.10)
- BrM universe: 46 genes (brm_universe_v1.json)
- Locked SL signals: ZEB1→ITGAV (delta=−0.72, FDR=0.013), SPP1→NFE2L2 (delta=−0.73, FDR=0.0001), VIM→FERMT2 (delta=−0.49, FDR=0.004), POSTN→BACE1 (r=−0.29, p=0.004 — co-dependency only, not SL)

## Progress
**COMPLETE — all 5 steps done.**

- Step 1: DepMap 24Q4 downloaded (Expression 506.6 MB, CRISPR 428.7 MB, Model 0.6 MB)
- Step 2: Vectorized build completed in 333s — all 46 genes NOVEL
- Step 3: ZEB1/SPP1/VIM confirmed; POSTN→BACE1 correctly absent; 117/117 BrM tests passing
- Step 4: brm_targetability_matrix_v3.json frozen, committed as 64258b6 to fjkiani/Synthetic-Lethality
- Step 5: _SL_HANDOFF_MAP provenance updated, committed as 7d65571 to fjkiani/crispro-backend-v3

## Data, Methods & Parameters
- Expression: OmicsExpressionProteinCodingGenesTPMLogp1.csv (1673 × 19193, 506.6 MB)
- CRISPR: CRISPRGeneEffect.csv (1178 × 17916, 428.7 MB)
- Model: Model.csv (2105 × 46, 0.6 MB)
- Column format: "GENE (ENTREZ_ID)" → stripped to plain symbol before use
- OncotreeLineage patch: 24Q4 uses "Lung" in OncotreeLineage; engine expects "Non-Small Cell Lung Cancer". Patched via OncotreePrimaryDisease column (adapter-layer fix, no engine code change).
- Vectorized stratifier: numpy pre-filter (delta < 0, |delta| >= 0.10) before Mann-Whitney loop → ~10x speedup vs original Python loop. Identical results.
- novelty_score formula: (n_fdr10 × 2 + n_fdr10_25 × 1) / 50 (POL-002, locked v2)

## Key Findings
- All 46 BrM universe genes have DepMap SL signal → all NOVEL (0 CALIBRATION, 0 NEGATIVE)
- Receipt-backed genes (BACE1, MMP9, CCL2, ICAM1, MMP2, CLDN5, TWIST1) → score_basis=HYBRID
- BACE1 is rows[0] with highest confidence (0.887) due to HYBRID score (receipt + DepMap)
- novelty_score (depmap_novelty_index) range: 2.46–16.06 (all >1.0 due to denominator-50 calibration vs observed 120–500 partners)
- POSTN→BACE1: delta=−0.012 in POSTN-stratified NSCLC lines (below DELTA_MIN=0.10) → correctly absent from SL partners. CALIBRATION_ONLY label in crispro-backend-v3 remains valid.
- ZEB1→ITGAV: PASS | SPP1→NFE2L2: PASS | VIM→FERMT2: PASS

## Decisions & Rationale
**A+A decision (approved by Fahad Kiani):**
- Issue 1 (0 CALIBRATION): Accept. score_basis=HYBRID distinguishes receipt-backed genes. No exploit_router patch.
- Issue 2 (novelty_score >1.0): Accept. Formula locked at POL-002 (v2). No clipping, no renormalization.
  - display_label: "depmap_novelty_index"
  - scale: "unbounded positive real; not capped at 1.0"
  - interpretation (approved wording): "Unbounded weighted partner-count index derived from the locked v2 DepMap SL formula. Values greater than 1.0 are valid outputs of the formula and indicate more threshold-passing partners relative to the v2 denominator baseline. Do not interpret as a probability, percentile, or universally comparable cross-context score."
  - observed_range_v3: "2.46 – 16.06"
  - rank_order: "preserved exactly as computed; no clipping or renormalization applied"

**Wording tweak (approved by Fahad Kiani):** Added "percentile, or universally comparable cross-context score" to the interpretation string to future-proof semantics against cross-context misuse.

## Commits
- fjkiani/Synthetic-Lethality: 64258b6 — brm_targetability_matrix_v3.json + run_brm_matrix_build_v3.py
- fjkiani/crispro-backend-v3: 7d65571 — _SL_HANDOFF_MAP provenance update (sl_rationale for ITGAV/FERMT2/NFE2L2)

## Manuscript Scaffold (2026-05-07)
- File: `/mnt/results/manuscript_resurrection_scaffold_v1.md`
- Type: Perspective / Commentary, ~5,400 words
- Central claim: SL dependencies cannot be immunoedited away (ACT IV as proof)
- Six GBM trial classes decoded: cilengitide, bevacizumab, EGFR TKIs, rindopepimut, checkpoint inhibitors, mTOR/PI3K inhibitors
- Key numbers used: ZEB1→ITGAV delta=−0.7184, FDR=0.001203, n=24/24; SPP1→NFE2L2 delta=−0.7326, FDR=8e-6; VIM→FERMT2 delta=−0.4863, FDR=0.000345
- Proposed journals: Nature Medicine (Perspective), Cancer Cell (Commentary), Clinical Cancer Research (Viewpoint)
- Status: Scaffold complete — needs author list, figure generation, reference formatting, journal-specific word count trim

## Open Questions / Validation Checks
- Future: novelty_score denominator recalibration (POL-002 v3) — not in this commit; requires new policy decision
- Future: NEGATIVE rows — currently 0 because all 46 genes pass pre-filter. Consider adding genes that fail DELTA_MIN as NEGATIVE in a future universe expansion.
- Future: UI/output layer should enforce that mixed score_basis tables are not sorted by raw novelty_score alone (metadata note points in this direction but consuming code must respect it)
- crispro-backend-v3 pre-existing failures (5 tests in holistic_score + io_engine) — unrelated to our changes, pre-date this session
