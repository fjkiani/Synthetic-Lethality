# ATR/WEE1 manuscript ↔ canonical rerun reconciliation

**Canonical artifacts (frozen):**

- `canonical_atr_wee1_rerun.csv`
- `canonical_atr_wee1_rerun.json`

**Rerun script:** `rxiv/scripts/canonical_atr_wee1_rerun.py`

**WT rule (Methods):** Wild-type = no somatic MBD4 in `OmicsSomaticMutations` (any MBD4 call excluded from WT). LOF = `LikelyLoF=True` for MBD4 only.

**Alternate WT row** in CSV: WT = all models not True-LOF (MBD4 missense/passenger may sit in WT pool) — labeled in `cohort_definition`.

Rounding: “match” for p/d/Δ when manuscript rounds to the printed precision; **n** must match exactly for a strict match.

| manuscript value (location) | canonical rerun value | match? | needs edit? |
|-----------------------------|------------------------|--------|-------------|
| Ceralasertib LN_IC50 n_WT=914 (Table, Results) | n_WT=914 (Methods WT) | **Yes** | No |
| Ceralasertib LN_IC50 n_LOF=14 | n_LOF=14 | Yes | No |
| Ceralasertib LN_IC50 Δ=−0.74 | Δ=−0.732452988 | Borderline | Yes if strict −0.74 text |
| Ceralasertib LN_IC50 p=0.021 | p=0.021484496737088882 | Yes (0.021 rounded) | No |
| Ceralasertib LN_IC50 d=−0.50 | d=−0.5032867186922607 | Yes (2 dp) | No |
| Ceralasertib AUC row (Table) | AUC row in CSV row 3 | Yes (Δ/p/d within print precision) | No |
| Ceralasertib Z_SCORE row (Table) | Z_SCORE row in CSV row 4 | Yes (Δ/p/d within print precision) | No |
| Fig2 caption: “−0.74 … d=−0.50” | Same primary row as table (LN_IC50 d=−0.503) | Yes for d; Δ wording | Yes if caption insists −0.74 exact |
| Adavosertib n_LOF=15, n_WT=920 | n_LOF=15, n_WT=920 | **Yes** | No |
| Adavosertib Δ=−0.51, p=0.074, d=−0.36 | −0.508, 0.074457, −0.359 | Yes (print) | No |
| MSI purge n=10 vs n=906 | n=10 vs n=906 | **Yes** | No |
| MSI purge Δ=−0.91, p=0.015, d=−0.62 | −0.910, 0.015329, −0.623 | Yes (print) | No |
| TP53 LN_IC50 n=11 vs n=619 | n=11 vs n=619 | **Yes** | No |
| TP53 LN_IC50 Δ=−1.07, p=0.003, d=−0.74 | −1.069, 0.003003, −0.740 | Yes (print) | No |
| TP53 AUC p=0.001, d=−0.89 | p=0.000873, d=−0.889 | Yes (0.001 / −0.89 rounded) | No |
| Bowel n=5 vs n=41, p=0.13 | n=5 vs n=41, p≈0.126 | **Yes** (manuscript rounds p to 0.13) | No |
| Non-Bowel n=9 vs n=873 | n=9 vs n=873 | **Yes** | No |
| Non-Bowel Δ=−0.87, p=0.025, d=−0.60 | −0.871, 0.025330, −0.599 | Yes (print) | No |
| LOO weakest p=0.045 | LOO max p=0.045165724128583974 | Yes (0.045) | No |
| Discussion: primary ceralasertib (p=0.021, d=−0.50); TP53-stratified (p=0.003, d=−0.74) | Primary LN_IC50 d≈−0.503; TP53 LN_IC50 d≈−0.740 | **Yes** | No |

**Summary:** Under the Methods WT rule (“no somatic MBD4” in WT), manuscript denominators **914 / 920 / 906 / 619 / 41 / 873** match the canonical frozen rerun. Alternate WT rows in `canonical_atr_wee1_rerun.csv` (e.g. ceralasertib n_WT=922, adavosertib n_WT=929 under `Alternate_non_LOF_only_WT_includes_MBD4_nonLoF_mutants`) are labeled separately and are **not** the manuscript primary cohort.

**Note:** CSV row `Alternate_non_LOF_only_WT_includes_MBD4_nonLoF_mutants` retains different n_WT for sensitivity analysis; the manuscript reports Methods-primary numbers only.
