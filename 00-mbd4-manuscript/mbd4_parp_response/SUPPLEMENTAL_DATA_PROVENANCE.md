## MBD4 → PARP inhibitor response — Data Provenance (receipt-only)

### Cross-check rule (MANDATORY before “Verified”)

- **Verified**: Any live or extracted receipt MUST include:
  - a raw minimized excerpt copied from the underlying record(s), and
  - a path→value table derived from that same excerpt.

If the path→value table contradicts the raw excerpt for the same request identity, the section is **FAIL** and must be regenerated.

---

### 1) Axis A receipts (MBD4 cytidine-analog synthetic lethality)

**Verified** artifacts:
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/axis_a/axis_a_ic50_receipts_frozen.csv`
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/axis_a/axis_a_ic50_fold_change_summary.csv`
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/npj_2022_mbd4_cytidine_sl_receipts.json`
- Decision rule: `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/axis_a/AXIS_A_DECISION_RULE.md`

---

### 2) Axis B receipts (MBD4 hypermutator / ICI outcomes)

**Verified** artifacts (hypermutator signature):
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/hypermutator/hypermutator_signature_moesm4_clean.csv`
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/hypermutator/hypermutator_signature_summary.json`

**Partial** (case-level ICI receipt only; not a cohort table yet):
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/axis_b_ici/axis_b_ici_patient_level_receipts.csv`
- `oncology-coPilot/oncology-backend-minimal/data/literature_mbd4/out/axis_b_ici/axis_b_ici_extraction_summary.json`

---

### 3) Axis C receipts (MBD4 / RNF144A / PARP1 → PARPi response)

**Current status**: **Hardened Submission Metadata** (GDSC2/DepMap 24Q2)

**1. PARP1-PARPi Pan-Cancer Correlation (Axis 1 Context)**
- **Sample Size**: n=481 cell lines
- **Spearman Correlation**: ρ=−0.416
- **Significance**: p=1.36×10⁻²¹
- **Source**: DepMap 24Q2 Logp1TPM vs GDSC2 IC50 Z-Scores

**2. ATR Inhibition (Axis 2 ATR Checkpoint Dependency)**
- **Baseline Cohort**: n=14 LOF vs 914 WT (Methods WT rule: no somatic MBD4 in WT pool)
- **Statistics**: p=0.021, Δ=−0.74, Cohen's d=−0.50
- **MSI-H Purge (Confound Test 1)**: n=10 MSS/LOF vs 906 WT; p=0.015, Δ=−0.91, Cohen's d=−0.62
- **TP53-Stratified (Confound Test 2)**: n=11 LOF/TP53 vs 619 WT/TP53; p=0.003, Δ=−1.07, Cohen's d=−0.74

**3. MBD4-PARP1/RNF144A Mechanistic Axis (Nullified/Exploratory)**
- **PARP1 Expression**: p=0.605 (n=19 LOF vs 1,498 non-LOF expression pool; two-sided MWU; DepMap expression drop)
- **RNF144A Expression**: p≈0.48 (same cohort definition; two-sided MWU)
- **Interpretation**: MBD4-LOF does not selectively drive PARP1 transcriptional upregulation; PARPi sensitivity in MBD4-LOF must be explained by non-transcriptional trapping substrate availability or alternative mechanisms.

**Artifacts & Receipts** (paths relative to this manuscript bundle root `publications/00-mbd4-manuscript/mbd4_parp_response/`):

**Authoritative (manuscript numbers):**
- PARP1/RNF144A expression (two-sided MWU, n=19 vs 1,498 comparator, medians and p as in text): `artifacts/axis_c_preclinical/parp_axis_expression_MANUSCRIPT_RECEIPT.json`
- PARP1 vs PARPi Spearman (n=481, ρ=−0.416, p=1.36×10⁻²¹): `artifacts/axis_c_preclinical/parp1_parpi_spearman_MANUSCRIPT_RECEIPT.json`
- Ceralasertib LN_IC50 leave-one-out (14 iterations, one-sided MWU; max p matches canonical rerun): `artifacts/canonical_atr_wee1_rerun_20260405/ceralasertib_ln_ic50_leave_one_out_MANUSCRIPT_RECEIPT.json`

**Deprecated / exploratory (do not use for manuscript citation):**
- `artifacts/axis_c_preclinical/DEPRECATED_mbd4_rnf144a_expression_n8_vs_1665_exploratory.json` — alternate cohort (8 vs 1665), not the n=19 / n=1,498 expression block
- `artifacts/axis_c_preclinical/DEPRECATED_parp1_parpi_spearman_n488_exploratory_join.json` — alternate join (n=488), not the n=481 Spearman block

**Other axis C artifacts:**
- IC50/AUC summary: `artifacts/axis_c_preclinical/gdsc2_parpi_mbd4_summary.json`
- Raw excerpts: `artifacts/axis_c_preclinical/raw_excerpt_gdsc2_parpi_rows_mbd4_lof.csv`
- Canonical ATR/WEE1 rerun (frozen): `artifacts/canonical_atr_wee1_rerun_20260405/canonical_atr_wee1_rerun.json` (and `.csv`)
- **Not bundled in this package:** exploratory clinical cohort JSON (`axis_c_clinical/parpi_cohort_analysis.json`) was referenced in an earlier draft; no `artifacts/axis_c_clinical/` directory is shipped here. Prospective cohort validation remains out of scope for the frozen artifact set.
