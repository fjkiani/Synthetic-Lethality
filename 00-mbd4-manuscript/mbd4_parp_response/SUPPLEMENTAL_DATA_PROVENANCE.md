## MBD4 → PARP inhibitor response — Data Provenance (receipt-only)

### Cross-check rule (MANDATORY before "Verified")

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
**Denominator version**: v3 (reconciled 2026-04-27 against `canonical_atr_wee1_rerun.json`)

#### Cohort Accounting Table (single source of truth for all manuscript WT counts)

| Stage | n_LOF | n_WT | Definition | Artifact cohort_definition |
|---|---|---|---|---|
| Primary (ceralasertib) | 14 | **914** | True-LOF vs no-somatic-MBD4 WT, post-GDSC2 ceralasertib intersection | `WT_rule=Methods_no_somatic_MBD4` |
| MSI purge (Stress Test 1) | 10 | **906** | 914 WT minus ModelSubtypeFeatures MSI-H lines | `Stress MSI purge` |
| TP53 stratification (Stress Test 2) | 11 | **619** | 914 WT intersected with TP53-mut | `Stress TP53` |
| Lineage — Bowel (Stress Test 4) | 5 | **41** | Bowel/Colorectal subset of 914 WT | `Stress Lineage Bowel` |
| Lineage — non-Bowel (Stress Test 4) | 9 | **873** | Complement of Bowel within 914 WT | `Stress Lineage non-Bowel` |
| Alternate WT (sensitivity only) | 14 | **922** | Includes MBD4 non-LOF mutants; NOT primary comparator | `WT_rule=Alternate_non_LOF_only` |
| Adavosertib (primary) | 15 | **920** | Same mutation rule, different drug intersection | adavosertib primary entry |

**Retired numbers (do not use):**
- **942** = pre-GDSC2-intersection mutation-filtered pool. Not a pharmacological comparator. Removed from manuscript_v3.md.
- **934** = stale, no artifact support. Replaced with 906.
- **625** = stale, no artifact support. Replaced with 619.
- **900** = stale, no artifact support. Replaced with 873.
- **929** = alternate WT pool for adavosertib. Replaced with 920 (primary).
- **42** (Bowel WT) = stale, no artifact support. Replaced with 41.

#### Key statistics

**1. PARP1-PARPi Pan-Cancer Correlation (Axis 1 Context)**
- **Sample Size**: n=481 cell lines
- **Spearman Correlation**: ρ=−0.416
- **Significance**: p=1.36×10⁻²¹
- **Source**: DepMap 24Q2 Logp1TPM vs GDSC2 IC50 Z-Scores

**2. ATR Inhibition (Axis 2 ATR Checkpoint Dependency)**
- **Primary cohort**: n=14 LOF vs **914** WT (post-GDSC2 ceralasertib intersection; Methods WT rule)
- **Statistics**: p=0.021, Δ=−0.73 (full precision: −0.7325), Cohen's d=−0.50
  - *Note: manuscript v1 reported Δ=−0.74 (incorrect rounding); corrected to −0.73 in manuscript_v2.md and retained in v3.*
- **MSI-H Purge (Stress Test 1)**: n=10 MSS/LOF vs **906** WT; p=0.015, Δ=−0.910, Cohen's d=−0.623
- **TP53-Stratified (Stress Test 2)**: n=11 LOF/TP53 vs **619** WT/TP53; p=0.003, Δ=−1.07, Cohen's d=−0.739
- **Lineage Bowel (Stress Test 4)**: n=5 LOF vs **41** WT; Δ=−0.72, p=0.114 (underpowered)
- **Lineage non-Bowel (Stress Test 4)**: n=9 LOF vs **873** WT; Δ=−0.871, p=0.025, Cohen's d=−0.599
- **Adavosertib (primary)**: n=15 LOF vs **920** WT; Δ=−0.508, p=0.074, Cohen's d=−0.359
- **WT-Pool Sensitivity**: n_WT=**922** (includes MBD4 non-LOF mutants); p=0.022, d=−0.501 — frozen in `canonical_atr_wee1_rerun.json` (cohort_definition: "WT_rule=Alternate_non_LOF_only_WT_includes_MBD4_nonLoF_mutants")

**3. MBD4-PARP1/RNF144A Mechanistic Axis (Nullified/Exploratory)**
- **PARP1 Expression**: p=0.605 (n=19 LOF vs 1,498 non-LOF expression pool; two-sided MWU; DepMap expression drop)
- **RNF144A Expression**: p≈0.48 (same cohort definition; two-sided MWU)
- **Interpretation**: MBD4-LOF does not selectively drive PARP1 transcriptional upregulation; PARPi sensitivity in MBD4-LOF must be explained by non-transcriptional trapping substrate availability or alternative mechanisms.

#### Artifacts & Receipts

(paths relative to manuscript bundle root `publications/00-mbd4-manuscript/mbd4_parp_response/`)

**Authoritative (manuscript numbers):**
- PARP1/RNF144A expression (two-sided MWU, n=19 vs 1,498 comparator): `artifacts/axis_c_preclinical/parp_axis_expression_MANUSCRIPT_RECEIPT.json`
- PARP1 vs PARPi Spearman (n=481, ρ=−0.416, p=1.36×10⁻²¹): `artifacts/axis_c_preclinical/parp1_parpi_spearman_MANUSCRIPT_RECEIPT.json`
- Ceralasertib LN_IC50 leave-one-out (14 iterations, max p=0.045): `artifacts/canonical_atr_wee1_rerun_20260405/ceralasertib_ln_ic50_leave_one_out_MANUSCRIPT_RECEIPT.json`
- All stress tests + WT-pool sensitivity + adavosertib (full precision, all n_wt values): `artifacts/canonical_atr_wee1_rerun_20260405/canonical_atr_wee1_rerun.json`

**Deprecated / exploratory (do not use for manuscript citation):**
- `artifacts/axis_c_preclinical/DEPRECATED_mbd4_rnf144a_expression_n8_vs_1665_exploratory.json`
- `artifacts/axis_c_preclinical/DEPRECATED_parp1_parpi_spearman_n488_exploratory_join.json`

**Other axis C artifacts:**
- IC50/AUC summary: `artifacts/axis_c_preclinical/gdsc2_parpi_mbd4_summary.json`
- Raw excerpts: `artifacts/axis_c_preclinical/raw_excerpt_gdsc2_parpi_rows_mbd4_lof.csv`
- **Not bundled:** exploratory clinical cohort JSON (`axis_c_clinical/parpi_cohort_analysis.json`) — out of scope for frozen artifact set.

---

### 4) A9 GDSC2 Pan-Cancer Benchmark (ATM-LOF vs MBD4-LOF comparator)

**Current status**: **Verified — Discussion context only** (not a primary stress test; reported in Supplementary Table S-ATM)

**Purpose**: Establishes that the published ATM-LOF + ATRi synthetic lethality is non-resolving in the same GDSC2 pan-cancer dataset, contextualizing why pan-cancer observational comparisons may fail to detect context-specific synthetic lethalities.

**Data sources**:
- GDSC2 release 8.5 (Oct 2023), DRUG_ID=1917 (ceralasertib/AZD6738), DRUG_ID=1036 (gemcitabine)
- DepMap 23Q4 OmicsSomaticMutations (ATM LOF classification: LikelyLoF=True)
- DepMap 24Q2 ModelSubtypeFeatures (MSI-H annotation for MSS-only stratum)
- MSI enrichment: ATM-LOF 13/31 MSI-H (42%); MBD4-LOF 6/14 MSI-H (43%) by ModelSubtypeFeatures

**Statistical approach**:
- One-sided Mann-Whitney U test (LOF < WT)
- BH correction applied within each gene×stratum group (6 comparisons per gene)
- WT pool: no somatic mutation in the respective gene (Methods primary definition)

**Key results**:

| Gene | Drug | Stratum | n_LOF | n_WT | Δ_LN_IC50 | p (unadj) | padj (BH) | Cohen's d | Direction |
|---|---|---|---|---|---|---|---|---|---|
| ATM | Ceralasertib | All lines | 31 | 905 | +0.121 | 0.773 | 0.773 | +0.121 | Wrong |
| ATM | Ceralasertib | MSS only | 18 | 792 | +0.348 | 0.225 | 0.450 | +0.348 | Wrong |
| ATM | Gemcitabine | All lines | 31 | 905 | +0.412 | 0.033 | 0.099 | +0.289 | Wrong |
| ATM | Gemcitabine | MSS only | 18 | 792 | +0.287 | 0.118 | 0.354 | +0.201 | Wrong |
| MBD4 | Ceralasertib | All lines | 14 | 914 | −0.732 | 0.021 | 0.063 | −0.503 | Correct |
| MBD4 | Ceralasertib | MSS only | 10 | 906 | −0.910 | 0.015 | 0.060 | −0.623 | Correct |
| MBD4 | Gemcitabine | All lines | 14 | 914 | −0.618 | 0.038 | 0.076 | −0.441 | Correct |
| MBD4 | Gemcitabine | MSS only | 10 | 906 | −0.724 | 0.029 | 0.072 | −0.512 | Correct |

**MSI definition note**:
- Primary MSI annotation (manuscript): DepMap 24Q2 `ModelSubtypeFeatures` → 6/21 MBD4-LOF MSI-H
- Sensitivity MSI annotation (A9 run): DepMap 23Q4 `MSIScore > 10` → 20/21 MBD4-LOF MSI-H

**Local artifact**: `A9_gdsc2_pharmacology_benchmark.csv` (saved to `/mnt/results/mbd4_atr_evidence/`)

**Manuscript placement**: Discussion paragraph 4; Supplementary Table S-ATM.

**Version note**: A9 analysis used DepMap 23Q4 for ATM LOF classification; manuscript primary analyses use DepMap 24Q2. GDSC2 release 8.5 is the same underlying dataset in both cases.

---

### Denominator reconciliation history

| Version | Date | Change |
|---|---|---|
| v1 (original) | pre-2026-04-27 | Mixed denominators: 942 (pre-intersection), 934 (stale), 625 (stale), 900 (stale), 929 (alternate pool), 42 (stale Bowel WT) |
| v2 (strengthening pass) | 2026-04-27 | Rounding fix (−0.74→−0.73), MSI dual-definition note, A9 Discussion paragraph, Limitations paragraph added. Denominators not yet reconciled. |
| v3 (denominator reconciliation) | 2026-04-27 | All WT counts corrected to post-intersection artifact values. Stale numbers retired. Methods cohort-accounting footnote added. Freeze candidate. |
