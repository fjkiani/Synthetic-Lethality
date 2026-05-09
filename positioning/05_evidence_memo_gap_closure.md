# Evidence Memo: Validation Gap Closure
**Version:** 2.0 — 2026-05-09
**Prepared by:** CrisPRO Research Engine
**Purpose:** Structured evidence summary for two validation gaps. v1.0 informed briefs v3.0. v2.0 adds in silico upgrade results: DepMap mutation/CN stratification, cBioPortal BrM frequency, TCGA LUAD co-expression, GSE271259 BrM co-expression. Informs v4.0 brief updates.
**Status:** Final — ready for counsel review

---

## Executive Summary

| Gap | Hypothesis | DepMap Result (v1.0) | In Silico Upgrade (v2.0) | Verdict |
|---|---|---|---|---|
| Gap 1 | PTEN-null NSCLC → ITGAV dependency | **NULL** (mRNA proxy: delta=−0.040, p=0.73, n=95) | **NULL confirmed** (mutation+CN: delta=+0.023, p=1.000, n=4 PTEN-loss); PTEN-loss freq in 322 NSCLC BrM: 7.5% overall, LUSC 20.7% (cBioPortal) | **Brief #2 remains COMPUTATIONAL HYPOTHESIS. All DepMap proxies null. Target population defined (7.5% BrM). Isogenic models required.** |
| Gap 3 | ZEB1-high NSCLC BrM → αV integrin + PD-1 combination | N/A (combination not testable in DepMap) | TCGA LUAD ZEB1/CD274: r=0.253, p=6.97×10⁻⁹, n=510 (POSITIVE, primary); GSE271259 BrM: r=0.050, p=0.780, n=34 (NULL, BrM-specific) | **Brief #3 upgraded for primary NSCLC (TCGA quantitative). BrM-specific co-expression not detected in available bulk RNA-seq. Gap partially closed.** |

---

## GAP 1: PTEN-Null NSCLC → ITGAV Dependency

### Hypothesis
PTEN loss in NSCLC BrM creates convergent vulnerability to αV integrin-mediated survival signaling, making PTEN-null tumors differentially sensitive to αV integrin inhibition (cilengitide) compared to PTEN-intact tumors.

### Tier A: DepMap 24Q4 Internal Analysis (COMPLETED — NULL RESULT)

**Dataset:** DepMap 24Q4, figshare article 27993248
- Expression: OmicsExpressionProteinCodingGenesTPMLogp1.csv (figshare ID 51065489)
- CRISPR: CRISPRGeneEffect.csv (figshare ID 51064667)
- Model: Model.csv (figshare ID 51065297)
- NSCLC filter: OncotreePrimaryDisease == "Non-Small Cell Lung Cancer"
- Overlap (expression + CRISPR): n = 95 NSCLC lines

**Method:** Q25/Q75 PTEN mRNA expression split → Mann-Whitney U test (two-sided) on ITGAV CRISPR gene effect scores

**Results:**

| Stratum | n | ITGAV CRISPR mean | ITGAV CRISPR SD |
|---|---|---|---|
| PTEN-low (≤Q25) | 24 | −0.7406 | 0.4821 |
| PTEN-high (≥Q75) | 24 | −0.7008 | 0.4612 |
| **Delta (low − high)** | — | **−0.0398** | — |

- Mann-Whitney U p-value (two-sided): **0.7337** (NOT SIGNIFICANT)
- Pearson r (PTEN expr vs ITGAV CRISPR): −0.071 (p = 0.494)
- Spearman r (PTEN expr vs ITGAV CRISPR): +0.050 (p = 0.629)
- Delta threshold (DELTA_MIN): 0.10 — **NOT MET** (observed delta = 0.040)

**Positive control (same pipeline, same cohort):**
- ZEB1→ITGAV: delta = −0.7184, p = 0.000002 (one-sided) — **REPLICATES EXACTLY**
- Confirms the pipeline is functioning correctly; the null result for PTEN is real, not a technical failure.

**Interpretation:**
PTEN mRNA expression level does NOT stratify ITGAV CRISPR dependency in DepMap 24Q4 NSCLC lines (n=95). Both PTEN-low and PTEN-high strata show substantial ITGAV dependency (means −0.74 and −0.70, respectively), suggesting ITGAV dependency in NSCLC is driven by ZEB1/EMT status rather than PTEN status. PTEN loss does not add incremental ITGAV dependency on top of the baseline.

**Important caveat:** This analysis uses PTEN mRNA expression as a proxy for PTEN functional status. PTEN is frequently lost by protein-level mechanisms (post-translational degradation, nuclear exclusion, promoter methylation) without corresponding mRNA reduction. A Q25 mRNA split may not capture the truly PTEN-null lines. A more rigorous test would require:
1. PTEN protein IHC data or somatic mutation/deletion status (OmicsSomaticMutationsMatrixDamaging.csv)
2. Isogenic PTEN-KO NSCLC cell lines (e.g., PTEN KO in H1299, A549) with ITGAV CRISPR KO or cilengitide treatment

**Verdict:** The DepMap mRNA-based stratification does NOT support the PTEN-null→ITGAV hypothesis. This is the strongest available computational evidence on this question, and it is null. The hypothesis remains unvalidated.

---

### Tier B: Published Cell Line Data — PTEN Loss + αV Integrin Sensitivity

**Search conducted:** 2026-05-09. Queries: "PTEN loss ITGAV integrin sensitivity NSCLC", "PTEN null cilengitide sensitivity cancer cell lines", "PTEN loss integrin αV dependency lung cancer"

| Paper | Finding | Relevance | Limitation |
|---|---|---|---|
| Girnius N et al. (2024). *Breast Cancer Res* 26:142. DOI: 10.1186/s13058-024-01942-2 | DepMap-identified ITGAV vulnerability + cilengitide sensitivity in cancer cell lines; low pan-integrin expression as sensitivity biomarker | Cilengitide sensitivity biomarker concept | **Breast cancer (TNBC), not NSCLC. Biomarker is pan-integrin expression level, not PTEN status.** |
| Expósito F et al. (2023). *Cancer Res* 83(14):2355–2370. DOI: 10.1158/0008-5472.can-22-3023 | PTEN loss in >25% NSCLC; higher PD-L1 and TGFβ; anti-PD-1 resistance | PTEN loss prevalence in NSCLC confirmed | Immunotherapy resistance endpoint, not integrin dependency |

**Verdict:** No published paper directly demonstrates that PTEN-null NSCLC cells are specifically sensitive to ITGAV inhibition or cilengitide. The Girnius 2024 paper is the closest analog but is in breast cancer and uses a different biomarker (pan-integrin expression, not PTEN status).

---

### Tier C: Mechanistic Literature — PTEN Loss → FAK/Integrin Crosstalk

| Paper | Finding | Relevance | Limitation |
|---|---|---|---|
| Cavazzoni A et al. (2017). *Oncotarget* 8(43):74107–74123. DOI: 10.18632/oncotarget.18087 | PTEN loss correlates with FAK phosphorylation in NSCLC SCC patient specimens; AKT+FAK combined inhibition is synergistic only in PTEN-reduced cells | **STRONGEST** — direct PTEN→FAK link in NSCLC patient tissue | FAK inhibitor endpoint, not ITGAV inhibitor. FAK is downstream of integrin signaling but is not the same target as ITGAV. |
| Peng W et al. (2016). *Cancer Cell* 29(5):636–649. DOI: 10.1016/j.ccell.2016.04.003 | PTEN loss → immune exclusion via PI3K/AKT → β-catenin → CCL4 suppression | PTEN loss → immune evasion mechanism | Immunotherapy resistance, not integrin |

**Mechanistic chain (supported):** PTEN loss → PI3K/AKT constitutive activation → FAK phosphorylation (Cavazzoni 2017) → integrin-mediated survival signaling (canonical)

**Mechanistic chain (gap):** FAK phosphorylation → ITGAV-specific dependency (NOT demonstrated in NSCLC; FAK is downstream of multiple integrins, not ITGAV-specific)

**Verdict:** The mechanistic rationale for PTEN→FAK→integrin signaling is supported in NSCLC (Cavazzoni 2017). However, the specific step from FAK activation to ITGAV-specific dependency is not demonstrated. FAK is activated by multiple integrin heterodimers; PTEN loss does not specifically upregulate ITGAV.

---

### Gap 1 Summary and Brief #2 Recommendation

**Evidence tier (updated):** COMPUTATIONAL HYPOTHESIS — DepMap NULL result, mechanistic rationale only

**What to add to Brief #2 specification:**
1. DepMap null result must be disclosed: "DepMap 24Q4 NSCLC stratification (n=95, Q25/Q75 PTEN mRNA expression split) shows no differential ITGAV CRISPR dependency (delta=−0.040, p=0.73). The PTEN-null→ITGAV dependency hypothesis is NOT supported by mRNA-based computational validation."
2. Cavazzoni 2017 should be added as mechanistic support for PTEN→FAK link in NSCLC SCC.
3. Validation path must be stated: isogenic PTEN-KO NSCLC models with ITGAV CRISPR KO or cilengitide treatment.
4. The method claim remains patentable — the null DepMap result does not invalidate the claim, but the specification must accurately represent the evidence tier.

**Filing recommendation:** File with Brief #1 as planned. The method claim is independent of validation status. The specification must not represent the PTEN-null claim as computationally validated.

---

## GAP 3: ZEB1-High NSCLC BrM → αV Integrin + PD-1 Combination

### Hypothesis
ZEB1-high NSCLC BrM tumors are simultaneously: (1) ITGAV-dependent (ZEB1→ITGAV SL, CRISPR-validated), AND (2) PD-L1-high (ZEB1→PD-L1 transcriptional activation). One biomarker — ZEB1 expression — identifies patients who should respond to both αV integrin inhibition and PD-1 checkpoint blockade. The combination should be synergistic: cilengitide disrupts integrin-mediated survival signaling while PD-1 blockade reinvigorates T cell killing.

### Tier A: ZEB1→PD-L1 in NSCLC — Evidence Map

| Paper | Finding | Evidence Type | Relevance | Limitation |
|---|---|---|---|---|
| Chen L et al. (2014). *Nat Commun* 5:5241. DOI: 10.1038/ncomms6241 | ZEB1 relieves miR-200 repression of PD-L1 in lung cancer datasets; ZEB1→PD-L1 axis drives metastasis and immune evasion | Mechanistic + lung cancer datasets | **LANDMARK** — ZEB1→PD-L1 in lung cancer | Mechanism via miR-200 de-repression, not direct promoter binding |
| Kim HR et al. (2016). *Hum Pathol* 58:45–52. DOI: 10.1016/j.humpath.2016.07.007 | EMT phenotype (ZEB1+) correlates with PD-L1 positivity in 477 lung adenocarcinoma patient specimens | NSCLC patient specimens (n=477) | **NSCLC patient data** — ZEB1/PD-L1 co-expression confirmed | ZEB1 as EMT marker, not direct transcriptional driver in this study |
| Wirbel J et al. (2025). *Cancer Immunol Immunother* 74:78. DOI: 10.1007/s00262-025-03978-5 | ZEB1 directly binds CD274 (PD-L1) promoter; drives PD-L1 transcription; ZEB1/PD-L1 co-expression in patient specimens | Direct promoter binding (ChIP-seq) | Strongest mechanistic evidence for direct ZEB1→PD-L1 transcription | **Melanoma, not NSCLC** |
| Bouillez A et al. (2017). *Oncogene* 36(32):4617–4628. DOI: 10.1038/onc.2017.47 | MUC1-C activates NF-κB→ZEB1→PD-L1 pathway in NSCLC cell lines | NSCLC cell lines | ZEB1→PD-L1 in NSCLC via NF-κB pathway | Indirect (MUC1-C upstream); not direct ZEB1 overexpression experiment |

**Synthesis:** The ZEB1→PD-L1 link in NSCLC is supported by two independent lines of evidence: (1) Chen 2014 (mechanistic, lung cancer datasets, miR-200 pathway), and (2) Kim 2016 (NSCLC patient specimens, n=477, ZEB1/PD-L1 co-expression). Wirbel 2025 provides the strongest mechanistic evidence (direct promoter binding) but is in melanoma. Bouillez 2017 provides NSCLC cell line support via NF-κB pathway. The combination of these papers constitutes multi-tier support for the ZEB1→PD-L1 axis in NSCLC.

**Evidence tier:** STRONG — multi-paper, multi-mechanism, NSCLC-specific patient data available

---

### Tier B: ZEB1 Upregulation in NSCLC BrM Patient Specimens

| Paper | Finding | Evidence Type | Relevance | Limitation |
|---|---|---|---|---|
| Lee JY et al. (2020). *Cancers* 12(12):3632. DOI: 10.3390/cancers12123632 | ZEB1 immunoreactivity: 15% in primary lung adenocarcinoma vs 55.9% in matched BrM (p=0.002, n=46 paired specimens) | Paired patient specimens (n=46) | **ZEB1 upregulation in NSCLC BrM confirmed** — strongest available BrM-specific data | Does NOT measure PD-L1 co-expression in same specimens |
| Kim S et al. (2025). *Neuro-Oncology* 27(Suppl 1):i17. DOI: 10.1093/neuonc/noaf201.0017 | ZEB1 upregulation in NSCLC BrM confirmed in independent cohort | Independent validation | Replication of Lee 2020 finding | Conference abstract; full paper pending |

**Synthesis:** ZEB1 upregulation in NSCLC BrM is confirmed in paired patient specimens (n=46, p=0.002). This is the strongest available evidence that ZEB1-high is a relevant biomarker in the NSCLC BrM context. **Critical gap:** Neither paper measures PD-L1 co-expression in the same BrM specimens. The ZEB1+PD-L1 co-expression in NSCLC BrM is inferred from (a) ZEB1 upregulation in BrM + (b) ZEB1→PD-L1 correlation in primary NSCLC — not directly measured in BrM tissue.

**Evidence tier:** MODERATE — ZEB1 upregulation in BrM confirmed; PD-L1 co-expression in BrM inferred, not directly measured

---

### Tier C: αV Integrin Inhibitor + PD-1 Blockade Combination — Preclinical Evidence

| Paper | Finding | Evidence Type | Relevance | Limitation |
|---|---|---|---|---|
| Pan Y et al. (2022). *Bioengineered* 13(2):3928–3940. DOI: 10.1080/21655979.2022.2029236 | Cilengitide + anti-PD-1 reduces tumor growth, extends survival in murine melanoma; cilengitide decreases PD-L1 expression via STAT3 pathway | Murine syngeneic model | **Cilengitide + anti-PD-1 combination preclinical proof-of-concept** | **Melanoma model, not NSCLC or BrM** |
| Ellert-Miklaszewska A et al. (2025). *J Exp Clin Cancer Res* 44:103. DOI: 10.1186/s13046-025-03393-9 | RGD integrin-blocking peptide + anti-PD-1 reduces glioma growth, restores hot TME (increased CD8+ T cell infiltration) | Murine glioma model | Integrin + checkpoint combination in CNS tumor context | **Glioma (primary brain tumor), not NSCLC BrM** |

**Synthesis:** The cilengitide + anti-PD-1 combination has preclinical proof-of-concept in melanoma (Pan 2022). The integrin + checkpoint combination has CNS-context proof-of-concept in glioma (Ellert-Miklaszewska 2025). Neither paper tests the combination in NSCLC BrM. The mechanistic rationale for the combination (cilengitide reduces PD-L1 via STAT3; integrin blockade restores T cell infiltration) is supported by these papers but requires validation in NSCLC BrM models.

**Evidence tier:** MODERATE — combination preclinical proof-of-concept in adjacent tumor types; no NSCLC BrM-specific combination data

---

### Gap 3 Summary and Brief #3 Recommendation

**Evidence tier (updated):** MECHANISTICALLY GROUNDED — multi-tier support, no NSCLC BrM-specific combination data

**What to add to Brief #3 specification:**

*ZEB1→PD-L1 axis (upgrade from v2.0):*
- Replace generic "ZEB1→PD-L1 transcriptional activation" with specific citations:
  - Chen 2014 (Nat Commun, DOI: 10.1038/ncomms6241): ZEB1→PD-L1 via miR-200 in lung cancer datasets
  - Kim 2016 (Hum Pathol, DOI: 10.1016/j.humpath.2016.07.007): ZEB1/PD-L1 co-expression in 477 NSCLC patient specimens
  - Bouillez 2017 (Oncogene, DOI: 10.1038/onc.2017.47): ZEB1→PD-L1 via NF-κB in NSCLC cell lines
  - Wirbel 2025 (Cancer Immunol Immunother, DOI: 10.1007/s00262-025-03978-5): Direct ZEB1→CD274 promoter binding (melanoma — note cancer type)

*ZEB1 upregulation in BrM (add):*
  - Lee 2020 (Cancers, DOI: 10.3390/cancers12123632): ZEB1 15% → 55.9% in matched NSCLC BrM (p=0.002, n=46 paired)

*Combination preclinical (add with honest cancer-type labels):*
  - Pan 2022 (Bioengineered, DOI: 10.1080/21655979.2022.2029236): Cilengitide + anti-PD-1 in murine melanoma (not NSCLC)
  - Ellert-Miklaszewska 2025 (J Exp Clin Cancer Res, DOI: 10.1186/s13046-025-03393-9): Integrin peptide + anti-PD-1 in murine glioma (CNS context, not NSCLC BrM)

*Honest gap statement (add):*
  - "ZEB1 and PD-L1 co-expression has not been directly measured in the same NSCLC BrM specimens. The co-expression is inferred from (a) ZEB1 upregulation in NSCLC BrM (Lee 2020, n=46) and (b) ZEB1/PD-L1 correlation in primary NSCLC (Kim 2016, n=477). Direct measurement in NSCLC BrM tissue is a validation priority."
  - "The αV integrin inhibitor + PD-1 blockade combination has not been tested in NSCLC BrM models. Validation in ZEB1-high NSCLC BrM syngeneic or PDX models is required before clinical translation."

**Filing recommendation:** File Q2 2026 as planned. The combination claim is mechanistically grounded with NSCLC-specific ZEB1→PD-L1 evidence (Chen 2014, Kim 2016). The specification must accurately represent the evidence tier and the gap in BrM-specific combination data.

---

## Validation Roadmap

### Gap 1 (PTEN→ITGAV) — Required to Upgrade Brief #2

| Priority | Experiment | Status | Expected Timeline | Outcome if Positive |
|---|---|---|---|---|
| 1 | DepMap PTEN mutation/CN stratification (OmicsSomaticMutationsMatrixDamaging + OmicsAbsoluteCNGene) | **COMPLETED (v2.0) — NULL** | — | n=4 PTEN-loss; all proxies null; pipeline confirmed |
| 1B | PTEN frequency in NSCLC BrM (cBioPortal bm_nsclc_mskcc_2023) | **COMPLETED (v2.0) — POSITIVE** | — | 7.5% overall; LUSC 20.7%; target population defined |
| 2 | Isogenic PTEN-KO NSCLC cell lines (H1299, A549) + cilengitide dose-response | Pending | 3–6 months | Upgrade to CRISPR-validated SL |
| 3 | PTEN-null NSCLC patient-derived organoids + cilengitide | Pending | 6–12 months | Clinical translation evidence |

### Gap 3 (ZEB1→PD-L1 + Combination) — Required to Upgrade Brief #3

| Priority | Experiment | Status | Expected Timeline | Outcome if Positive |
|---|---|---|---|---|
| 1 | TCGA LUAD ZEB1/CD274 co-expression (cBioPortal, n=510) | **COMPLETED (v2.0) — POSITIVE** | — | r=0.253, p=6.97×10⁻⁹; quantitative RNA-seq upgrade |
| 1B | GSE271259 ZEB1/CD274 in NSCLC BrM (n=34 MET-amp) | **COMPLETED (v2.0) — NULL** | — | r=0.050, p=0.780; MET-amp cohort confound; underpowered |
| 1C | GSE223503 NSCLC BrM Atlas (snRNA-seq) | **DROPPED — snRNA-seq, not executable** | — | Dataset is snRNA-seq; requires cell-type deconvolution |
| 2 | ZEB1/PD-L1 co-IHC in NSCLC BrM tissue microarray | Pending | 3–6 months | Direct co-expression evidence in BrM |
| 3 | ZEB1-high NSCLC BrM syngeneic model + cilengitide + anti-PD-1 | Pending | 6–12 months | Combination preclinical proof in NSCLC BrM |
| 4 | ZEB1 overexpression in NSCLC cell lines → PD-L1 upregulation (Western, flow) | Pending | 1–3 months | Direct mechanistic validation in NSCLC |

---

*Evidence Memo v2.0 — 2026-05-09 | CrisPRO Research Engine | Confidential — Attorney-Client Privileged when transmitted to counsel*
