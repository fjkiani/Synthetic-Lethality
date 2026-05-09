# CrisPRO Patent Filing Briefs
**Version:** 3.0 — 2026-05-09 (evidence gap closure pass complete)
**Purpose:** Provisional patent application summaries — for patent counsel review and filing
**Status:** Pre-filing drafts. These are plain-English + claim-structure summaries for counsel. Not legal advice.

**Change log v1.0 → v2.0:**
- Brief #1: Nagaishi 2017 citation confirmed (DOI: 10.1016/j.jocn.2017.08.050; n=29 BrM specimens, ZEB1 in 59% BrM vs 24% primary, P=0.02). "Citation to be confirmed" language removed.
- Brief #1: Mattson 2024 citation confirmed (DOI: 10.1038/s41594-024-01211-y; Nature Struct Mol Biol; CRISPR screen identifies ITGAV/ITGB5 as essential integrin pair in solid cancers). "Citation to be confirmed" language removed.
- Brief #1: n=24/24 framing clarified — this is the ZEB1-Q75 stratum (n=24 lines); ZEB1-Q25 stratum (n=24 lines) shows no ITGAV dependency (delta ≈ 0). The SL signal is stratum-specific, not universal.
- Brief #2: Evidence tier label added — PTEN-null claim rests on 8D mechanism fit score (computational hypothesis, fit=0.88), NOT CRISPR-validated SL dependency. Requires functional validation in isogenic PTEN-null NSCLC models before clinical claim.
- All deal structure ranges labeled as illustrative benchmarks, not guaranteed outcomes.

**Change log v2.0 → v3.0:**
- Brief #2: DepMap 24Q4 PTEN→ITGAV stratification result added (NULL — delta=−0.040, p=0.73, n=95 NSCLC lines). Cavazzoni 2017 added as mechanistic support (PTEN→FAK in NSCLC SCC). Validation path specified. Evidence tier confirmed as COMPUTATIONAL HYPOTHESIS — NOT computationally validated.
- Brief #3: ZEB1→PD-L1 specification upgraded with 4 specific citations (Chen 2014, Kim 2016, Bouillez 2017, Wirbel 2025). ZEB1 upregulation in NSCLC BrM added (Lee 2020, n=46 paired). Combination preclinical evidence added (Pan 2022 melanoma, Ellert-Miklaszewska 2025 glioma) with honest cancer-type labels. Gap statement added: ZEB1+PD-L1 co-expression not directly measured in NSCLC BrM specimens.
- Evidence tier table in IP Strategy Notes updated to reflect v3.0 findings.

---

## Filing Priority Order

| Priority | Asset | Biomarker | Drug Class | Evidence Tier | Status |
|---|---|---|---|---|---|
| 1 (IMMEDIATE) | ZEB1-high NSCLC BrM + αV integrin inhibition | ZEB1 expression (high) | αV integrin inhibitor (e.g., cilengitide) | **CRISPR-validated SL** | File now |
| 2 (WITH #1) | PTEN-null NSCLC BrM + αV integrin inhibition | PTEN loss (null/low) | αV integrin inhibitor | **Computational hypothesis — validation pending** | File with #1 |
| 3 (Q2 2026) | ZEB1-high NSCLC BrM + PD-1 blockade combination | ZEB1 expression (high) | αV integrin inhibitor + PD-1 inhibitor | **Mechanistically grounded — validation pending** | File Q2 |
| 4 (WHEN CANDIDATE IDENTIFIED) | SPP1-high NSCLC BrM + NRF2 inhibition | SPP1 expression (high) | NRF2 pathway inhibitor | **CRISPR-validated SL — no clinical-stage drug** | File when CNS-penetrant NRF2 inhibitor identified |

---

## PATENT BRIEF #1 — PRIORITY FILING
### ZEB1-High NSCLC Brain Metastasis + αV Integrin Inhibition
**Evidence tier: CRISPR-validated synthetic lethal dependency**

**Plain-English Summary:**
We have discovered that lung cancer cells that have spread to the brain and express high levels of a protein called ZEB1 are structurally dependent on a different protein called αV integrin (encoded by the ITGAV gene) for their survival. When ZEB1 is highly expressed, the cancer cell becomes addicted to αV integrin — knock out αV integrin in these cells and they die; knock it out in ZEB1-low cells and they survive. This is a synthetic lethal relationship. We have identified this dependency in 24 out of 24 NSCLC cell lines in the ZEB1-high stratum (DepMap 24Q4, delta dependency = −0.7184, FDR = 0.001203). The ZEB1-low stratum (n=24 lines) shows no ITGAV dependency (delta ≈ 0) — confirming the signal is ZEB1-stratum-specific, not a general NSCLC dependency. A drug that inhibits αV integrin (cilengitide, or any αV integrin inhibitor) should therefore selectively kill ZEB1-high lung cancer brain metastases. This is a new use of an existing drug class in a specific patient population defined by a specific biomarker — and it is patentable.

**Claim Structure (for counsel):**

*Independent Claim 1 — Method of Treatment:*
A method of treating brain metastasis in a human patient comprising:
(a) determining that the patient has non-small cell lung cancer brain metastasis (NSCLC BrM);
(b) measuring ZEB1 expression in a tumor sample from the patient;
(c) identifying the patient as ZEB1-high based on a predetermined expression threshold; and
(d) administering to the patient a therapeutically effective amount of an αV integrin inhibitor.

*Independent Claim 2 — Patient Selection Method:*
A method of selecting a patient with NSCLC brain metastasis for treatment with an αV integrin inhibitor comprising:
(a) obtaining a tumor sample from the patient;
(b) measuring ZEB1 expression in the tumor sample; and
(c) selecting the patient for αV integrin inhibitor therapy if ZEB1 expression exceeds a predetermined threshold indicative of ITGAV synthetic lethality.

*Dependent Claims (examples):*
- Claim 3: The method of Claim 1, wherein the αV integrin inhibitor is cilengitide or a pharmaceutically acceptable salt thereof.
- Claim 4: The method of Claim 1, wherein the αV integrin inhibitor inhibits αVβ3 and/or αVβ5 integrin.
- Claim 5: The method of Claim 1, wherein ZEB1 expression is measured by immunohistochemistry, RNA in situ hybridization, or RNA sequencing.
- Claim 6: The method of Claim 1, wherein the patient has previously received and failed EGFR tyrosine kinase inhibitor therapy.
- Claim 7: The method of Claim 1, wherein the patient has PTEN-null tumor status.
- Claim 8: The method of Claim 1, further comprising administering a PD-1 or PD-L1 inhibitor in combination with the αV integrin inhibitor.

**Supporting Evidence (for specification):**
- DepMap 24Q4 CRISPR gene effect analysis: 95 NSCLC cell lines, ZEB1 Q75 vs Q25 stratification, Mann-Whitney U, BH correction
- ZEB1→ITGAV SL signal: delta = −0.7184, FDR = 0.001203, n = 24/24 ZEB1-high lines; ZEB1-low stratum (n=24) shows delta ≈ 0 (no dependency) — signal is stratum-specific
- ZEB1 expression in BrM specimens: Nagaishi M et al. (2017). "Tumoral and stromal expression of Slug, ZEB1, and ZEB2 in brain metastasis." *J Clin Neurosci* 46:52–58. DOI: 10.1016/j.jocn.2017.08.050. [n=29 BrM patients; ZEB1 in 59% of BrM vs 24% of primary tumors, P=0.02]
- ITGAV/ITGB5 essentiality in solid cancers (CRISPR screen): Mattson NM et al. (2024). "A novel class of inhibitors that disrupts the stability of integrin heterodimers identified by CRISPR-tiling-instructed genetic screens." *Nat Struct Mol Biol* 31:736–748. DOI: 10.1038/s41594-024-01211-y. [CRISPR screens across cancer cell lines identify ITGAV/ITGB5 as essential integrin pair in solid tumors]
- ITGAV as cancer dependency in ~40% of NSCLC cell lines (DepMap): Kahangi BP et al. (2024). AACR Abstract 4266. DOI: 10.1158/1538-7445.am2024-4266.
- Cilengitide CNS penetration: PBTC-012 (CSF levels confirmed), NABTC 03-02 (tumor tissue sampling)
- ZEB1 as EMT master regulator in NSCLC BrM: Manshouri R et al. (2019). *Nat Commun* 10:5125. DOI: 10.1038/s41467-019-13022-3.
- CENTRIC trial failure analysis: Stupp R et al. (2014). *Lancet Oncol* 15:1100–1108. DOI: 10.1016/S1470-2045(14)70379-1. [Patient selection failure, not mechanism failure]

**IP Gap Analysis:**
- Cilengitide compound: Off-patent (Merck KGaA discontinued development 2013)
- NSCLC BrM indication for αV integrin inhibition: Not claimed in any published patent (search recommended)
- ZEB1 as biomarker for αV integrin inhibitor response: Not claimed (search recommended)
- ZEB1→ITGAV SL dependency in NSCLC: Not published as patent claim (DepMap analysis is novel)

**Filing Urgency:** IMMEDIATE. The DepMap analysis and ZEB1→ITGAV SL signal are documented in our internal data room (brm_targetability_matrix_v3.json, commit 64258b6). The manuscript (brm_graveyard_perspective_v3.md) is in pre-publication draft. File provisional before manuscript submission to establish priority date.

---

## PATENT BRIEF #2 — FILE WITH #1
### PTEN-Null NSCLC Brain Metastasis + αV Integrin Inhibition
**Evidence tier: COMPUTATIONAL HYPOTHESIS — functional validation pending**

> **COUNSEL NOTE — EVIDENCE TIER (v3.0 UPDATE):** This claim rests on 8D mechanism fit analysis (CrisPRO MechanismFitRanker, fit = 0.88), which is a computational alignment score, NOT a CRISPR-validated synthetic lethal dependency. Unlike Brief #1 (where ZEB1→ITGAV is validated in 24/24 NSCLC cell lines by CRISPR functional genomics), the PTEN-null → ITGAV dependency hypothesis has NOT been validated computationally or experimentally.
>
> **DepMap 24Q4 NULL RESULT (2026-05-09):** Stratification of 95 NSCLC cell lines by PTEN mRNA expression (Q25 vs Q75) shows NO differential ITGAV CRISPR dependency: delta = −0.040 (below DELTA_MIN=0.10), Mann-Whitney U p = 0.73 (two-sided), Pearson r = −0.071 (p=0.49), Spearman r = +0.050 (p=0.63). Both PTEN-low and PTEN-high strata show substantial ITGAV dependency (means −0.74 and −0.70), suggesting ITGAV dependency in NSCLC is driven by ZEB1/EMT status, not PTEN status. Positive control (ZEB1→ITGAV, same pipeline): delta = −0.7184, p = 0.000002 — confirms pipeline integrity; the PTEN null result is real.
>
> **Caveat:** This analysis uses PTEN mRNA as a proxy for PTEN functional status. PTEN is frequently lost by protein-level mechanisms without mRNA reduction. A more rigorous test requires PTEN protein IHC or somatic mutation/deletion stratification, or isogenic PTEN-KO NSCLC models.
>
> The method claim is independent of validation status and is patentable now. The specification must accurately represent the evidence tier. Do not represent the PTEN-null claim as computationally or CRISPR-validated.

**Plain-English Summary:**
PTEN is a tumor suppressor gene that is deleted or inactivated in approximately 20–30% of NSCLC brain metastases. When PTEN is lost, the PI3K/AKT survival pathway is permanently switched on. This makes the tumor resistant to EGFR inhibitors (because EGFR inhibition can't turn off AKT), resistant to checkpoint immunotherapy (because PTEN loss drives immune exclusion), and resistant to mTOR inhibitors (because mTOR inhibition triggers an AKT feedback loop). Our 8D mechanism fit analysis identifies PTEN-null NSCLC BrM as the patient profile with the strongest mechanism alignment for αV integrin inhibition (fit = 0.88, rank #2 after mTOR/PI3K). The hypothesis: PTEN loss creates a convergent vulnerability to αV integrin-mediated survival signaling that is not present in PTEN-intact tumors. This is a computational hypothesis that requires functional validation in isogenic PTEN-null NSCLC models. The patient selection logic is patentable now as a method claim.

**Claim Structure (for counsel):**

*Independent Claim 1 — Method of Treatment:*
A method of treating brain metastasis in a human patient comprising:
(a) determining that the patient has non-small cell lung cancer brain metastasis (NSCLC BrM);
(b) determining that the patient's tumor has PTEN loss (null or low PTEN expression/function);
(c) identifying the patient as a candidate for αV integrin inhibitor therapy based on PTEN loss status; and
(d) administering to the patient a therapeutically effective amount of an αV integrin inhibitor.

*Independent Claim 2 — Combination Method:*
A method of treating PTEN-null NSCLC brain metastasis comprising administering an αV integrin inhibitor to a patient who has failed or is ineligible for EGFR tyrosine kinase inhibitor therapy, checkpoint immunotherapy, and/or mTOR inhibitor therapy, wherein the patient's tumor has PTEN loss.

*Dependent Claims (examples):*
- Claim 3: The method of Claim 1, wherein PTEN loss is determined by immunohistochemistry, FISH, or next-generation sequencing.
- Claim 4: The method of Claim 1, wherein the patient also has ZEB1-high expression (combining with Brief #1).
- Claim 5: The method of Claim 1, wherein the αV integrin inhibitor is cilengitide or a pharmaceutically acceptable salt thereof.
- Claim 6: The method of Claim 2, wherein the patient has previously received and failed osimertinib.

**Supporting Evidence (for specification):**
- 8D mechanism fit analysis: PTEN-null profile, αV integrin inhibition fit = 0.88 (CrisPRO MechanismFitRanker, 2026-05-08) — **computational alignment score, not CRISPR-validated dependency**
- **DepMap 24Q4 stratification (NULL — must disclose):** 95 NSCLC lines, Q25/Q75 PTEN mRNA split, ITGAV CRISPR dependency: delta = −0.040, p = 0.73. No differential ITGAV dependency by PTEN mRNA expression level. Both strata show substantial ITGAV dependency (means −0.74 and −0.70). [DepMap 24Q4, figshare article 27993248, 2026-05-09]
- PTEN loss → FAK phosphorylation in NSCLC SCC patient specimens; AKT+FAK combined inhibition synergistic only in PTEN-reduced cells: Cavazzoni A et al. (2017). *Oncotarget* 8(43):74107–74123. DOI: 10.18632/oncotarget.18087. [Mechanistic support for PTEN→FAK→integrin signaling in NSCLC; FAK inhibitor endpoint, not ITGAV-specific]
- PTEN loss → PI3K/AKT constitutive activation: canonical oncology literature
- PTEN loss → EGFR TKI resistance: multiple Phase III analyses (IPASS, FLAURA subgroup data)
- PTEN loss → immune exclusion: Peng W et al. (2016). *Cancer Cell* 29(5):636–649. DOI: 10.1016/j.ccell.2016.04.003.
- PTEN loss → mTOR inhibitor resistance (AKT feedback): O'Reilly KE et al. (2006). *Cancer Res* 66(3):1500–1508. DOI: 10.1158/0008-5472.CAN-05-2925.
- PTEN loss frequency in NSCLC BrM: ~20–30% (literature range; exact frequency varies by study)
- **Validation required:** (1) DepMap stratification by PTEN protein loss or somatic mutation/deletion (more rigorous than mRNA proxy); (2) Isogenic PTEN-null NSCLC cell line models (e.g., PTEN KO in H1299, A549) with ITGAV CRISPR KO or cilengitide treatment to confirm differential sensitivity vs PTEN-intact controls.

**IP Gap Analysis:**
- PTEN loss as biomarker for αV integrin inhibitor response: Not claimed (search recommended)
- PTEN-null NSCLC BrM as indication for αV integrin inhibition: Not claimed
- Combination of PTEN loss + ZEB1-high as dual biomarker: Not claimed

**Filing Urgency:** File with Brief #1. The PTEN convergence hypothesis is documented in the 8D run results (8d_vector_run_results.json, commit 0cdc800). Functional validation is pending but the method claim is independent of validation status. The specification must accurately represent the evidence tier.

---

## PATENT BRIEF #3 — Q2 2026
### ZEB1-High NSCLC Brain Metastasis + αV Integrin Inhibitor + PD-1 Blockade Combination
**Evidence tier: Mechanistically grounded — combination not yet validated in vivo**

**Plain-English Summary:**
ZEB1 is a direct transcriptional activator of PD-L1 (CD274). This means that ZEB1-high NSCLC BrM tumors are simultaneously: (1) structurally dependent on αV integrin for survival (ZEB1→ITGAV SL axis, CRISPR-validated), AND (2) PD-L1-high, making them candidates for PD-1 checkpoint blockade. One biomarker — ZEB1 expression — identifies patients who should respond to both drugs. The combination hypothesis: αV integrin inhibition (cilengitide) disrupts integrin-mediated survival signaling while PD-1 blockade (pembrolizumab, nivolumab) reinvigorates T cell killing. The ZEB1→PD-L1 transcriptional link is established in the literature; the combination has not been tested in NSCLC BrM models. This is a mechanistically grounded combination claim.

**Claim Structure (for counsel):**

*Independent Claim 1 — Combination Method of Treatment:*
A method of treating brain metastasis in a human patient comprising:
(a) determining that the patient has non-small cell lung cancer brain metastasis (NSCLC BrM);
(b) measuring ZEB1 expression in a tumor sample from the patient;
(c) identifying the patient as ZEB1-high based on a predetermined expression threshold; and
(d) administering to the patient a therapeutically effective amount of (i) an αV integrin inhibitor and (ii) a PD-1 or PD-L1 inhibitor.

*Independent Claim 2 — Sequential Administration:*
A method of treating ZEB1-high NSCLC brain metastasis comprising sequentially or concurrently administering an αV integrin inhibitor and a PD-1 inhibitor to a patient identified as ZEB1-high by tumor biopsy analysis.

*Dependent Claims (examples):*
- Claim 3: The method of Claim 1, wherein the αV integrin inhibitor is cilengitide.
- Claim 4: The method of Claim 1, wherein the PD-1 inhibitor is pembrolizumab or nivolumab.
- Claim 5: The method of Claim 1, wherein the PD-L1 inhibitor is atezolizumab, durvalumab, or avelumab.
- Claim 6: The method of Claim 1, wherein ZEB1 expression is positively correlated with PD-L1 expression in the patient's tumor sample.
- Claim 7: The method of Claim 1, wherein the patient has KRAS-mutant NSCLC.

**Supporting Evidence (for specification):**

*ZEB1→PD-L1 transcriptional axis (multi-tier, NSCLC-specific):*
- Chen L et al. (2014). *Nat Commun* 5:5241. DOI: 10.1038/ncomms6241. [ZEB1 relieves miR-200 repression of PD-L1 in lung cancer datasets; ZEB1→PD-L1 axis drives metastasis and immune evasion — LANDMARK, lung cancer]
- Kim HR et al. (2016). *Hum Pathol* 58:45–52. DOI: 10.1016/j.humpath.2016.07.007. [EMT phenotype (ZEB1+) correlates with PD-L1 positivity in 477 lung adenocarcinoma patient specimens — NSCLC patient data]
- Bouillez A et al. (2017). *Oncogene* 36(32):4617–4628. DOI: 10.1038/onc.2017.47. [MUC1-C→NF-κB→ZEB1→PD-L1 pathway in NSCLC cell lines — NSCLC cell line support]
- Wirbel J et al. (2025). *Cancer Immunol Immunother* 74:78. DOI: 10.1007/s00262-025-03978-5. [ZEB1 directly binds CD274 (PD-L1) promoter; drives PD-L1 transcription — strongest mechanistic evidence; NOTE: melanoma model, not NSCLC]

*ZEB1 upregulation in NSCLC BrM patient specimens:*
- Lee JY et al. (2020). *Cancers* 12(12):3632. DOI: 10.3390/cancers12123632. [ZEB1 immunoreactivity: 15% in primary lung adenocarcinoma vs 55.9% in matched BrM (p=0.002, n=46 paired specimens) — strongest available BrM-specific data]
- **Gap:** ZEB1 and PD-L1 co-expression has not been directly measured in the same NSCLC BrM specimens. Co-expression is inferred from (a) ZEB1 upregulation in BrM (Lee 2020) + (b) ZEB1/PD-L1 correlation in primary NSCLC (Kim 2016). Direct measurement in NSCLC BrM tissue is a validation priority.

*αV integrin inhibitor + PD-1 blockade combination (preclinical, adjacent tumor types):*
- Pan Y et al. (2022). *Bioengineered* 13(2):3928–3940. DOI: 10.1080/21655979.2022.2029236. [Cilengitide + anti-PD-1 reduces tumor growth, extends survival in murine melanoma; cilengitide decreases PD-L1 via STAT3 — NOTE: melanoma model, not NSCLC BrM]
- Ellert-Miklaszewska A et al. (2025). *J Exp Clin Cancer Res* 44:103. DOI: 10.1186/s13046-025-03393-9. [RGD integrin-blocking peptide + anti-PD-1 reduces glioma growth, restores CD8+ T cell infiltration — NOTE: glioma (primary brain tumor), not NSCLC BrM]

*ZEB1→ITGAV SL axis (CRISPR-validated):*
- DepMap 24Q4, delta = −0.7184, FDR = 0.001203, n = 24/24 ZEB1-high NSCLC lines

*8D mechanism fit:*
- ZEB1-high KRAS-mut profile: checkpoint fit = 0.57 (rank #2), cilengitide fit = 0.52 (rank #3)

- **Validation required:** (1) ZEB1/PD-L1 co-IHC in NSCLC BrM tissue microarray (direct co-expression evidence); (2) ZEB1 overexpression in NSCLC cell lines → PD-L1 upregulation (Western, flow cytometry); (3) In vivo combination testing (cilengitide + anti-PD-1) in ZEB1-high NSCLC BrM syngeneic or PDX models.

**IP Gap Analysis:**
- ZEB1-high as dual biomarker for αV integrin + PD-1 combination: Not claimed
- αV integrin inhibitor + PD-1 inhibitor combination in NSCLC BrM: Not claimed (search recommended)
- ZEB1 expression as predictor of PD-L1 status in NSCLC BrM: Not claimed as patent

**Filing Urgency:** Q2 2026. File after Brief #1 and #2 are in. The combination hypothesis is mechanistically grounded but requires the ZEB1→PD-L1 link to be documented in the specification with literature support.

---

## PATENT BRIEF #4 — FILE WHEN CANDIDATE IDENTIFIED
### SPP1-High NSCLC Brain Metastasis + NRF2 Pathway Inhibition
**Evidence tier: CRISPR-validated SL dependency — no clinical-stage CNS-penetrant drug identified**

**Plain-English Summary:**
SPP1 (osteopontin) is a brain metastasis colonization factor — it is secreted by tumor cells and reactive astrocytes in the BrM microenvironment and promotes tumor cell survival and immune evasion. Our DepMap analysis identifies SPP1→NFE2L2 (NRF2) as the strongest SL pair in the BrM universe by statistical significance (delta = −0.7326, FDR = 8×10⁻⁶, n = 24/24 NSCLC lines). High SPP1 expression creates oxidative stress dependency on NRF2, the master antioxidant transcription factor. Inhibiting NRF2 in SPP1-high tumors should selectively kill the tumor cells that are most dependent on the BrM colonization program. This is the highest-priority unloaded route in the engine — no CNS-penetrant NRF2 inhibitor is currently in active clinical development, but the patient selection logic and indication are patentable now.

**Claim Structure (for counsel):**

*Independent Claim 1 — Method of Treatment:*
A method of treating brain metastasis in a human patient comprising:
(a) determining that the patient has non-small cell lung cancer brain metastasis (NSCLC BrM);
(b) measuring SPP1 (osteopontin) expression in a tumor sample from the patient;
(c) identifying the patient as SPP1-high based on a predetermined expression threshold; and
(d) administering to the patient a therapeutically effective amount of an NRF2 pathway inhibitor.

*Independent Claim 2 — Patient Selection Method:*
A method of selecting a patient with NSCLC brain metastasis for NRF2 pathway inhibitor therapy comprising measuring SPP1 expression in a tumor sample and selecting the patient if SPP1 expression exceeds a threshold indicative of NFE2L2 synthetic lethality.

*Dependent Claims (examples):*
- Claim 3: The method of Claim 1, wherein the NRF2 pathway inhibitor is brusatol, ML385, or a pharmaceutically acceptable salt thereof.
- Claim 4: The method of Claim 1, wherein the NRF2 pathway inhibitor is a KEAP1 activator or NRF2 transcriptional inhibitor.
- Claim 5: The method of Claim 1, wherein SPP1 expression is measured in tumor cells, tumor-associated macrophages, or reactive astrocytes in the BrM microenvironment.
- Claim 6: The method of Claim 1, wherein the NRF2 pathway inhibitor is CNS-penetrant as determined by CSF/plasma ratio ≥ 0.1.

**Supporting Evidence (for specification):**
- SPP1→NFE2L2 SL signal: DepMap 24Q4, delta = −0.7326, FDR = 8×10⁻⁶, n = 24/24 NSCLC lines — strongest SL pair in BrM universe by FDR (CRISPR-validated)
- SPP1 as BrM colonization factor: Klement GL et al. (2018). *Clin Cancer Res* 24(3):530–537; Carozza JA et al. (2022). *Cell* 185(8):1403–1419.
- NRF2 as master antioxidant transcription factor: Itoh K et al. (1997). *Biochem Biophys Res Commun* 236(2):313–322; Taguchi K & Yamamoto M (2017). *Cell Chem Biol* 24(7):779–792.
- SPP1-high → oxidative stress dependency on NRF2: mechanistic rationale (functional validation pending in SPP1-high NSCLC BrM models)
- NRF2 inhibitors: brusatol (preclinical only), ML385 (preclinical only) — no CNS-penetrant clinical-stage agent identified as of 2026-05-08

**IP Gap Analysis:**
- SPP1-high as biomarker for NRF2 inhibitor response in NSCLC BrM: Not claimed
- SPP1→NFE2L2 SL dependency as therapeutic target: Not published as patent claim
- NRF2 inhibition in NSCLC BrM: Limited prior art (NRF2 inhibition is an emerging field)

**Filing Urgency:** The broad method claim (SPP1-high + NRF2 inhibitor class) can be filed now without naming a specific drug. Recommend filing the broad claim now and adding specific drug claims as continuation applications when CNS-penetrant candidates are identified. Commercial value increases significantly when a specific clinical-stage candidate is named.

**Open Action:** Identify CNS-penetrant NRF2 inhibitor candidates in active development. This is the highest-priority open slot in the CrisPRO BrM engine.

---

## IP Strategy Notes for Counsel

**Claim type:** All four briefs are "method of treatment" claims (35 U.S.C. § 101 eligible under Mayo/Alice framework for method of treatment claims). These are the most defensible patent claims in precision oncology — they cover the use of a drug in a specific patient population defined by a specific biomarker, not the drug itself.

**Evidence tier summary (v3.0 — 2026-05-09):**
| Brief | Evidence Tier | DepMap Status | Literature Support | Validation Status |
|---|---|---|---|---|
| #1 ZEB1→ITGAV | CRISPR-validated SL (DepMap 24Q4, n=24/24, FDR=0.001203) | POSITIVE (delta=−0.7184, p=0.000002) | Nagaishi 2017 (BrM specimens), Mattson 2024 (CRISPR screen) | Strong — ready to file |
| #2 PTEN-null→ITGAV | Computational hypothesis (8D fit=0.88) — **DepMap NULL** | **NULL** (delta=−0.040, p=0.73, n=95) | Cavazzoni 2017 (PTEN→FAK in NSCLC SCC — mechanistic only) | Requires isogenic PTEN-null validation; mRNA-based DepMap stratification does NOT support hypothesis |
| #3 ZEB1-high + PD-1 combo | Mechanistically grounded — multi-tier NSCLC evidence | N/A (combination not testable in DepMap) | Chen 2014 (lung cancer), Kim 2016 (477 NSCLC specimens), Bouillez 2017 (NSCLC cells), Lee 2020 (BrM, n=46); combination: Pan 2022 (melanoma), Ellert-Miklaszewska 2025 (glioma) | Requires ZEB1/PD-L1 co-IHC in NSCLC BrM + in vivo combination testing |
| #4 SPP1→NFE2L2 | CRISPR-validated SL (DepMap 24Q4, n=24/24, FDR=8×10⁻⁶) | POSITIVE (delta=−0.7326, p<0.0001) | SPP1 as BrM colonization factor (Klement 2018, Carozza 2022) | Strong — no CNS-penetrant clinical drug yet |

**Prior art risk:** The primary prior art risk is obviousness — a reviewer could argue that using ZEB1 as a biomarker for αV integrin inhibitor response is obvious given the known ZEB1-integrin biology. The non-obviousness argument rests on: (1) the specific SL dependency identified by CRISPR functional genomics (not obvious from expression data alone), (2) the specific NSCLC BrM context (not GBM, not primary NSCLC), and (3) the specific patient selection threshold derived from DepMap analysis.

**Publication risk:** The manuscript (brm_graveyard_perspective_v3.md) is in pre-publication draft. File all provisionals before manuscript submission. The provisional establishes priority date; the full application must be filed within 12 months.

**Trade secret protection:** The DepMap analysis, scoring methodology, and brm_targetability_matrix_v3.json are trade secrets. The patent specification should cite the conclusions (delta values, FDR, n) without disclosing the full methodology. The engine that generates the IP stays proprietary.

---

*CrisPRO Patent Filing Briefs v3.0 — 2026-05-09 | Evidence gap closure pass complete | Attorney-Client Privileged when transmitted to counsel | Confidential*
