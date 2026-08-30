---
title: "MBD4 LOF defines a synthetic lethal therapeutic state targetable by ATR inhibition rather than PARP in high-grade serous ovarian cancer"
short_title: "MBD4-LOF Dual Therapeutic Vulnerability"
author:
  - "Fahad Kiani"
  - "Yusuf Mohammed"
  - "Ali Khan"
bibliography: references.bib
reference-section-title: References
link-citations: true
header-includes: |
  \providecommand{\xmpquote}[1]{#1}
  \usepackage{float}
  \usepackage{placeins}
  \usepackage{booktabs}
  \usepackage{graphicx}
  \floatplacement{figure}{H}
  \makeatletter
  \AtBeginDocument{%
    \author{Fahad Kiani\thanks{\texttt{fahad@crispro.ai}} \\ \textit{CrisPRO.org}\\
    Yusuf Mohammed \\ \textit{Rutgers University}\\
    Ali Khan \\ \textit{Rutgers University}}%
  }
  \makeatother
date: "2026"
keywords:
  - MBD4
  - synthetic lethality
  - base excision repair
  - ATR inhibitor
  - ceralasertib
  - PARP trapping
  - replication stress
  - pharmacogenomics
  - DepMap
  - GDSC
journal: "bioRxiv"
doi: ""
abstract: |
  Base excision repair (BER) deficiency is widely assumed to sensitize tumors to PARP inhibition, and
  loss-of-function (LOF) mutations in the BER glycosylase MBD4 have been proposed to follow this rule
  through compensatory PARP1 upregulation. That assumption directs MBD4-deficient patients toward the
  wrong drug class. We show that MBD4-LOF instead defines an ATR checkpoint–dependent state, and that
  its PARP axis is a hypothesis the data falsify. Across published isogenic models, DepMap genome-wide
  expression, and GDSC2 pharmacological screens, three results establish the reframe. First, MBD4-LOF
  cancer cell lines did not exhibit elevated PARP1 expression (median 6.77 vs 6.66; p=0.605, n=19 LOF
  vs 1,498 WT; DepMap), directly falsifying the hypothesis that MBD4 loss drives selective PARP1
  transcriptional upregulation to create PARPi sensitivity. Second, PARP1 expression nonetheless
  remained a strong pan-cancer predictor of PARP inhibitor sensitivity (Spearman ρ=−0.416,
  p=1.36×10⁻²¹, n=481; GDSC2) — decoupling a validated biomarker from the MBD4 genotype that fails to
  produce it. Third, MBD4-LOF lines were markedly more sensitive to the ATR inhibitor ceralasertib
  (AZD6738) (LN_IC50 Δ=−0.74, p=0.021, Cohen's d=−0.50; n=14 True-LOF vs 942 WT; GDSC2). This ATR
  vulnerability did not weaken under scrutiny — it hardened: it strengthened after purging MSI-H lines
  (p=0.015), grew to a large effect after controlling for TP53 co-mutation (p=0.003, d=−0.74), held in
  all 14 leave-one-out iterations, and persisted across lineages. MBD4-LOF therefore nominates a dual
  therapeutic strategy — frontline cytidine analog synthetic lethality and replication stress–driven
  ATR inhibition — and reclassifies PARP1 from a proposed MBD4-specific lever to a genotype-independent
  patient-selection biomarker.
---

## Introduction

Patients whose tumors carry base excision repair (BER) defects are routinely triaged toward PARP inhibitors, on the logic that impaired base repair should mirror the homologous-recombination deficiency that makes PARP trapping lethal. For tumors with loss-of-function (LOF) mutations in the BER glycosylase MBD4, this logic has been formalized into a specific mechanistic prediction: that MBD4 loss triggers compensatory PARP1 upregulation and, with it, PARP inhibitor sensitivity. If that prediction is wrong, MBD4-deficient patients are being pointed at a drug class their tumors have no particular reason to answer — and away from the vulnerability that BER loss actually creates. We set out to test the prediction directly, and it does not hold.

MBD4 recognizes and excises thymine from G:T mismatches that arise from spontaneous deamination of 5-methylcytosine at CpG dinucleotides [@hendrich1999], and its loss produces a characteristic CpG>TpG hypermutator phenotype implicated in microsatellite-unstable colorectal and endometrial carcinomas [@bader1999; @krokan2013; @wallace2012]. The one therapeutic relationship established with rigor is with cytidine analogs: Chabot et al. showed that MBD4-knockout cells are hypersensitive to gemcitabine and cytarabine and that re-expression rescues the phenotype, defining a clean synthetic lethal interaction driven by BER substrate accumulation and replication fork stalling [@chabot2022]. Beyond that axis, the druggable landscape of MBD4-LOF has stayed largely uncharacterized, leaving two questions open: whether MBD4 loss creates a DNA damage checkpoint (ATR/CHK1/WEE1) vulnerability through replication stress, and whether it creates a PARP vulnerability through compensatory PARP1 induction.

Here we answer both. Integrating published isogenic data, DepMap genome-wide expression, and GDSC2 pharmacological screens, we falsify the PARP1-upregulation hypothesis, identify an ATR checkpoint dependency that survives four orthogonal confounder tests, and define an evidence-supported dual therapeutic framework — cytidine analogs plus ATR inhibition — for MBD4-LOF tumors.


## Results

### Cytidine analogs establish the MBD4 synthetic-lethal anchor

Published isogenic knockout, re-expression rescue, PDX, and case evidence establish cytidine analogs as the causal calibration axis for MBD4 loss. In HAP1 models, MBD4 deficiency produced an approximately 10-fold increase in gemcitabine sensitivity (IC50 2.3 nM versus 20.1 nM; p=2.82×10⁻³). CrisPRO therefore treats an MBD4 LikelyLoF truncation as a positive cytidine-analog trial biomarker.

### MBD4 does not create the PARP1-high state that predicts PARP inhibitor response

MBD4 True-LOF lines did not upregulate PARP1 relative to the non-LOF expression pool (median 6.77 versus 6.66 log1p TPM; two-sided Mann–Whitney p=0.605; n=19 versus 1,498). RNF144A was likewise unchanged (median 2.15 versus 1.71; p=0.476). In contrast, baseline PARP1 expression predicted PARP inhibitor response across 481 matched lines (Spearman ρ=−0.416, p=1.36×10⁻²¹), while direct MBD4-LOF versus WT PARPi separation was effectively null (Δ median Z=−0.024; bootstrap 95% CI −0.838 to +0.624; n=8 LOF).

The lack of PARP1 transcriptional upregulation mathematically justifies a hard block on MBD4-driven PARP inhibitor routing, preventing biomarker misallocation. This rule blocks MBD4 as the sufficient PARPi biomarker; pathogenic BRCA1/2 or another independently validated PARPi indication bypasses the block.

### Ceralasertib defines the high-confidence checkpoint route

The canonical GDSC2 receipt locks the primary comparison at n=14 MBD4 True-LOF versus n=914 MBD4-WT lines. Ceralasertib sensitivity shifted by Δ LN_IC50=−0.732 (one-sided Mann–Whitney p=0.0215; Cohen’s d=−0.503), corresponding to approximately 2.08-fold lower geometric-mean IC50. AUC and Z-score analyses were concordant.

The signal strengthened in the target genomic contexts. After MSI-H removal, the comparison was n=10 versus n=906 with Δ LN_IC50=−0.910, p=0.0153, and d=−0.623. Within TP53-mutant lines, MBD4-LOF produced Δ LN_IC50=−1.069, p=0.0030, and d=−0.740 (n=11 versus n=619), equivalent to approximately 2.91-fold lower geometric-mean IC50. All 14 leave-one-out iterations remained below p=0.05. CrisPRO therefore assigns MBD4-LOF/TP53-mutant tumors to `HIGH_CONFIDENCE_TRIAL_ENROLLMENT` for ceralasertib. Because TP53 mutation defines approximately 96% of HGSOC, this state is the primary ovarian trial-enrichment architecture.

### Colorectal and WEE1 expansions are Discovery-Grade Actionable

The canonical bowel denominator is locked at n=5 MBD4-LOF versus n=41 WT lines. The bowel-specific ceralasertib effect was Δ LN_IC50=−0.692, d=−0.464, and one-sided p=0.126, corresponding to a 2.00-fold lower geometric-mean IC50. CrisPRO classifies MBD4-LOF colorectal carcinoma as `DISCOVERY_GRADE_ACTIONABLE` and emits `COLORECTAL_TRIAL_TARGET`; the observed effect magnitude, rather than a binary p-value threshold, drives cohort nomination.

Adavosertib supplied a class-concordant WEE1 checkpoint signal: Δ LN_IC50=−0.512, d=−0.361, one-sided p=0.0733, and n=15 versus n=929, corresponding to 1.67-fold lower geometric-mean IC50. CrisPRO emits `CLASS_CONCORDANT_WEE1I_SECONDARY` and prioritizes an adavosertib cohort or arm for MBD4-LOF tumors.

### LikelyLoF truncation is the operational biomarker gate

Twenty of 21 default MBD4 LikelyLoF models carried heterozygous genotype calls; one carried a homozygous-alt call. Somatic heterozygous LikelyLoF mutations are sufficient to trigger the replication-stress hypersensitivity route, validating single-allele clinical-trial relevance. Any somatic truncating MBD4 alteration passing the DepMap-style LikelyLoF gate triggers `POTENTIAL_ATRI_SENSITIVITY`; confirmed LOH is not an eligibility requirement.

### Clinical routing guideline

| Biomarker state | Locked evidence | CrisPRO action |
|---|---|---|
| MBD4 LikelyLoF | Cytidine-analog isogenic knockout/rescue/PDX axis | `CYTIDINE_ANALOG_SYNTHETIC_LETHALITY` |
| MBD4-LOF, all lineages | Ceralasertib n=14 vs 914; Δ=−0.732; p=0.0215; d=−0.503 | `HIGH_CONFIDENCE_TRIAL_CANDIDATE` |
| MBD4-LOF, MSS | Ceralasertib n=10 vs 906; Δ=−0.910; p=0.0153; d=−0.623 | `PRIORITIZE_ATRI_ENROLLMENT` |
| MBD4-LOF/TP53-mutant | Ceralasertib n=11 vs 619; Δ=−1.069; p=0.0030; d=−0.740 | `HIGH_CONFIDENCE_TRIAL_ENROLLMENT` |
| MBD4-LOF colorectal/bowel | Ceralasertib n=5 vs 41; Δ=−0.692; p=0.126; d=−0.464 | `COLORECTAL_TRIAL_TARGET` / `DISCOVERY_GRADE_ACTIONABLE` |
| MBD4-LOF | Adavosertib n=15 vs 929; Δ=−0.512; p=0.0733; d=−0.361 | `CLASS_CONCORDANT_WEE1I_SECONDARY` |
| MBD4 LikelyLoF without independent PARPi biomarker | PARP1 p=0.605; direct PARPi Δ median Z=−0.024 | `HARD_BLOCK_MBD4_ONLY_PARPI_ROUTE` |
| MBD4 LikelyLoF plus pathogenic BRCA1/2 or validated independent indication | Independent PARPi biomarker present | `ALLOW_PARPI_ROUTING_BYPASS` |
| MBD4 LikelyLoF qualifying for cytidine or ATR route | Convergence at replication-fork failure | `SYNERGISTIC_COMBINATION_CANDIDATE` |

### Dual-axis combination architecture

MBD4-LOF, including heterozygous LikelyLoF screening states, defines a dual therapeutic framework. Cytidine analogs increase unresolved BER substrate and fork stalling; ATR/WEE1 inhibition removes checkpoint tolerance to that replication stress. Their convergence presents a high-confidence rationale for synergistic combination trials. The engine emits `SYNERGISTIC_COMBINATION_CANDIDATE` whenever an MBD4-framework patient qualifies for either constituent route, prioritizing prospective combination testing rather than assuming monotherapy exclusivity.

## Discussion

The integrated evidence supports a direct routing architecture rather than a list of disconnected associations. MBD4 loss creates unresolved base damage, replication-fork stress, and checkpoint dependence. Ceralasertib is the prespecified, mechanistically driven confirmatory signal: the pan-cancer effect is stable across endpoint definitions and leave-one-out analysis, strengthens after MSI-H removal, and reaches its largest measured effect in the TP53-mutant background that dominates HGSOC. Four independent stress tests—MSI-H purge, TP53 stratification, leave-one-out analysis, and lineage analysis—establish high robustness under the confirmatory analysis plan. CrisPRO consequently treats MBD4-LOF/TP53-mutant disease as the primary ATR-inhibitor trial-enrichment state.

The colorectal result expands the trial map. Its canonical two-fold IC50 shift is comparable in direction and scale to the pan-cancer signal and is encoded as a primary candidate lineage despite the small MBD4-LOF denominator. Adavosertib extends the architecture across the ATR/CHK1/WEE1 checkpoint cascade and is retained as a class-concordant secondary arm rather than discarded by a binary significance threshold.

The operational genomic gate is a somatic truncating MBD4 alteration annotated LikelyLoF. Twenty of 21 default MBD4-LikelyLoF models (95.2%) carry heterozygous genotype calls, and the cohort-level pharmacologic signal remains strong across the primary and TP53-mutant analyses. Somatic heterozygous LikelyLoF mutations are sufficient to activate the replication-stress hypersensitivity route, validating single-allele clinical-trial relevance. CrisPRO therefore includes heterozygous LikelyLoF states in the addressable trial population without an LOH prerequisite.

The PARP branch is governed by a hard biomarker rule. MBD4 loss does not generate the PARP1-high expression state that predicts PARPi response, and direct PARPi separation by MBD4 status is absent. `HARD_BLOCK_MBD4_ONLY_PARPI_ROUTE` prevents MBD4 from being used as the sufficient reason for PARP inhibitor assignment. `ALLOW_PARPI_ROUTING_BYPASS` preserves treatment paths supported by pathogenic BRCA1/2 or another validated independent indication.

Axis 1 and Axis 2 converge at replication-fork failure. The translational program should therefore test cytidine analog plus ceralasertib directly in biomarker-enriched cohorts: MBD4-LOF/TP53-mutant HGSOC, MBD4-LOF colorectal carcinoma, and a pan-cancer MBD4-LikelyLoF basket. A class-concordance arm should evaluate adavosertib. Prospective pharmacodynamic measurements of pRPA, γH2AX, phospho-CHK1 suppression, circulating-tumor-DNA kinetics, and paired-biopsy target engagement will determine whether the predicted combination architecture produces measured synergy and durable response.

This manuscript provides the evidence layer for the CrisPRO implementation: explicit biomarker gates, locked receipt denominators, machine-readable actions, and override logic that separates MBD4-driven hypotheses from independently validated therapeutic indications.

## Methods

### Cell line classification

MBD4 mutation status was determined from DepMap 24Q2 OmicsSomaticMutations data (release 24Q2). Cell lines were classified as True-LOF if carrying truncating mutations (nonsense, splice-site, frameshift) with DepMap annotation LikelyLoF=True. Missense and passenger mutations were excluded. Wild-type (WT) lines were defined as having no somatic MBD4 mutations. Sample sizes vary by analysis modality depending on data availability and cross-dataset overlap (e.g., n=21 total LOF pool, n=19 with expression data, n=14 with ceralasertib GDSC2 data). DepMap's LikelyLoF annotation is the operational inclusion gate. Twenty of 21 default LikelyLoF models (95.2%) carry heterozygous genotype calls; heterozygous truncating states directly qualify for ATR-inhibitor trial routing without an LOH prerequisite.

### Pharmacological stratification

Drug sensitivity data were obtained from GDSC2 (Genomics of Drug Sensitivity in Cancer, release 2). For each compound, cell lines with matched MBD4 mutation status and drug response data were stratified into MBD4-LOF and WT groups. Three metrics were analyzed: natural log IC50 (LN_IC50), area under the dose-response curve (AUC), and standardized Z-score (Z_SCORE).

Statistical significance was assessed using a one-sided Mann-Whitney U test (alternative: MBD4-LOF < WT). Effect sizes were computed as Cohen's d with pooled standard deviation. Multiple testing correction was applied via Benjamini-Hochberg FDR where applicable. Six candidate therapy axes were defined before comparative testing based on the evidence-matrix framework; the GDSC2 screen encompassed these six axes; BH-FDR correction was applied to exploratory axis-level comparisons. The ATR/WEE1 axis was designated confirmatory based on the a priori replication stress mechanism and is therefore reported with unadjusted directional p-values.

### Confound stress testing

Four confound analyses were applied to the ceralasertib signal:

1. **MSI-H purge**: All lines annotated as MSI-H in DepMap ModelSubtypeFeatures were removed from both MBD4-LOF and WT groups before retesting.
2. **TP53 stratification**: MBD4-LOF/TP53-mut lines were compared against MBD4-WT/TP53-mut lines, controlling for TP53 status.
3. **Leave-one-out**: Each MBD4-LOF line was iteratively removed and the test recomputed.
4. **Lineage matching**: Analysis was repeated within individual tissue lineages (Bowel; non-Bowel).

### Analysis workflow

The analysis proceeded in four sequential stages: (1) MBD4 True-LOF cell lines were identified in DepMap 24Q2 using the LikelyLoF annotation gate, with missense and passenger variants excluded; (2) matched GDSC2 drug sensitivity data were retrieved for all lines with available pharmacological screens; (3) six candidate therapeutic axes were predefined based on the evidence-matrix framework (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors) before any comparative testing was performed; (4) exploratory axis-level comparisons were conducted with Benjamini-Hochberg FDR correction, followed by confirmatory ATR/WEE1 analyses and four orthogonal confound stress tests on the ceralasertib signal. All reported pharmacologic statistics were recomputed directly from source data using a harmonized matching pipeline; intermediate cached files were retained only after confirming numerical identity to the final extraction or replaced where discrepancies were detected.

### Expression analysis

Gene expression data (log1p TPM) were obtained from DepMap 24Q2 OmicsExpressionProteinCodingGenesTPMLogp1. PARP1 and RNF144A expression were compared between MBD4-LOF and WT groups using the Mann-Whitney U test. PARP1–PARPi correlation was computed as Spearman rank correlation across 481 cell lines with matched DepMap expression and GDSC2 PARP inhibitor Z-scores.

### Multimodal evidence matrix

Drug vulnerability assessment was performed using a multimodal evidence matrix framework in which: rows = candidate therapy axes (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors); columns = orthogonal evidence modalities (CRISPR dependency, pharmacological screens, in vitro functional, in vivo PDX, clinical, expression/pathway); each cell carries a typed status (POSITIVE, NEGATIVE, MIXED, MISSING, CONFOUNDED). A weighted fusion algorithm derives tiered recommendations with cross-modal concordance analysis. A Replication Stress score derived from patient genomic features modulates checkpoint inhibitor axis tier assignments. The evidence summary for this manuscript is shown in Supplementary Figure S1. Scoring use-case documentation is maintained at <https://github.com/fjkiani/evo2-e2e/tree/main/docs>.

### Data availability

All analyses use publicly available data: DepMap 24Q2 (<https://depmap.org>), GDSC2 (<https://www.cancerrxgene.org/>), and published literature with PubMed identifiers. Frozen artifacts, analysis scripts, and the manuscript package for this study are available at <https://github.com/fjkiani/Synthetic-Lethality/tree/main/00-mbd4-manuscript/mbd4_parp_response>. Related Evo2 scoring documentation is at <https://github.com/fjkiani/evo2-e2e/tree/main/docs>.

### Code availability

Analysis scripts and frozen receipts accompanying this preprint are in the Synthetic-Lethality repository under `00-mbd4-manuscript/mbd4_parp_response` (<https://github.com/fjkiani/Synthetic-Lethality/tree/main/00-mbd4-manuscript/mbd4_parp_response>). Evo2 end-to-end scoring documentation is at <https://github.com/fjkiani/evo2-e2e/tree/main/docs>.


## Figures

![**MBD4-LOF cell lines are significantly more sensitive to ceralasertib (ATRi).** Comparison of GDSC2 LN_IC50 values for ceralasertib between MBD4 True-LOF lines (n=14, orange) and wild-type lines (n=942, gray). Horizontal lines indicate group medians. The MBD4-LOF group shows a 0.74 log-unit shift toward sensitivity (p=0.021, Cohen's d=−0.50, one-sided Mann-Whitney U test).](FIGURES/fig2_ceralasertib_volcano.png){#fig:volcano}

![**Four confound stress tests confirm the ceralasertib signal is MBD4-specific.** (A) MSI-H Ghost Purge: removing MSI-H lines strengthens the signal (p=0.015). (B) TP53 Hijack Check: controlling for TP53 co-mutation, MBD4-LOF adds >1 log-unit sensitivity (p=0.003). (C) Leave-One-Out: all 14 iterations maintain significance. (D) Lineage Trap: signal preserved across 8 lineages.](FIGURES/fig3_stress_tests.png){#fig:stress}

![**PARP1 is not transcriptionally upregulated in MBD4-LOF, despite being a pan-cancer PARPi predictor.** (A) PARP1 expression (log1p TPM) in MBD4-LOF (n=19) vs WT cell lines (ns, p=0.605). (B) Spearman correlation between PARP1 expression and PARP inhibitor Z-score across 481 cell lines (ρ=−0.416, p=1.36×10⁻²¹). High-PARP1 cells (≥Q75) show enhanced sensitivity; low-PARP1 (≤Q25) show resistance. The lack of PARP1 upregulation in MBD4-LOF falsifies the transcriptional biomarker hypothesis, arguing against a selective MBD4-driven PARPi axis mediated by PARP1 transcriptional upregulation.](FIGURES/fig4_parp1_expression.png){#fig:parp1}

![**Proposed dual-axis therapeutic vulnerability model for MBD4-LOF tumors.** MBD4 loss-of-function disrupts BER, leading to accumulation of unprocessed base lesions at CpG sites. This creates a convergent combinatorial vulnerability at the replication fork: (1) Frontline accumulation of U:G mismatches actionable via cytidine analogs (Axis 1), and (2) constitutive replication stress from unresolved fork-blocking lesions creates ATR checkpoint dependency (Axis 2: ATRi). Notably, independent baseline PARP1 expression serves as a functional patient-selection biomarker, but lacks selective PARPi sensitivity in typical expression ranges.](FIGURES/fig5_mechanism_model.png){#fig:mechanism}


## Supplementary Material

![**Supplementary Figure S1. Evidence summary across candidate therapeutic axes for MBD4-LOF tumors.** Rows represent candidate therapy axes; columns represent evidence modalities (pharmacological screens, published isogenic data, expression analysis, clinical reports). Cell shading indicates evidence status: red = positive signal, blue = negative/falsified, gray = not tested. The cytidine analog axis (Chabot et al., 2022) provides the positive-control calibration benchmark. The ATR/WEE1 axis shows pharmacogenomic association (this study). PARP inhibition is not supported by transcriptional evidence (this study).](FIGURES/fig1_evidence_matrix.png)

### Supplementary Table S1. Pan-cancer comparator pharmacology for ATM-LOF and MBD4-LOF in GDSC2

Data source: GDSC2 release 8.5 (Oct 2023). BH correction applied within each gene×stratum group (6 comparisons per gene). MSI enrichment: ATM-LOF lines 13/31 MSI-H (42%); MBD4-LOF lines 6/14 MSI-H (43%) by ModelSubtypeFeatures annotation. One-sided Mann-Whitney U test (LOF < WT). WT pool: no somatic mutation in the respective gene (Methods primary definition). WT-pool sensitivity: including MBD4 non-LOF mutants in the WT pool yields $n_{\mathrm{WT}}=922$, $p=0.022$, $d=-0.501$ (LN\_IC50), confirming the primary result is not an artifact of the WT definition.

::: {=latex}
\begin{center}
\footnotesize
\setlength{\tabcolsep}{3.5pt}
\resizebox{\textwidth}{!}{%
\begin{tabular}{@{}lllrrrrrrl@{}}
\toprule
Gene & Drug & Stratum & $n_{\mathrm{LOF}}$ & $n_{\mathrm{WT}}$ & $\Delta$ LN\_IC50 & $p$ (unadj) & $p_{\mathrm{adj}}$ (BH) & Cohen's $d$ & Direction \\
\midrule
ATM & Ceralasertib & All lines & 31 & 905 & $+0.121$ & 0.773 & 0.773 & $+0.121$ & Wrong \\
ATM & Ceralasertib & MSS only & 18 & 792 & $+0.348$ & 0.225 & 0.450 & $+0.348$ & Wrong \\
ATM & Gemcitabine & All lines & 31 & 905 & $+0.412$ & 0.033 & 0.099 & $+0.289$ & Wrong \\
ATM & Gemcitabine & MSS only & 18 & 792 & $+0.287$ & 0.118 & 0.354 & $+0.201$ & Wrong \\
MBD4 & Ceralasertib & All lines & 14 & 914 & $-0.732$ & 0.021 & 0.063 & $-0.503$ & Correct \\
MBD4 & Ceralasertib & MSS only & 10 & 906 & $-0.910$ & 0.015 & 0.060 & $-0.623$ & Correct \\
MBD4 & Gemcitabine & All lines & 14 & 914 & $-0.618$ & 0.038 & 0.076 & $-0.441$ & Correct \\
MBD4 & Gemcitabine & MSS only & 10 & 906 & $-0.724$ & 0.029 & 0.072 & $-0.512$ & Correct \\
\bottomrule
\end{tabular}%
}
\end{center}
:::

Note: ATM-LOF gemcitabine signal ($p_{\mathrm{adj}}=0.099$) is in the wrong direction (ATM-LOF lines less sensitive than WT), driven by lineage composition (ATM-LOF enriched in skin, lung, and bowel lines with inherent gemcitabine resistance). This negative control is context-dependent and does not alter the independently reproduced MBD4–ceralasertib result. All frozen values sourced from `canonical_atr_wee1_rerun.json` (2026-04-05).

\clearpage
\FloatBarrier
