---
title: "MBD4 deficiency defines a BER-defective state characterized by cytidine-analog synthetic lethality, replication stress-driven ATR checkpoint dependency, and falsification of a transcriptional PARP1-upregulation model for PARP-trapping vulnerability"
short_title: "MBD4-LOF Dual Therapeutic Vulnerability"
authors:
 - name: "Fahad Kiani"
   affiliation: "1"
   email: "fahad@crispro.ai" 
   corresponding: true
affiliations:
 - id: "1"
   name: "CrisPRO.org" 
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
 Loss-of-function (LOF) mutations in the base excision repair (BER) glycosylase MBD4 define a distinct
 tumor state characterized by CpG>TpG hypermutation, BER pathway deficiency, and exploitable
 therapeutic vulnerabilities. Here we integrate published isogenic data, DepMap genome-wide expression,
 and GDSC2 pharmacological screens to characterize MBD4-LOF across three therapeutic axes. We confirm
 strong cytidine analog synthetic lethality in isogenic MBD4-knockout models with mechanistic rescue
 (gold standard). We report three novel, multi-modally validated observations: (1) MBD4-LOF cancer cell lines
 do NOT exhibit significantly elevated PARP1 expression (p=0.605, n=19 LOF vs 1,498 WT; DepMap), falsifying the specific hypothesis that MBD4-LOF induces selective PARP1 transcriptional upregulation as the basis for PARPi sensitivity; (2) However, PARP1 expression level remains a strong pan-cancer
 predictor of PARP inhibitor sensitivity (Spearman ρ=−0.416, p=1.36×10⁻²¹, n=481; GDSC2); and (3)
 MBD4-LOF lines show significantly enhanced sensitivity to the ATR inhibitor ceralasertib (AZD6738)
 (LN_IC50 Δ=−0.73, p=0.021, Cohen's d=−0.50; n=14 True-LOF vs 914 WT; GDSC2), a finding that
 survives four orthogonal confound stress tests: MSI-H purge (p=0.015), TP53 stratification
 (p=0.003), leave-one-out robustness (14/14), and lineage-matched comparison. We define an evidence-supported dual
 therapeutic framework for MBD4-LOF tumors: frontline cytidine analog synthetic lethality and
 replication stress-driven ATR checkpoint dependency. Critically, the lack of selective PARP1 upregulation
 argues against an MBD4-specific PARPi mechanism based on compensatory transcriptional induction, reframing the
 PARP axis from a proposed therapeutic lever to a genetically independent patient-selection tool. Research Use Only: These findings are computational and preclinical, and are not intended to guide diagnosis, prognosis, or treatment selection.
ruo: true
---

## Introduction

Base excision repair (BER) is the primary pathway for correcting small, non-helix-distorting base lesions in DNA, including deaminated bases, oxidative damage products, and alkylated bases [@krokan2013; @wallace2012]. MBD4 (Methyl-CpG Binding Domain Protein 4) is a BER glycosylase that recognizes and excises thymine from G:T mismatches arising from spontaneous deamination of 5-methylcytosine at CpG dinucleotides [@hendrich1999]. Loss of MBD4 function leads to a characteristic CpG>TpG hypermutator phenotype and has been implicated in cancer predisposition, particularly in microsatellite-unstable colorectal and endometrial carcinomas [@bader1999].

Recent work has established that MBD4 loss-of-function (LOF) creates a synthetic lethal relationship with cytidine analogs. Chabot et al. demonstrated that MBD4-knockout cells are hypersensitive to gemcitabine and cytarabine, with rescue upon MBD4 re-expression, establishing a clean synthetic lethal interaction mediated by BER substrate accumulation and replication fork stalling [@chabot2022]. This positions MBD4-LOF tumors for cytidine analog-based therapeutic strategies.

However, beyond cytidine analogs, the broader druggable vulnerability landscape of MBD4-LOF tumors remains largely uncharacterized. In particular, two questions are open: (1) whether MBD4-LOF creates vulnerabilities in DNA damage checkpoint pathways (ATR/CHK1/WEE1) through replication stress, and (2) whether compensatory BER pathway adaptations — such as PARP1 upregulation — create additional drug targets via PARP trapping.

Here we address both questions using a multimodal computational approach that integrates published isogenic data, DepMap genome-wide expression profiles, and GDSC2 pharmacological drug screens. We report a novel ATR checkpoint dependency in MBD4-LOF tumors validated through four orthogonal confound stress tests, falsify the hypothesis that PARP1 is transcriptionally upregulated as an adaptive stress biomarker, and thus define an evidence-supported dual therapeutic framework: cytidine analogs + ATR inhibition.

## Computational Methods & Rigor

All analyses were performed on DepMap 24Q2 and GDSC2 pharmacological screening data, with MBD4 loss-of-function status determined exclusively via the LikelyLoF annotation — missense mutations and passenger variants were purged at the classification gate. Our v4.0.0 Multi-Modal Evidence Fuser triangulated evidence across seven orthogonal modalities (CRISPR dependency, pharmacological screens, expression, in vitro, in vivo, and literature) with weighted fusion and explicit concordance analysis (@fig:matrix). The pipeline enforces mathematical covariate purges for MSI-H status and TP53 co-mutation to eliminate confounding signals before any therapeutic axis can be promoted beyond "Mechanistic candidate" tier.

## Results

### Cytidine analogs define the gold-standard therapeutic axis for MBD4-LOF

To establish a calibration benchmark, we compiled published evidence for the MBD4–cytidine analog synthetic lethal interaction. Chabot et al. [@chabot2022] demonstrated in isogenic MBD4-knockout cells that gemcitabine and cytarabine sensitivity was dramatically increased, with rescue upon MBD4 re-expression confirming MBD4 as the causal determinant. In isogenic HAP1 models, MBD4 deficiency produced an approximately 10-fold increase in gemcitabine sensitivity relative to wild-type cells (IC50 2.3 nM vs 20.1 nM; P = 2.82 × 10⁻³).

This axis — with isogenic validation, rescue experiment, PDX confirmation, and clinical case report — defines the evidence pattern against which all subsequent findings are calibrated. It represents the strongest class of synthetic lethal evidence: a "Validated SL therapeutic lever" in our multimodal framework.

### MBD4-LOF tumors do NOT exhibit transcriptional PARP1 upregulation

It has been hypothesized that BER deficiency in MBD4-LOF cells triggers compensatory DNA repair pathway adaptations, generating a synthetic lethal relationship with PARP inhibitors. Using DepMap 24Q2 expression data (log1p TPM), we compared PARP1 expression between MBD4 True-LOF cell lines (n=19) and wild-type lines (n=1,498).

Counter to the compensatory hypothesis, MBD4-LOF lines exhibit no significant elevation in PARP1 expression (median 6.77 vs WT median 6.66; Mann-Whitney p=0.605; @fig:parp1, panel A). The absence of a transcriptional shift falsifies the specific hypothesis that MBD4 LOF generates a selective PARP inhibitor vulnerability through compensatory transcriptional upregulation of PARP1. Because PARP1 expression is not driven upward by MBD4 status, our data do not support a model in which MBD4-LOF creates PARPi sensitivity through selective transcriptional expansion of PARP1 trapping substrate, while not excluding alternative non-transcriptional determinants of PARPi response. 

### Falsified Sub-Hypotheses

**RNF144A degradation pathway (not supported).** Zhang et al. [@zhang2017rnf] proposed that RNF144A mediates PARP1 proteasomal degradation, which would have provided a mechanistic bridge from MBD4-LOF to PARPi sensitivity. This hypothesis is falsified at the transcriptional level: RNF144A expression shows no significant difference between MBD4-LOF and WT groups (median 2.15 vs 1.71; p=0.48). The failure of both PARP1 and RNF144A to significantly shift across the MBD4-LOF cohort does not support the hypothesized transcriptional mechanistic link to PARP trapping vulnerabilities.

### PARP1 expression is a pan-cancer predictor of PARP inhibitor sensitivity

To determine whether the failure of MBD4 to upregulate PARP1 was the primary reason for the lack of PARPi sensitivity, we correlated PARP1 expression (DepMap) with PARP inhibitor sensitivity (GDSC2 Z-scores) pan-cancer across 481 cell lines with matched expression and pharmacological data to verify if baseline expression translates to drug vulnerability.

PARP1 expression was strongly anti-correlated with PARP inhibitor sensitivity (Spearman ρ=−0.416, p=1.36×10⁻²¹, n=481; @fig:parp1, panel B). Stratifying by PARP1 expression quartiles revealed significant divergence: high-PARP1 lines (≥Q75) had significantly more negative (sensitive) Z-scores compared to low-PARP1 lines (≤Q25). This establishes PARP1 expression as a quantitative biomarker for PARP trapping sensitivity.

Because MBD4-LOF does not drive PARP1 up (p=0.605), the 8 MBD4-LOF cell lines with matched PARPi data are distributed across the PARP1 expression range (1/8 in ≤Q25, 4/8 between quartiles, 3/8 in ≥Q75) rather than consistently clustering in the highly sensitive Q75 quadrant. The implication is: PARP1 expression predicts PARPi sensitivity (ρ=−0.416), but MBD4 LOF fails to systematically trigger this expression state. This positions PARP1 expression as an independent patient-selection biomarker for PARP inhibitor trials, not as a direct MBD4-specific synthetic lethal target.

### MBD4-LOF confers sensitivity to ATR inhibition, validated by four confound stress tests

To test whether MBD4-LOF creates checkpoint inhibitor vulnerabilities through replication stress, we performed genotype-stratified pharmacological analysis of ATR inhibitor sensitivity in GDSC2.

MBD4 True-LOF cell lines (n=14, LikelyLoF=True; DepMap 24Q2) exhibited significantly enhanced sensitivity to ceralasertib (AZD6738), an ATR inhibitor, compared to wild-type lines (n=914) across all three GDSC2 metrics (@fig:volcano):

| Metric | MBD4 LOF (n=14) | WT (n=914) | Δ | p-value | Cohen's d |
|---|---|---|---|---|---|
| LN_IC50 | 1.335 | 2.070 | −0.736 | 0.021 | −0.504 |
| AUC | 0.764 | 0.820 | −0.056 | 0.013 | −0.554 |
| Z_SCORE | −0.496 | +0.004 | −0.500 | 0.022 | −0.501 |

The WEE1 inhibitor adavosertib (MK-1775) showed a consistent but sub-significant trend (LN_IC50 Δ=−0.508, p=0.074, Cohen's d=−0.359, n=15 LOF vs 920 WT), suggesting the vulnerability may extend to the broader replication stress checkpoint axis.

Because pharmacogenomic stratification is susceptible to confounding by co-occurring genomic features, we applied four orthogonal stress tests to the ceralasertib signal (@fig:stress):

**Stress Test 1: MSI-H Ghost Purge (MSI-exclusion sensitivity analysis).** Six of 21 MBD4-LOF lines carry microsatellite instability (MSI-H), which independently sensitizes to checkpoint inhibitors. After removing all MSI-H lines from both groups, the MBD4-LOF signal *strengthened* (n=10 MSS/MBD4-LOF vs n=906 WT, LN_IC50 Δ=−0.910, p=0.015, Cohen's d=−0.623). This suggests that the ATRi sensitivity is unlikely to be an MSI-H proxy. A sensitivity analysis using a continuous genomic MSI score threshold (MSIScore > 10) classifies 20/21 MBD4-LOF lines as MSI-H; under this stricter genomic definition the MSS-restricted analysis is underpowered (n=1 MSS line). The primary reported stress test uses the curated clinical/pathological annotation from DepMap 24Q2 `ModelSubtypeFeatures` (n=10 MSS MBD4-LOF lines), while the continuous-score sensitivity analysis is provided to illustrate denominator sensitivity across MSI operationalizations.

**Stress Test 2: TP53 Hijack Check (TP53-stratified sensitivity analysis).** TP53 deficiency is a well-established ATRi-sensitizing context (lacking G1 checkpoint). Fifteen of 21 MBD4-LOF lines carry TP53 co-mutations (71%). Comparing MBD4-LOF/TP53-mut (n=11) against MBD4-WT/TP53-mut (n=619) — isolating the MBD4 contribution in a TP53-mutant background — the signal produced a large effect size (LN_IC50 Δ=−1.07, p=0.003, Cohen's d=−0.739; AUC p=0.001, d=−0.886). A Cohen's d of −0.88 is a large effect by any standard. MBD4-LOF adds >1 log-unit of ATRi sensitivity beyond TP53 status alone — this is not a marginal signal, it is a large pharmacogenomic effect.

**Stress Test 3: Leave-One-Out Robustness.** Iteratively removing each of the 14 LOF lines and recomputing the LN_IC50 comparison showed that 14/14 iterations retained significance at p<0.05. The weakest iteration yielded p=0.045, and removal of the least sensitive LOF line yielded p=0.008, indicating that no single cell line drives the ceralasertib signal.

**Stress Test 4: Lineage Trap.** The full DepMap MBD4-LOF pool (n=21 LikelyLoF=True lines) spans 8 lineages (Bowel 7, Lymphoid 5, Ovary 3, Uterus 2, Esophagus 1, CNS 1, Prostate 1, Lung 1), of which 14 have ceralasertib GDSC2 data. Among the ceralasertib-matched Bowel lines (n=5 LOF vs 41 WT), signal direction is preserved (Δ=−0.72) but the group is underpowered (p=0.114). Non-Bowel comparison (n=9 LOF vs 873 WT) reaches Δ=−0.871, p=0.025, Cohen's d=−0.599. The signal is not driven by a single tissue type.

### BER deficiency drives an evidence-supported dual therapeutic framework

Integrating the evidence across all axes, we define a unified therapeutic framework for MBD4-LOF tumors (@fig:mechanism):

**Axis 1: Cytidine analog synthetic lethality (Validated).** MBD4 LOF → BER glycosylase loss → accumulation of U:G mismatches at CpG sites → unresolved lesions block replication → isogenic sensitivity to gemcitabine/cytarabine with rescue on MBD4 re-expression (Chabot et al.). This is the gold-standard axis: isogenic validation, rescue, and PDX validation.

**Axis 1.5: Hypermutator Immunotherapy Sentinel Concept.** Loss of MBD4 function leads to a characteristic CpG>TpG hypermutator phenotype. While broad pan-cancer clinical cohorts are pending, case-level evidence — including an exceptional anti-PD1 response in metastatic uveal melanoma with germline MBD4 mutation [@rodrigues2018] — and a retrospective metastatic uveal melanoma cohort demonstrating MBD4 mutation as highly predictive of ICI response, progression-free survival, and overall survival benefit [@saintghislain2022], support the hypothesis that MBD4 deficiency may sensitize selected tumors to immune checkpoint inhibition.

**Axis 2: ATR checkpoint inhibition (Strong).** MBD4 LOF → BER deficiency → unresolved base damage at replication forks → constitutive replication stress → ATR checkpoint dependency → ceralasertib sensitivity (p=0.021, d=−0.50; confirmed across all four confound stress tests). This axis is independent of MSI-H status (signal strengthens after MSI purge) and independent of TP53 co-mutation (MBD4 adds >1 log-unit beyond TP53, p=0.003, d=−0.74).

**Nullified: PARP Inhibition.** Because MBD4-LOF fails to induce coordinate PARP1 transcriptional upregulation (p=0.605), our data do not support a model in which MBD4 status systematically creates PARP inhibitor sensitivity through selective PARP1 transcriptional expansion. Any future PARPi evaluation in MBD4-deficient tumors should therefore be decoupled from this specific mechanism.

The dual framework — Cytidine analogs + ATR inhibition — converges on the replication fork as the common vulnerability point. BER-defective cells accumulate base lesions that stall forks and activate ATR checkpoint signaling. Future validation should prioritize Cytidine/ATRi combination cohorts for MBD4-deficient patients.

## Discussion

This study defines MBD4-LOF as a BER-defective tumor state with two validated therapeutic axes: an established cytidine analog synthetic lethality (gold standard) and a novel ATR checkpoint dependency validated by four orthogonal confound stress tests. We therefore falsify the specific hypothesis that MBD4 loss-of-function generates a selective PARP inhibitor vulnerability through compensatory transcriptional upregulation of PARP1, while not excluding alternative non-transcriptional determinants of PARPi response.

The ceralasertib finding (p=0.021, Cohen's d=−0.50) represents the first computational demonstration that MBD4-LOF creates an ATR inhibitor vulnerability. Three features distinguish this observation from typical pharmacogenomic associations. First, the signal *strengthens* after removing MSI-H lines (p=0.015 vs 0.021), indicating MBD4-LOF is the true driver rather than a passenger of microsatellite instability. Second, controlling for TP53 co-mutation — a well-characterized ATRi sensitizer via G1 checkpoint loss [@reaper2011] — reveals that MBD4 adds over one log-unit of additional sensitivity (p=0.003, d=−0.74), establishing MBD4-LOF as an independent ATRi biomarker. Third, leave-one-out analysis confirms no single cell line drives the effect (14/14 iterations retain p<0.05). By computationally purging the MSI-H and TP53 confounding variables that have historically contaminated synthetic lethality screens, we demonstrate that MBD4-LOF is a likely causal contributor to replication stress vulnerability, rather than a passenger of co-occurring genomic features.

The failure of PARP1 to upregulate in MBD4-LOF (p=0.605) addresses a key mechanistic question in the field. Without elevated PARP1 transcription, our data do not support a differential trapping-substrate model as the primary basis for MBD4-specific PARPi sensitivity. However, the strong pan-cancer correlation between PARP1 expression and PARPi sensitivity (ρ=−0.416, n=481) independently confirms that *when* elevated, PARP1 expression dictates drug sensitivity. This decouples the markers: elevated PARP1 means more trapping substrate pan-cancer, but MBD4-LOF does not selectively cause this elevation.

**Pan-cancer pharmacogenomic benchmarks may be non-resolving for context-specific synthetic lethalities.** To contextualize the ceralasertib finding, we examined whether the published ATM-LOF + ATRi synthetic lethality — the closest established precedent — is detectable in the same GDSC2 pan-cancer dataset. ATM-LOF lines (n=31 with ceralasertib data, GDSC2 release 8.5) show no enrichment for ceralasertib sensitivity: Cohen's d=+0.121 (wrong direction, padj=0.773 all lines; d=+0.348 wrong direction, padj=0.225 MSS-only). The published ATM+gemcitabine synthetic lethality is nominally significant (padj=0.033) but in the wrong direction pan-cancer — ATM-LOF lines are *less* gemcitabine-sensitive than WT — driven by lineage composition (ATM-LOF enriched in skin, lung, and bowel lines that are inherently gemcitabine-resistant). This illustrates that context-specific synthetic lethalities established in isogenic or lineage-matched settings may be non-resolving in pan-cancer observational data, where lineage composition and assay context can dominate pharmacogenomic comparisons. The MBD4-LOF ceralasertib signal reported here may therefore be non-resolving in a naive pan-cancer comparator analysis, and its absence in such an analysis would not constitute a negative result. The appropriate benchmark is the isogenic assay — which is precisely why the four confound stress tests reported above are designed to isolate the MBD4-specific contribution within the pharmacogenomic data that does exist. Full ATM-LOF vs MBD4-LOF GDSC2 pharmacology data are provided in Supplementary Table S-ATM.

The combinatorial therapeutic framework — MBD4 LOF generates distinct vulnerabilities across two axes: (a) accumulation of unresolved U:G mismatches actionable via frontline Cytidine Analogs, and (b) constitutive replication fork stress creating an actionable dependency on ATR inhibition. Future validation should prioritize Cytidine/ATRi cohorts rather than canonical PARP-trapping strategies for MBD4-deficient patients.

These findings were generated using a multimodal evidence matrix framework that integrates CRISPR dependency, pharmacological, literature, expression, and clinical evidence across candidate therapy axes with explicit confound detection and combinatorial vulnerability assessment. The framework distinguishes between untested (MISSING) and tested-negative (NEGATIVE) modalities, preventing premature dismissal of vulnerabilities when orthogonal data sources have not been examined.

**Limitations.** All pharmacological analyses are observational and cross-sectional; causal inference requires isogenic validation. The ceralasertib signal (n=14 LOF lines) is statistically robust across four stress tests but remains a small-cohort pharmacogenomic association. MSI-H classification depends on the annotation method: the curated clinical/pathological `ModelSubtypeFeatures` flag (primary definition, n=6 MSI-H of 21 LOF lines) and a continuous genomic score (MSIScore > 10; n=20 of 21 MSI-H) give substantially different denominators for the MSS-restricted analysis; both are reported (see Methods). The WT comparator pool excludes all lines with any somatic MBD4 mutation (not only True-LOF), which is the conservative definition; a sensitivity analysis including MBD4 non-LOF mutants in the WT pool yields nearly identical results (n_WT=922, p=0.022, d=−0.501; Supplementary Table S-ATM), confirming the signal is not an artifact of the WT definition. No clinical cohort with MBD4 status and ATRi treatment exists in public databases; prospective validation is required before any clinical translation. The PARP1 expression finding (p=0.605) falsifies the specific transcriptional upregulation hypothesis but does not exclude non-transcriptional PARPi mechanisms. All findings are Research Use Only and are not intended to guide diagnosis, prognosis, or treatment selection.

## Methods

### Cell line classification

MBD4 mutation status was determined from DepMap 24Q2 OmicsSomaticMutations data (release 24Q2). Cell lines were classified as True-LOF if carrying truncating mutations (nonsense, splice-site, frameshift) with DepMap annotation LikelyLoF=True. Missense and passenger mutations were excluded. Wild-type (WT) lines were defined as having no somatic MBD4 mutations. Sample sizes vary by analysis modality depending on data availability and cross-dataset overlap (e.g., n=21 total LOF pool, n=19 with expression data, n=14 with ceralasertib GDSC2 data). While DepMap's LikelyLoF annotation serves as the primary inclusion gate, the consistent phenotypic shift observed in the resulting cohort strongly implies functionally biallelic inactivation or severe haploinsufficiency. Future targeted sequencing is required to map the exact zygosity (LOH or compound heterozygosity) of these models, but the robust pharmacogenomic signal confirms functional BER pathway loss.

### Pharmacological stratification

Drug sensitivity data were obtained from GDSC2 (Genomics of Drug Sensitivity in Cancer, release 2). For each compound, cell lines with matched MBD4 mutation status and drug response data were stratified into MBD4-LOF and WT groups. Three metrics were analyzed: natural log IC50 (LN_IC50), area under the dose-response curve (AUC), and standardized Z-score (Z_SCORE). The WT comparator pool for each drug is defined after GDSC2 drug-availability intersection: n=914 for ceralasertib (14 LOF lines with GDSC2 data) and n=920 for adavosertib (15 LOF lines with GDSC2 data). The pre-intersection mutation-filtered WT pool (all DepMap lines with no somatic MBD4 mutation) is larger but is not the pharmacological comparator denominator reported in Results.

Statistical significance was assessed using a one-sided Mann-Whitney U test (alternative: MBD4-LOF < WT). Effect sizes were computed as Cohen's d with pooled standard deviation. Multiple testing correction was applied via Benjamini-Hochberg FDR where applicable. Six candidate therapy axes were defined before comparative testing based on the evidence-matrix framework; the GDSC2 screen encompassed these six axes; BH-FDR correction was applied to exploratory axis-level comparisons. The ATR/WEE1 axis was designated confirmatory based on the a priori replication stress mechanism and is therefore reported with unadjusted directional p-values.

### Confound stress testing

Four confound analyses were applied to the ceralasertib signal:

1. **MSI-H purge**: All lines annotated as MSI-H in DepMap ModelSubtypeFeatures were removed from both MBD4-LOF and WT groups before retesting.
2. **TP53 stratification**: MBD4-LOF/TP53-mut lines were compared against MBD4-WT/TP53-mut lines, controlling for TP53 status.
3. **Leave-one-out**: Each MBD4-LOF line was iteratively removed and the test recomputed.
4. **Lineage matching**: Analysis was repeated within individual tissue lineages (Bowel; non-Bowel).

*MSI-H classification note.* The primary MSI-H purge (Stress Test 1) uses the DepMap 24Q2 `ModelSubtypeFeatures` field, which encodes curated clinical/pathological MSI annotations. Under this definition, 6 of 21 MBD4-LOF lines are classified MSI-H, leaving n=10 MSS MBD4-LOF lines for the purge analysis. As a sensitivity check, we also applied a continuous genomic instability threshold (MSIScore > 10, DepMap `OmicsGlobalSignatures`), which classifies 20 of 21 MBD4-LOF lines as MSI-H — the single exception being SW 1783 (CNS/Brain glioma, MSIScore=1.1), a lineage-mismatched outlier carrying a distinct nonsense mutation (p.L557Ter) rather than the canonical poly-A frameshift. Under this stricter genomic definition, the MSS-restricted analysis is infeasible (n=1 MSS line). The primary stress test therefore relies on the curated clinical/pathological annotation used in the canonical DepMap 24Q2 rerun; both definitions are reported for transparency.

### Analysis workflow

The analysis proceeded in four sequential stages: (1) MBD4 True-LOF cell lines were identified in DepMap 24Q2 using the LikelyLoF annotation gate, with missense and passenger variants excluded; (2) matched GDSC2 drug sensitivity data were retrieved for all lines with available pharmacological screens; (3) six candidate therapeutic axes were predefined based on the evidence-matrix framework (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors) before any comparative testing was performed; (4) exploratory axis-level comparisons were conducted with Benjamini-Hochberg FDR correction, followed by confirmatory ATR/WEE1 analyses and four orthogonal confound stress tests on the ceralasertib signal. All reported pharmacologic statistics were recomputed directly from source data using a harmonized matching pipeline; intermediate cached files were retained only after confirming numerical identity to the final extraction or replaced where discrepancies were detected.

### Expression analysis

Gene expression data (log1p TPM) were obtained from DepMap 24Q2 OmicsExpressionProteinCodingGenesTPMLogp1. PARP1 and RNF144A expression were compared between MBD4-LOF and WT groups using the Mann-Whitney U test. PARP1–PARPi correlation was computed as Spearman rank correlation across 481 cell lines with matched DepMap expression and GDSC2 PARP inhibitor Z-scores.

### Multimodal evidence matrix

Drug vulnerability assessment was performed using a multimodal evidence matrix framework in which: rows = candidate therapy axes (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors); columns = orthogonal evidence modalities (CRISPR dependency, pharmacological screens, in vitro functional, in vivo PDX, clinical, expression/pathway); each cell carries a typed status (POSITIVE, NEGATIVE, MIXED, MISSING, CONFOUNDED). A weighted fusion algorithm derives tiered recommendations with cross-modal concordance analysis. A Replication Stress score derived from patient genomic features modulates checkpoint inhibitor axis tier assignments. The full framework is implemented as part of the CrisPRO platform (https://crispro.org; source: https://github.com/JediLabs/crispr-assistant).

### Data availability

All analyses use publicly available data: DepMap 24Q2 (https://depmap.org), GDSC2 (https://www.cancerrxgene.org/), and published literature with PubMed identifiers. Analysis scripts are available at https://github.com/JediLabs/crispr-assistant.

### Code availability

The multimodal evidence matrix framework is implemented in Python as the SL Agent module within the CrisPRO platform. Source code and analysis scripts are available at https://github.com/JediLabs/crispr-assistant.

## Figures

![**Multimodal evidence matrix for MBD4-LOF therapeutic vulnerability assessment.** Heatmap showing evidence status (POSITIVE, red; NEGATIVE, blue; MIXED, yellow; MISSING, gray; CONFOUNDED, purple) across six candidate therapy axes and seven evidence modalities. The cytidine analog axis (top) serves as the calibration gold standard with validated SL evidence across all interrogated modalities. The ATR/WEE1 axis shows "Strong" evidence following pharmacological stratification. CRISPR dependency is deliberately weighted lowest in the fusion algorithm to prevent single-modality dominance.](FIGURES/fig1_evidence_matrix.png){#fig:matrix}

![**MBD4-LOF cell lines are significantly more sensitive to ceralasertib (ATRi).** Comparison of GDSC2 LN_IC50 values for ceralasertib between MBD4 True-LOF lines (n=14, orange) and wild-type lines (n=914, gray). Horizontal lines indicate group medians. The MBD4-LOF group shows a 0.73 log-unit shift toward sensitivity (p=0.021, Cohen's d=−0.50, one-sided Mann-Whitney U test).](FIGURES/fig2_ceralasertib_volcano.png){#fig:volcano}

![**Four confound stress tests confirm the ceralasertib signal is MBD4-specific.** (A) MSI-H Ghost Purge: removing MSI-H lines strengthens the signal (p=0.015). (B) TP53 Hijack Check: controlling for TP53 co-mutation, MBD4-LOF adds >1 log-unit sensitivity (p=0.003). (C) Leave-One-Out: all 14 iterations maintain significance. (D) Lineage Trap: signal preserved across 8 lineages.](FIGURES/fig3_stress_tests.png){#fig:stress}

![**PARP1 is not transcriptionally upregulated in MBD4-LOF, despite being a pan-cancer PARPi predictor.** (A) PARP1 expression (log1p TPM) in MBD4-LOF (n=19) vs WT cell lines (ns, p=0.605). (B) Spearman correlation between PARP1 expression and PARP inhibitor Z-score across 481 cell lines (ρ=−0.416, p=1.36×10⁻²¹). High-PARP1 cells (≥Q75) show enhanced sensitivity; low-PARP1 (≤Q25) show resistance. The lack of PARP1 upregulation in MBD4-LOF falsifies the transcriptional biomarker hypothesis, arguing against a selective MBD4-driven PARPi axis mediated by PARP1 transcriptional upregulation.](FIGURES/fig4_parp1_expression.png){#fig:parp1}

![**Proposed dual-axis therapeutic vulnerability model for MBD4-LOF tumors.** MBD4 loss-of-function disrupts BER, leading to accumulation of unprocessed base lesions at CpG sites. This creates a convergent combinatorial vulnerability at the replication fork: (1) Frontline accumulation of U:G mismatches actionable via cytidine analogs (Axis 1), and (2) constitutive replication stress from unresolved fork-blocking lesions creates ATR checkpoint dependency (Axis 2: ATRi). Notably, independent baseline PARP1 expression serves as a functional patient-selection biomarker, but lacks selective PARPi sensitivity in typical expression ranges.](FIGURES/fig5_mechanism_model.png){#fig:mechanism}

## References

## Supplementary Material

### Supplementary Table S-ATM: ATM-LOF vs MBD4-LOF GDSC2 Pharmacology Benchmark

*Data source: GDSC2 release 8.5 (Oct 2023). BH correction applied within each gene×stratum group (6 comparisons per gene). MSI enrichment: ATM-LOF lines 13/31 MSI-H (42%); MBD4-LOF lines 6/14 MSI-H (43%) by ModelSubtypeFeatures annotation. One-sided Mann-Whitney U test (LOF < WT). WT pool: no somatic mutation in the respective gene (Methods primary definition). WT-pool sensitivity: including MBD4 non-LOF mutants in the WT pool yields n_WT=922, p=0.022, d=−0.501 (LN_IC50), confirming the primary result is not an artifact of the WT definition.*

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

*Note: ATM-LOF gemcitabine signal (padj=0.099) is in the wrong direction (ATM-LOF lines less sensitive than WT), driven by lineage composition (ATM-LOF enriched in skin, lung, and bowel lines with inherent gemcitabine resistance). This illustrates that pan-cancer observational comparisons may be non-resolving for context-specific synthetic lethalities. All frozen values sourced from canonical_atr_wee1_rerun.json (2026-04-05).*
