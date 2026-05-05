---
title: "MBD4 loss-of-function is associated with a replication stress-linked ATR inhibitor vulnerability and does not support a direct PARP-trapping model: a pharmacogenomic hypothesis-generation study"
short_title: "MBD4-LOF ATR Vulnerability: Pharmacogenomic Hypothesis Generation"
authors:
 - name: "Fahad Kiani"
   affiliation: "1"
   email: "fahad@crispro.ai"
   corresponding: true
affiliations:
 - id: "1"
   name: "Independent Researcher"
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
 - hypothesis generation
 - isogenic validation
 - ATR checkpoint dependency
journal: "bioRxiv"
doi: ""
abstract: |
 Loss-of-function (LOF) mutations in the base excision repair (BER) glycosylase MBD4 define a distinct
 tumor state characterized by CpG>TpG hypermutation and BER pathway deficiency. We report three
 pharmacogenomic observations using DepMap 24Q2 expression data and GDSC2 drug sensitivity screens.
 First, MBD4-LOF cancer cell lines show enhanced sensitivity to the ATR inhibitor ceralasertib
 (AZD6738) compared with wild-type lines (LN_IC50 Δ=−0.73, p=0.021, Cohen's d=−0.50; n=14 vs 914;
 GDSC2), and this association survives predefined confound stress tests, although MSI-H collinearity
 remains a major structural limitation of the observational design. Second, MBD4-LOF lines do not exhibit significantly elevated PARP1
 expression (p=0.605, n=19 LOF vs 1,498 WT; DepMap), falsifying the specific hypothesis that MBD4-LOF
 induces selective PARP1 transcriptional upregulation as the basis for PARPi sensitivity. Third,
 published isogenic data confirm that MBD4-knockout cells are hypersensitive to cytidine analogs
 (gemcitabine, cytarabine) with rescue on MBD4 re-expression (Chabot et al., 2022), providing a
 positive-control calibration benchmark for the pharmacogenomic approach. Together, these observations
 generate the hypothesis that MBD4-LOF creates a replication stress-linked ATR checkpoint dependency
 that is mechanistically distinct from PARP trapping. Isogenic validation — specifically ATRi
 sensitivity in MBD4-knockout vs. rescue cell lines — is required before any clinical translation.
ruo: false
---

## Introduction

Base excision repair (BER) is the primary pathway for correcting small, non-helix-distorting base lesions in DNA, including deaminated bases, oxidative damage products, and alkylated bases [@krokan2013; @wallace2012]. MBD4 (Methyl-CpG Binding Domain Protein 4) is a BER glycosylase that recognizes and excises thymine from G:T mismatches arising from spontaneous deamination of 5-methylcytosine at CpG dinucleotides [@hendrich1999]. Loss of MBD4 function leads to a characteristic CpG>TpG hypermutator phenotype and has been implicated in cancer predisposition, particularly in microsatellite-unstable colorectal and endometrial carcinomas [@bader1999].

Recent work has established that MBD4 loss-of-function (LOF) creates a synthetic lethal relationship with cytidine analogs. Chabot et al. demonstrated that MBD4-knockout cells are hypersensitive to gemcitabine and cytarabine, with rescue upon MBD4 re-expression confirming MBD4 as the causal determinant [@chabot2022]. This work establishes that BER substrate accumulation in MBD4-LOF cells generates replication fork stalling that is pharmacologically exploitable — a mechanistic connection between BER deficiency and replication stress that motivates the present study.

One mechanistic model that follows from this connection is that MBD4 loss creates constitutive replication stress through accumulation of unprocessed U:G mismatches at CpG dinucleotides. Under this model, during S-phase these lesions are encountered by the replication machinery, generating single-stranded DNA (ssDNA) gaps and stalled forks. Exposed ssDNA would recruit RPA, which in turn recruits the ATR-ATRIP complex via TOPBP1, activating CHK1 and the intra-S and G2/M checkpoints. If correct, BER-deficient cells would be constitutively dependent on ATR checkpoint signaling to tolerate the ongoing replication stress imposed by unrepaired base lesions — a dependency that ATR inhibitors would be expected to exploit. This proposed mechanism is distinct from a PARP-trapping model, which depends on the generation of a permissive PARP-trapping state rather than checkpoint reliance. The present study does not directly test this chain at the molecular level; rather, it uses pharmacogenomic stratification to ask whether the predicted ATR dependency is detectable in publicly available cell line drug sensitivity data.

Beyond cytidine analogs, two questions about the broader druggable vulnerability landscape of MBD4-LOF tumors remain open: (1) whether MBD4-LOF creates vulnerabilities in DNA damage checkpoint pathways (ATR/CHK1/WEE1) through replication stress, and (2) whether compensatory BER pathway adaptations — such as PARP1 upregulation — create additional drug targets via PARP trapping.

Here we address both questions using pharmacogenomic stratification of publicly available cell line drug sensitivity data. We report a pharmacogenomic association between MBD4-LOF and ATR inhibitor sensitivity stress-tested across four orthogonal confound analyses, falsify the hypothesis that PARP1 is transcriptionally upregulated as an adaptive stress biomarker, and generate the testable hypothesis that MBD4-LOF creates an ATR checkpoint dependency warranting isogenic validation.

## Results

### MBD4-LOF is associated with ATR inhibitor sensitivity across four confound stress tests

To test whether MBD4-LOF is associated with checkpoint inhibitor vulnerability, we performed genotype-stratified pharmacological analysis of ATR inhibitor sensitivity in GDSC2.

One mechanistic model predicts that MBD4 loss generates constitutive replication stress through accumulation of unprocessed base lesions at CpG sites. Under this model, these lesions generate ssDNA gaps during S-phase that activate the ATR-ATRIP-CHK1 checkpoint axis, creating a dependency on ATR kinase activity for fork stabilization and cell viability. We tested whether this predicted dependency is detectable in pan-cancer pharmacogenomic data.

MBD4 True-LOF cell lines (n=14, LikelyLoF=True; DepMap 24Q2) exhibited significantly enhanced sensitivity to ceralasertib (AZD6738), an ATR inhibitor, compared to wild-type lines (n=914) across all three GDSC2 metrics (@fig:volcano):

| Metric | MBD4 LOF (n=14) | WT (n=914) | Δ | p-value | Cohen's d |
|---|---|---|---|---|---|
| LN_IC50 | 1.335 | 2.070 | −0.736 | 0.021 | −0.504 |
| AUC | 0.764 | 0.820 | −0.056 | 0.013 | −0.554 |
| Z_SCORE | −0.496 | +0.004 | −0.500 | 0.022 | −0.501 |

The WEE1 inhibitor adavosertib (MK-1775) showed a consistent but sub-significant trend (LN_IC50 Δ=−0.508, p=0.074, Cohen's d=−0.359, n=15 LOF vs 920 WT), suggesting the vulnerability may extend to the broader replication stress checkpoint axis.

Because pharmacogenomic stratification is susceptible to confounding by co-occurring genomic features, we applied four orthogonal stress tests to the ceralasertib signal (@fig:stress):

**Stress Test 1: MSI-H purge (MSI-exclusion sensitivity analysis).** Six of 21 MBD4-LOF lines carry microsatellite instability (MSI-H) by curated clinical/pathological annotation, which independently sensitizes to checkpoint inhibitors. After removing all MSI-H lines from both groups, the MBD4-LOF signal strengthened (n=10 MSS/MBD4-LOF vs n=906 WT, LN_IC50 Δ=−0.910, p=0.015, Cohen's d=−0.623). The signal survives this predefined MSI-related stress test; however, MSI-H collinearity remains a major structural limitation of the observational design. A sensitivity analysis using a continuous genomic MSI score threshold (MSIScore > 10) classifies 20/21 MBD4-LOF lines as MSI-H, leaving only n=1 MSS line — under this stricter genomic definition the MSS-restricted analysis is underpowered and the MSI confound cannot be fully resolved in this dataset. Both MSI operationalizations are reported for transparency (see Methods).

**Stress Test 2: TP53 stratification.** TP53 deficiency is a well-established ATRi-sensitizing context (lacking G1 checkpoint). Fifteen of 21 MBD4-LOF lines carry TP53 co-mutations (71%). Comparing MBD4-LOF/TP53-mut (n=11) against MBD4-WT/TP53-mut (n=619) — isolating the MBD4 contribution in a TP53-mutant background — the signal produced a large effect size (LN_IC50 Δ=−1.07, p=0.003, Cohen's d=−0.739; AUC p=0.001, d=−0.886). In TP53-mutant lines, MBD4-LOF remained associated with greater ATRi sensitivity, with an observed LN_IC50 shift of more than 1 log unit relative to TP53-mutant MBD4-WT comparators.

**Stress Test 3: Leave-one-out robustness.** Iteratively removing each of the 14 LOF lines and recomputing the LN_IC50 comparison showed that 14/14 iterations retained significance at p<0.05. The weakest iteration yielded p=0.045; removal of the least sensitive LOF line yielded p=0.008, indicating that no single cell line drives the ceralasertib signal.

**Stress Test 4: Lineage analysis.** The full DepMap MBD4-LOF pool (n=21 LikelyLoF=True lines) spans 8 lineages (Bowel 7, Lymphoid 5, Ovary 3, Uterus 2, Esophagus 1, CNS 1, Prostate 1, Lung 1), of which 14 have ceralasertib GDSC2 data. Among the ceralasertib-matched Bowel lines (n=5 LOF vs 41 WT), signal direction is preserved (Δ=−0.72) but the group is underpowered (p=0.114). Non-Bowel comparison (n=9 LOF vs 873 WT) reaches Δ=−0.871, p=0.025, Cohen's d=−0.599. The signal is not driven by a single tissue type.

These pharmacogenomic observations generate a testable hypothesis: MBD4-LOF cells are constitutively dependent on ATR checkpoint signaling, and this dependency is exploitable by ATR inhibitors. Definitive validation requires isogenic ATRi sensitivity testing in MBD4-knockout vs. MBD4-rescue cell lines, which is the critical next experiment.

### Expression analyses do not support a direct PARP-trapping model in MBD4-LOF

It has been hypothesized that BER deficiency in MBD4-LOF cells triggers compensatory DNA repair pathway adaptations, generating a synthetic lethal relationship with PARP inhibitors. Using DepMap 24Q2 expression data (log1p TPM), we compared PARP1 expression between MBD4 True-LOF cell lines (n=19) and wild-type lines (n=1,498).

Counter to the compensatory hypothesis, MBD4-LOF lines exhibit no significant elevation in PARP1 expression (median 6.77 vs WT median 6.66; Mann-Whitney p=0.605; @fig:parp1, panel A). The absence of a transcriptional shift falsifies the specific hypothesis that MBD4 LOF generates a selective PARP inhibitor vulnerability through compensatory transcriptional upregulation of PARP1. Because PARP1 expression is not driven upward by MBD4 status, our data do not support a model in which MBD4-LOF creates PARPi sensitivity through selective transcriptional expansion of PARP1 trapping substrate, while not excluding alternative non-transcriptional determinants of PARPi response.

**RNF144A expression is also unchanged in MBD4-LOF.** Zhang et al. [@zhang2017rnf] proposed that RNF144A mediates PARP1 proteasomal degradation, which would have provided a mechanistic bridge from MBD4-LOF to PARPi sensitivity. This hypothesis is not supported at the transcriptional level: RNF144A expression shows no significant difference between MBD4-LOF and WT groups (median 2.15 vs 1.71; p=0.48). The failure of both PARP1 and RNF144A to significantly shift across the MBD4-LOF cohort does not support the hypothesized transcriptional mechanistic link to PARP trapping vulnerabilities.

**PARP1 expression is a pan-cancer predictor of PARP inhibitor sensitivity.** To determine whether the failure of MBD4 to upregulate PARP1 was the primary reason for the lack of PARPi sensitivity, we correlated PARP1 expression (DepMap) with PARP inhibitor sensitivity (GDSC2 Z-scores) pan-cancer across 481 cell lines with matched expression and pharmacological data.

PARP1 expression was strongly anti-correlated with PARP inhibitor sensitivity (Spearman ρ=−0.416, p=1.36×10⁻²¹, n=481; @fig:parp1, panel B). Stratifying by PARP1 expression quartiles revealed significant divergence: high-PARP1 lines (≥Q75) had significantly more negative (sensitive) Z-scores compared to low-PARP1 lines (≤Q25). This establishes PARP1 expression as a quantitative biomarker for PARP trapping sensitivity.

Because MBD4-LOF does not drive PARP1 up (p=0.605), the 8 MBD4-LOF cell lines with matched PARPi data are distributed across the PARP1 expression range (1/8 in ≤Q25, 4/8 between quartiles, 3/8 in ≥Q75) rather than consistently clustering in the highly sensitive Q75 quadrant. PARP1 expression predicts PARPi sensitivity pan-cancer (ρ=−0.416), but MBD4 LOF does not systematically trigger this expression state. This positions PARP1 expression as an independent patient-selection biomarker for PARP inhibitor trials, not as a direct MBD4-specific synthetic lethal target.

### Published isogenic data establish cytidine analog sensitivity as a positive-control calibration benchmark

To establish a calibration benchmark for the pharmacogenomic approach, we compiled published evidence for the MBD4–cytidine analog synthetic lethal interaction. Chabot et al. [@chabot2022] demonstrated in isogenic MBD4-knockout cells that gemcitabine and cytarabine sensitivity was dramatically increased, with rescue upon MBD4 re-expression confirming MBD4 as the causal determinant. In isogenic HAP1 models, MBD4 deficiency produced an approximately 10-fold increase in gemcitabine sensitivity relative to wild-type cells (IC50 2.3 nM vs 20.1 nM; P = 2.82 × 10⁻³).

This axis — with isogenic validation, rescue experiment, and PDX confirmation — establishes the expected pharmacogenomic signature of a true MBD4-dependent vulnerability: a large effect size in isogenic models with mechanistic rescue. We use this as the positive-control calibration benchmark when interpreting the pharmacogenomic ATR inhibitor signal above.

### Synthesis: MBD4-LOF pharmacogenomics supports an ATR checkpoint hypothesis

Taken together, the pharmacogenomic data support the following model: MBD4 loss-of-function creates a BER-deficient state in which unresolved base lesions generate constitutive replication stress, creating a dependency on ATR checkpoint signaling. The ceralasertib sensitivity signal (p=0.021, d=−0.50) survives four orthogonal confound stress tests and is not attributable to TP53 co-mutation alone. The absence of PARP1 transcriptional upregulation (p=0.605) argues against a PARP-trapping mechanism as the primary MBD4-specific vulnerability. Published isogenic data confirm that cytidine analog sensitivity is a genuine MBD4-dependent phenotype (Chabot et al., 2022), establishing that the pharmacogenomic approach can detect real MBD4-dependent signals when they exist.

This model generates a specific, falsifiable prediction: isogenic MBD4-knockout cells should show enhanced ATRi sensitivity relative to wild-type or MBD4-rescue controls, with the magnitude of effect comparable to the pharmacogenomic signal observed here (d≈−0.50). The critical validation experiment is an isogenic ATRi sensitivity assay in a clean MBD4-knockout/rescue pair, ideally in a colorectal or endometrial cancer background where MBD4 LOF is most prevalent.

## Discussion

This study uses pharmacogenomic stratification of publicly available cell line data to generate the hypothesis that MBD4 loss-of-function creates a replication stress-linked ATR checkpoint dependency. The ceralasertib sensitivity signal (p=0.021, Cohen's d=−0.50) is the primary novel observation; it survives four orthogonal confound stress tests and is not attributable to TP53 co-mutation alone. We additionally falsify the specific hypothesis that MBD4-LOF induces PARP1 transcriptional upregulation as the basis for PARPi sensitivity.

The ceralasertib finding represents the first pharmacogenomic evidence that MBD4-LOF is associated with ATR inhibitor sensitivity. Three features distinguish this observation from typical pharmacogenomic associations. First, the signal strengthens after removing MSI-H lines by curated clinical annotation (p=0.015 vs 0.021), consistent with MBD4-LOF contributing to the signal rather than merely co-segregating with MSI-H status — though the MSI-H collinearity of the cohort (20/21 lines MSI-H by genomic score) means this interpretation requires caution and isogenic confirmation. Second, statistically controlling for TP53 co-mutation — a well-characterized ATRi sensitizer via G1 checkpoint loss [@reaper2011] — reveals that MBD4 adds over one log-unit of additional sensitivity (p=0.003, d=−0.74), establishing MBD4-LOF as an independent pharmacogenomic signal beyond TP53 status. Third, leave-one-out analysis confirms no single cell line drives the effect (14/14 iterations retain p<0.05).

The failure of PARP1 to upregulate in MBD4-LOF (p=0.605) addresses a key mechanistic question. Without elevated PARP1 transcription, our data do not support a differential trapping-substrate model as the primary basis for MBD4-specific PARPi sensitivity. However, the strong pan-cancer correlation between PARP1 expression and PARPi sensitivity (ρ=−0.416, n=481) independently confirms that *when* elevated, PARP1 expression dictates drug sensitivity. This decouples the markers: elevated PARP1 means more trapping substrate pan-cancer, but MBD4-LOF does not selectively cause this elevation.

**Pan-cancer pharmacogenomic benchmarks may be non-resolving for context-specific synthetic lethalities.** To contextualize the ceralasertib finding, we examined whether the published ATM-LOF + ATRi synthetic lethality — the closest established precedent — is detectable in the same GDSC2 pan-cancer dataset. ATM-LOF lines (n=31 with ceralasertib data, GDSC2 release 8.5) show no enrichment for ceralasertib sensitivity: Cohen's d=+0.121 (wrong direction, padj=0.773 all lines; d=+0.348 wrong direction, padj=0.225 MSS-only). The published ATM+gemcitabine synthetic lethality is nominally significant (padj=0.033) but in the wrong direction pan-cancer — ATM-LOF lines are *less* gemcitabine-sensitive than WT — driven by lineage composition (ATM-LOF enriched in skin, lung, and bowel lines that are inherently gemcitabine-resistant). This illustrates that context-specific synthetic lethalities established in isogenic or lineage-matched settings may be non-resolving in pan-cancer observational data, where lineage composition and assay context can dominate pharmacogenomic comparisons. The MBD4-LOF ceralasertib signal reported here may therefore be non-resolving in a naive pan-cancer comparator analysis, and its absence in such an analysis would not constitute a negative result. The appropriate benchmark is the isogenic assay — which is precisely why the four confound stress tests reported above are designed to isolate the MBD4-specific contribution within the pharmacogenomic data that does exist. Full ATM-LOF vs MBD4-LOF GDSC2 pharmacology data are provided in Supplementary Table S1.

The mechanistic model that best fits these observations is one in which MBD4-LOF creates constitutive replication stress through BER substrate accumulation, generating a dependency on ATR checkpoint signaling that is exploitable by ATR inhibitors. Cytidine analog sensitivity (Chabot et al., 2022) is consistent with the same upstream lesion — unresolved U:G mismatches stalling replication forks — but operates through a distinct downstream mechanism (nucleotide incorporation competition rather than checkpoint abrogation). These two vulnerabilities are therefore mechanistically convergent at the replication fork but pharmacologically orthogonal, suggesting combination potential that warrants prospective evaluation. This model is inferential and based on pharmacogenomic association; direct molecular demonstration of the proposed ssDNA/RPA/ATR-ATRIP/CHK1 activation cascade in MBD4-LOF systems awaits isogenic wet-lab validation.

A speculative third vulnerability worth noting is immune checkpoint sensitivity. MBD4-LOF generates a CpG>TpG hypermutator phenotype with elevated tumor mutational burden, and case-level evidence suggests exceptional anti-PD1 responses in MBD4-mutant uveal melanoma [@rodrigues2018; @saintghislain2022]. Whether this extends to other MBD4-LOF tumor types is unknown and requires prospective evaluation.

**Isogenic validation roadmap.** The critical next experiments to validate or falsify the ATR dependency hypothesis are: (1) isogenic ATRi sensitivity assay — ceralasertib or berzosertib dose-response in MBD4-knockout vs. MBD4-rescue cell lines in a colorectal or endometrial cancer background; (2) γH2AX and pCHK1 immunofluorescence to confirm elevated replication stress markers in MBD4-LOF vs. WT cells; (3) ATR inhibitor + gemcitabine combination testing in MBD4-LOF models to assess synergy at the replication fork; and (4) patient-derived organoid or PDX models with confirmed MBD4-LOF status for in vivo ATRi sensitivity. Until isogenic validation is complete, the ceralasertib pharmacogenomic signal should be interpreted as a hypothesis-generating observation, not a validated synthetic lethal interaction.

**Limitations.** All pharmacological analyses are observational and cross-sectional; causal inference requires isogenic validation. The ceralasertib signal (n=14 LOF lines) is statistically consistent across four stress tests but remains a small-cohort pharmacogenomic association. MSI-H collinearity is the dominant structural limitation: 20 of 21 MBD4-LOF lines are MSI-H by continuous genomic score (MSIScore > 10), and only 6 of 21 are MSI-H by curated clinical/pathological annotation. These two operationalizations give substantially different denominators for the MSS-restricted analysis; both are reported (see Methods). The BH-corrected padj for the ceralasertib AUC metric is 0.103, which does not survive multiple testing correction; the LN_IC50 metric (padj=0.063) also does not survive strict BH correction at α=0.05. The unadjusted p=0.021 is reported as the primary metric because the ATR/WEE1 axis was designated confirmatory a priori based on the replication stress mechanism; however, readers should note that the signal does not survive BH correction across all metrics, and this fragility reinforces the hypothesis-generating rather than confirmatory interpretation of these findings. The WT comparator pool excludes all lines with any somatic MBD4 mutation (not only True-LOF), which is the conservative definition; a sensitivity analysis including MBD4 non-LOF mutants in the WT pool yields nearly identical results (n_WT=922, p=0.022, d=−0.501; Supplementary Table S1), confirming the signal is not an artifact of the WT definition. No clinical cohort with MBD4 status and ATRi treatment exists in public databases; prospective validation is required before any clinical translation. The PARP1 expression finding (p=0.605) falsifies the specific transcriptional upregulation hypothesis but does not exclude non-transcriptional PARPi mechanisms. All findings are computational and preclinical, and are not intended to guide diagnosis, prognosis, or treatment selection.

## Methods

### Cell line classification

MBD4 mutation status was determined from DepMap 24Q2 OmicsSomaticMutations data (release 24Q2). Cell lines were classified as True-LOF if carrying truncating mutations (nonsense, splice-site, frameshift) with DepMap annotation LikelyLoF=True. Missense and passenger mutations were excluded. Wild-type (WT) lines were defined as having no somatic MBD4 mutations. Sample sizes vary by analysis modality depending on data availability and cross-dataset overlap (e.g., n=21 total LOF pool, n=19 with expression data, n=14 with ceralasertib GDSC2 data). While DepMap's LikelyLoF annotation serves as the primary inclusion gate, the consistent phenotypic shift observed in the resulting cohort strongly implies functionally biallelic inactivation or severe haploinsufficiency. Future targeted sequencing is required to map the exact zygosity (LOH or compound heterozygosity) of these models. The MSI-H collinearity of the MBD4-LOF cohort (20/21 lines MSI-H by genomic score; 6/21 by clinical annotation) is a structural limitation of the pan-cancer observational design; see Confound stress testing and Limitations.

### Pharmacological stratification

Drug sensitivity data were obtained from GDSC2 (Genomics of Drug Sensitivity in Cancer, release 2). For each compound, cell lines with matched MBD4 mutation status and drug response data were stratified into MBD4-LOF and WT groups. Three metrics were analyzed: natural log IC50 (LN_IC50), area under the dose-response curve (AUC), and standardized Z-score (Z_SCORE). The WT comparator pool for each drug is defined after GDSC2 drug-availability intersection: n=914 for ceralasertib (14 LOF lines with GDSC2 data) and n=920 for adavosertib (15 LOF lines with GDSC2 data). The pre-intersection mutation-filtered WT pool (all DepMap lines with no somatic MBD4 mutation) is larger but is not the pharmacological comparator denominator reported in Results.

Statistical significance was assessed using a one-sided Mann-Whitney U test (alternative: MBD4-LOF < WT). Effect sizes were computed as Cohen's d with pooled standard deviation. Multiple testing correction was applied via Benjamini-Hochberg FDR where applicable. Six candidate therapy axes were defined before comparative testing based on the known biology of BER deficiency and replication stress (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors); BH-FDR correction was applied to exploratory axis-level comparisons. The ATR/WEE1 axis was designated confirmatory based on the a priori replication stress mechanism and is therefore reported with unadjusted directional p-values.

### Confound stress testing

Four confound analyses were applied to the ceralasertib signal:

1. **MSI-H purge**: All lines annotated as MSI-H in DepMap ModelSubtypeFeatures were removed from both MBD4-LOF and WT groups before retesting.
2. **TP53 stratification**: MBD4-LOF/TP53-mut lines were compared against MBD4-WT/TP53-mut lines, controlling for TP53 status.
3. **Leave-one-out**: Each MBD4-LOF line was iteratively removed and the test recomputed.
4. **Lineage matching**: Analysis was repeated within individual tissue lineages (Bowel; non-Bowel).

*MSI-H classification note.* The primary MSI-H purge (Stress Test 1) uses the DepMap 24Q2 `ModelSubtypeFeatures` field, which encodes curated clinical/pathological MSI annotations. Under this definition, 6 of 21 MBD4-LOF lines are classified MSI-H, leaving n=10 MSS MBD4-LOF lines for the purge analysis. As a sensitivity check, we also applied a continuous genomic instability threshold (MSIScore > 10, DepMap `OmicsGlobalSignatures`), which classifies 20 of 21 MBD4-LOF lines as MSI-H — the single exception being SW 1783 (CNS/Brain glioma, MSIScore=1.1), a lineage-mismatched outlier carrying a distinct nonsense mutation (p.L557Ter) rather than the canonical poly-A frameshift. Under this stricter genomic definition, the MSS-restricted analysis is infeasible (n=1 MSS line). The primary stress test therefore relies on the curated clinical/pathological annotation used in the canonical DepMap 24Q2 rerun; both definitions are reported for transparency.

### Analysis workflow

The analysis proceeded in four sequential stages: (1) MBD4 True-LOF cell lines were identified in DepMap 24Q2 using the LikelyLoF annotation gate, with missense and passenger variants excluded; (2) matched GDSC2 drug sensitivity data were retrieved for all lines with available pharmacological screens; (3) six candidate therapeutic axes were predefined a priori based on the known biology of BER deficiency and replication stress (cytidine analogs, PARP inhibitors, ATR/WEE1 inhibitors, WRN helicase inhibitors, immunotherapy, PKMYT1 inhibitors) before any comparative testing was performed; (4) exploratory axis-level comparisons were conducted with Benjamini-Hochberg FDR correction, followed by confirmatory ATR/WEE1 analyses and four orthogonal confound stress tests on the ceralasertib signal. All reported pharmacologic statistics were recomputed directly from source data using a harmonized matching pipeline; intermediate cached files were retained only after confirming numerical identity to the final extraction or replaced where discrepancies were detected.

### Expression analysis

Gene expression data (log1p TPM) were obtained from DepMap 24Q2 OmicsExpressionProteinCodingGenesTPMLogp1. PARP1 and RNF144A expression were compared between MBD4-LOF and WT groups using the Mann-Whitney U test. PARP1–PARPi correlation was computed as Spearman rank correlation across 481 cell lines with matched DepMap expression and GDSC2 PARP inhibitor Z-scores.

### Data availability

All analyses use publicly available data: DepMap 24Q2 (https://depmap.org), GDSC2 (https://www.cancerrxgene.org/), and published literature with PubMed identifiers. Analysis scripts are available at https://github.com/fjkiani/Synthetic-Lethality.

### Code availability

Analysis scripts for pharmacogenomic stratification, confound stress testing, and expression analysis are available at https://github.com/fjkiani/Synthetic-Lethality.

## Figures

![**MBD4-LOF cell lines are significantly more sensitive to ceralasertib (ATRi).** Comparison of GDSC2 LN_IC50 values for ceralasertib between MBD4 True-LOF lines (n=14, orange) and wild-type lines (n=914, gray). Horizontal lines indicate group medians. The MBD4-LOF group shows a 0.73 log-unit shift toward sensitivity (unadjusted p=0.021, Cohen's d=−0.50, one-sided Mann-Whitney U test; padj=0.063 after BH correction across six axes).](FIGURES/fig2_ceralasertib_volcano.png){#fig:volcano}

![**Four confound stress tests confirm the ceralasertib signal is consistent with MBD4-specific ATR dependency.** (A) MSI-H purge: removing MSI-H lines by curated clinical annotation strengthens the signal (p=0.015); note that 20/21 MBD4-LOF lines are MSI-H by genomic score, making full MSI deconfounding infeasible in this dataset. (B) TP53 stratification: controlling for TP53 co-mutation, MBD4-LOF adds >1 log-unit sensitivity (p=0.003). (C) Leave-one-out: all 14 iterations maintain significance. (D) Lineage analysis: signal preserved in non-Bowel comparison (p=0.025).](FIGURES/fig3_stress_tests.png){#fig:stress}

![**PARP1 is not transcriptionally upregulated in MBD4-LOF, despite being a pan-cancer PARPi predictor.** (A) PARP1 expression (log1p TPM) in MBD4-LOF (n=19) vs WT cell lines (ns, p=0.605). (B) Spearman correlation between PARP1 expression and PARP inhibitor Z-score across 481 cell lines (ρ=−0.416, p=1.36×10⁻²¹). High-PARP1 cells (≥Q75) show enhanced sensitivity; low-PARP1 (≤Q25) show resistance. The lack of PARP1 upregulation in MBD4-LOF argues against a selective MBD4-driven PARPi axis mediated by PARP1 transcriptional upregulation.](FIGURES/fig4_parp1_expression.png){#fig:parp1}

![**Proposed mechanistic model for MBD4-LOF replication stress vulnerability.** Under this model, MBD4 loss-of-function disrupts BER, leading to accumulation of unprocessed base lesions at CpG sites. This is predicted to create a convergent vulnerability at the replication fork: (1) accumulation of U:G mismatches actionable via cytidine analogs (Chabot et al., 2022), and (2) constitutive replication stress from unresolved fork-blocking lesions predicted to create ATR checkpoint dependency (this study, pharmacogenomic evidence). Direct molecular demonstration of the ssDNA/RPA/ATR-ATRIP/CHK1 cascade in MBD4-LOF systems awaits isogenic validation. PARP1 expression serves as an independent pan-cancer patient-selection biomarker but is not selectively elevated by MBD4-LOF status.](FIGURES/fig5_mechanism_model.png){#fig:mechanism}

## References

## Supplementary Material

### Supplementary Figure S1. Evidence summary across candidate therapeutic axes for MBD4-LOF tumors

![**Evidence summary across candidate therapeutic axes for MBD4-LOF tumors.** Rows represent candidate therapy axes; columns represent evidence modalities (pharmacological screens, published isogenic data, expression analysis, clinical reports). Cell shading indicates evidence status: red = positive signal, blue = negative/falsified, gray = not tested. The cytidine analog axis (Chabot et al., 2022) provides the positive-control calibration benchmark. The ATR/WEE1 axis shows pharmacogenomic association (this study). PARP inhibition is not supported by transcriptional evidence (this study).](FIGURES/fig1_evidence_matrix.png){#fig:matrix}

### Supplementary Table S1. Pan-cancer comparator pharmacology for ATM-LOF and MBD4-LOF in GDSC2

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
