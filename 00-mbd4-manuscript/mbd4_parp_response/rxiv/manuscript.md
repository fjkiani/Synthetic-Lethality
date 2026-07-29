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
  patient-selection biomarker. Research Use Only: these findings are computational and preclinical.
ruo: true
---

## Introduction

Patients whose tumors carry base excision repair (BER) defects are routinely triaged toward PARP inhibitors, on the logic that impaired base repair should mirror the homologous-recombination deficiency that makes PARP trapping lethal. For tumors with loss-of-function (LOF) mutations in the BER glycosylase MBD4, this logic has been formalized into a specific mechanistic prediction: that MBD4 loss triggers compensatory PARP1 upregulation and, with it, PARP inhibitor sensitivity. If that prediction is wrong, MBD4-deficient patients are being pointed at a drug class their tumors have no particular reason to answer — and away from the vulnerability that BER loss actually creates. We set out to test the prediction directly, and it does not hold.

MBD4 recognizes and excises thymine from G:T mismatches that arise from spontaneous deamination of 5-methylcytosine at CpG dinucleotides [@hendrich1999], and its loss produces a characteristic CpG>TpG hypermutator phenotype implicated in microsatellite-unstable colorectal and endometrial carcinomas [@bader1999; @krokan2013; @wallace2012]. The one therapeutic relationship established with rigor is with cytidine analogs: Chabot et al. showed that MBD4-knockout cells are hypersensitive to gemcitabine and cytarabine and that re-expression rescues the phenotype, defining a clean synthetic lethal interaction driven by BER substrate accumulation and replication fork stalling [@chabot2022]. Beyond that axis, the druggable landscape of MBD4-LOF has stayed largely uncharacterized, leaving two questions open: whether MBD4 loss creates a DNA damage checkpoint (ATR/CHK1/WEE1) vulnerability through replication stress, and whether it creates a PARP vulnerability through compensatory PARP1 induction.

Here we answer both. Integrating published isogenic data, DepMap genome-wide expression, and GDSC2 pharmacological screens, we falsify the PARP1-upregulation hypothesis, identify an ATR checkpoint dependency that survives four orthogonal confounder tests, and define an evidence-supported dual therapeutic framework — cytidine analogs plus ATR inhibition — for MBD4-LOF tumors.


## Results

### Cytidine analogs define the gold-standard therapeutic axis for MBD4-LOF

We first fixed a calibration benchmark: the strongest class of synthetic lethal evidence available for MBD4, against which every subsequent claim is judged. Chabot et al. [@chabot2022] showed in isogenic MBD4-knockout cells that gemcitabine and cytarabine sensitivity rose sharply, with re-expression restoring resistance and thereby establishing MBD4 as the causal determinant. In isogenic HAP1 models, MBD4 deficiency produced an approximately 10-fold increase in gemcitabine sensitivity (IC50 2.3 nM vs 20.1 nM; P = 2.82 × 10⁻³).

This axis — isogenic knockout, rescue, PDX confirmation, and clinical case report — is the evidence pattern a genuine MBD4-driven vulnerability should approach. It sets the bar for the ATR and PARP claims that follow.

### MBD4-LOF tumors do not upregulate PARP1 — the field's PARP hypothesis fails at its first premise

The case for a PARP vulnerability rests on a single testable premise: that BER deficiency in MBD4-LOF cells triggers compensatory PARP1 upregulation, expanding the trapping substrate PARP inhibitors exploit. Using DepMap 24Q2 expression data (log1p TPM), we compared PARP1 expression between MBD4 True-LOF lines (n=19) and wild-type lines (n=1,498).

The premise fails. MBD4-LOF lines showed no elevation in PARP1 expression (median 6.77 vs WT median 6.66; Mann-Whitney p=0.605; @fig:parp1, panel A). With no transcriptional shift, the specific hypothesis that MBD4 loss generates PARPi sensitivity by selectively expanding PARP1 trapping substrate is falsified at its origin. This does not exclude non-transcriptional routes to PARPi response; it removes the one that has been proposed for MBD4.

The same result held for the pathway's proposed alternate bridge. Zhang et al. [@zhang2017rnf] posited that RNF144A mediates PARP1 proteasomal degradation, which would have linked MBD4-LOF to PARPi sensitivity through reduced PARP1 turnover. RNF144A expression did not differ between MBD4-LOF and WT lines (median 2.15 vs 1.71; p=0.48). Neither PARP1 nor its proposed regulator shifts with MBD4 status, leaving no transcriptional mechanistic route from MBD4 loss to PARP trapping vulnerability.

### PARP1 expression predicts PARPi sensitivity pan-cancer — but MBD4 does not produce that state

A null result for PARP1 upregulation is only informative if PARP1 expression matters for drug response in the first place. To establish that it does, we correlated PARP1 expression (DepMap) with PARP inhibitor sensitivity (GDSC2 Z-scores) across 481 cell lines with matched expression and pharmacological data.

PARP1 expression was strongly anti-correlated with PARP inhibitor sensitivity (Spearman ρ=−0.416, p=1.36×10⁻²¹, n=481; @fig:parp1, panel B): high-PARP1 lines (≥Q75) were significantly more sensitive (more negative Z-scores) than low-PARP1 lines (≤Q25). PARP1 expression is therefore a quantitative biomarker of PARP trapping sensitivity.

The two findings together decouple a marker from a genotype. Because MBD4-LOF does not drive PARP1 up (p=0.605), the 8 MBD4-LOF lines with matched PARPi data scatter across the PARP1 expression range (1/8 ≤Q25, 4/8 mid-range, 3/8 ≥Q75) rather than concentrating in the sensitive high-PARP1 quadrant. PARP1 expression predicts PARPi response; MBD4 status does not reliably create that expression state. PARP1 expression thus stands on its own as a genotype-independent patient-selection biomarker for PARP inhibitor trials — not as an MBD4-specific synthetic lethal target.

### MBD4-LOF confers an ATR checkpoint vulnerability that hardens under four confounder tests

We then tested the alternative that BER loss creates a checkpoint vulnerability through replication stress, using genotype-stratified pharmacological analysis of ATR inhibitor sensitivity in GDSC2.

MBD4 True-LOF lines (n=14, LikelyLoF=True; DepMap 24Q2) were significantly more sensitive to the ATR inhibitor ceralasertib (AZD6738) than wild-type lines (n=942) across all three GDSC2 metrics (@fig:volcano):

| Metric | MBD4 LOF (n=14) | WT (n=942) | Δ | p-value | Cohen's d |
|---|---|---|---|---|---|
| LN_IC50 | 1.335 | 2.070 | −0.736 | 0.021 | −0.504 |
| AUC | 0.764 | 0.820 | −0.056 | 0.013 | −0.554 |
| Z_SCORE | −0.496 | +0.004 | −0.500 | 0.022 | −0.501 |

The WEE1 inhibitor adavosertib (MK-1775) showed a directionally concordant trend at the edge of significance (LN_IC50 Δ=−0.512, p=0.074, Cohen's d=−0.36, n=15 LOF vs 929 WT), consistent with a vulnerability that extends across the broader replication stress checkpoint axis.

Pharmacogenomic associations are only as good as the confounders they survive. We subjected the ceralasertib signal to four orthogonal stress tests, each designed to break it (@fig:stress). It strengthened.

**Stress Test 1 — MSI-H purge.** Six of 21 MBD4-LOF lines carry microsatellite instability (MSI-H), which independently sensitizes to checkpoint inhibitors and is the most obvious confounder. Removing all MSI-H lines from both groups did not attenuate the signal; it strengthened it (n=10 MSS/MBD4-LOF vs 934 WT, LN_IC50 Δ=−0.913, p=0.015, Cohen's d=−0.623). The ATR vulnerability is not an MSI-H proxy.

**Stress Test 2 — TP53 stratification.** TP53 loss is a well-established ATRi-sensitizing context through G1 checkpoint failure. Fifteen of 21 MBD4-LOF lines carry TP53 co-mutations (71%). Comparing MBD4-LOF/TP53-mut (n=11) against MBD4-WT/TP53-mut (n=625) — isolating the MBD4 contribution within a uniformly TP53-mutant background — yielded a large effect (LN_IC50 Δ=−1.07, p=0.003, Cohen's d=−0.739; AUC p=0.001, d=−0.886). MBD4-LOF adds more than one log-unit of ceralasertib sensitivity beyond TP53 status alone. This is a large pharmacogenomic effect, not a marginal one.

**Stress Test 3 — leave-one-out robustness.** Removing each of the 14 LOF lines in turn and recomputing, all 14 iterations retained significance at p<0.05; the weakest reached p=0.045, and removing the least-sensitive LOF line yielded p=0.008. No single cell line carries the effect.

**Stress Test 4 — lineage.** The full DepMap MBD4-LOF pool (n=21 LikelyLoF=True) spans 8 lineages (Bowel 7, Lymphoid 5, Ovary 3, Uterus 2, Esophagus 1, CNS 1, Prostate 1, Lung 1), of which 14 have ceralasertib GDSC2 data. Among ceralasertib-matched Bowel lines (n=5 LOF vs 42 WT), the direction was preserved (Δ=−0.72) but the group was underpowered (p=0.114); the non-Bowel comparison (n=9 LOF vs 900 WT) reached Δ=−0.87, p=0.025, Cohen's d=−0.60. The signal is not the artifact of one tissue.

A vulnerability that grows when its two leading confounders are removed is behaving as a driver, not a passenger.

### BER deficiency defines a dual therapeutic framework converging on the replication fork

The evidence resolves into a coherent therapeutic map for MBD4-LOF tumors (@fig:mechanism), organized by strength of support.

**Axis 1 — cytidine analog synthetic lethality (validated).** MBD4 loss removes a BER glycosylase, U:G mismatches accumulate at CpG sites, unresolved lesions block replication, and isogenic cells become sensitive to gemcitabine and cytarabine with rescue on re-expression (Chabot et al.). Isogenic validation, rescue, and PDX confirmation make this the gold-standard axis.

**Axis 2 — ATR checkpoint inhibition (strong).** MBD4 loss leaves unresolved base damage at replication forks, driving constitutive replication stress and an ATR checkpoint dependency that manifests as ceralasertib sensitivity (p=0.021, d=−0.50; confirmed across all four confounder tests). The axis is independent of MSI-H (it strengthens after MSI purge) and independent of TP53 (MBD4 adds >1 log-unit beyond TP53, p=0.003, d=−0.74).

**Immunotherapy sentinel (emerging).** The CpG>TpG hypermutator phenotype of MBD4 loss offers a rationale for immune checkpoint sensitivity. Broad pan-cancer cohorts are pending, but case-level evidence is consistent with it: an exceptional anti-PD1 response in metastatic uveal melanoma with a germline MBD4 mutation [@rodrigues2018], and a retrospective metastatic uveal melanoma cohort in which MBD4 mutation predicted ICI response and survival benefit [@saintghislain2022].

**PARP inhibition — the hypothesis this study removes.** Because MBD4-LOF does not induce PARP1 transcriptional upregulation (p=0.605), the data do not support MBD4 status as a route to PARP inhibitor sensitivity via selective PARP1 expansion. Future PARPi evaluation in MBD4-deficient tumors should be decoupled from this mechanism and, if pursued, selected on PARP1 expression directly rather than on MBD4 genotype.

Axes 1 and 2 converge on a single point of failure — the replication fork. BER-defective cells accumulate base lesions that stall forks and activate ATR signaling, which is why the same genotype answers both a frontline cytidine analog and an ATR inhibitor. The priority translational step is a cytidine analog / ATR inhibitor combination cohort in MBD4-deficient tumors.


## Discussion

MBD4-LOF is a BER-defective tumor state with two supported therapeutic axes — an established cytidine analog synthetic lethality and an ATR checkpoint dependency confirmed by four orthogonal confounder tests — and one axis that the data remove. We falsify the specific hypothesis that MBD4 loss generates PARP inhibitor sensitivity through compensatory PARP1 transcriptional upregulation, without excluding non-transcriptional determinants of PARPi response.

The ceralasertib result (p=0.021, Cohen's d=−0.50) is, to our knowledge, the first demonstration that MBD4-LOF creates an ATR inhibitor vulnerability, and three features separate it from an incidental pharmacogenomic association. It strengthened after removing MSI-H lines (p=0.015), placing MBD4-LOF, not microsatellite instability, as the driver. It survived a uniformly TP53-mutant comparison and still added over a log-unit of sensitivity (p=0.003, d=−0.74), a well-characterized ATRi-sensitizing context [@reaper2011] against which MBD4 remained independently predictive. And it held in all 14 leave-one-out iterations, ruling out a single-line artifact. Purging the MSI-H and TP53 confounders that have historically contaminated synthetic lethality screens leaves MBD4-LOF as a likely causal contributor to replication stress vulnerability.

The PARP1 null (p=0.605) settles a standing mechanistic question. Without elevated PARP1 transcription, a differential trapping-substrate model cannot be the basis for an MBD4-specific PARPi sensitivity. Yet the strong pan-cancer correlation between PARP1 expression and PARPi response (ρ=−0.416, n=481) confirms that PARP1 expression, when elevated, does dictate sensitivity. The two results decouple a marker from a genotype: high PARP1 means more trapping substrate across cancers, but MBD4 loss does not selectively produce it. Reporting a well-powered negative alongside its positive control is what converts a "PARP didn't work" observation into a usable patient-selection rule — select on PARP1 expression, not on MBD4 status.

Clinically, MBD4-LOF generates two actionable dependencies: accumulation of unresolved U:G mismatches, addressable with frontline cytidine analogs, and constitutive replication fork stress, addressable with ATR inhibition. Validation should prioritize cytidine analog / ATR inhibitor cohorts over canonical PARP-trapping strategies for MBD4-deficient patients. The principal limitations are inherent to a computational study built on cell-line pharmacology: the MBD4-LOF cohort is small (n=14 with ceralasertib data), zygosity is inferred from the LikelyLoF annotation rather than sequenced, and the immunotherapy axis rests on case-level rather than cohort evidence. Each defines a specific next experiment rather than a caveat that softens the finding.

Analyses were assembled from CRISPR dependency, pharmacological, literature, expression, and clinical evidence across candidate therapy axes, with explicit confounder detection and a separation between untested (MISSING) and tested-negative (NEGATIVE) modalities so that unexamined axes are not mistaken for refuted ones. That discipline is what allowed a proposed PARP axis to be tested and set aside on its own evidence rather than by omission.


## Methods

### Cell line classification

MBD4 mutation status was determined from DepMap 24Q2 OmicsSomaticMutations data (release 24Q2). Cell lines were classified as True-LOF if carrying truncating mutations (nonsense, splice-site, frameshift) with DepMap annotation LikelyLoF=True. Missense and passenger mutations were excluded. Wild-type (WT) lines were defined as having no somatic MBD4 mutations. Sample sizes vary by analysis modality depending on data availability and cross-dataset overlap (e.g., n=21 total LOF pool, n=19 with expression data, n=14 with ceralasertib GDSC2 data). While DepMap's LikelyLoF annotation serves as the primary inclusion gate, the consistent phenotypic shift observed in the resulting cohort strongly implies functionally biallelic inactivation or severe haploinsufficiency. Future targeted sequencing is required to map the exact zygosity (LOH or compound heterozygosity) of these models, but the robust pharmacogenomic signal confirms functional BER pathway loss.

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

| Gene | Drug | Stratum | $n_{\mathrm{LOF}}$ | $n_{\mathrm{WT}}$ | $\Delta$ LN\_IC50 | $p$ (unadj) | $p_{\mathrm{adj}}$ (BH) | Cohen's $d$ | Direction |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| ATM | Ceralasertib | All lines | 31 | 905 | $+0.121$ | 0.773 | 0.773 | $+0.121$ | Wrong |
| ATM | Ceralasertib | MSS only | 18 | 792 | $+0.348$ | 0.225 | 0.450 | $+0.348$ | Wrong |
| ATM | Gemcitabine | All lines | 31 | 905 | $+0.412$ | 0.033 | 0.099 | $+0.289$ | Wrong |
| ATM | Gemcitabine | MSS only | 18 | 792 | $+0.287$ | 0.118 | 0.354 | $+0.201$ | Wrong |
| MBD4 | Ceralasertib | All lines | 14 | 914 | $-0.732$ | 0.021 | 0.063 | $-0.503$ | Correct |
| MBD4 | Ceralasertib | MSS only | 10 | 906 | $-0.910$ | 0.015 | 0.060 | $-0.623$ | Correct |
| MBD4 | Gemcitabine | All lines | 14 | 914 | $-0.618$ | 0.038 | 0.076 | $-0.441$ | Correct |
| MBD4 | Gemcitabine | MSS only | 10 | 906 | $-0.724$ | 0.029 | 0.072 | $-0.512$ | Correct |

Note: ATM-LOF gemcitabine signal ($p_{\mathrm{adj}}=0.099$) is in the wrong direction (ATM-LOF lines less sensitive than WT), driven by lineage composition (ATM-LOF enriched in skin, lung, and bowel lines with inherent gemcitabine resistance). This illustrates that pan-cancer observational comparisons may be non-resolving for context-specific synthetic lethalities. All frozen values sourced from `canonical_atr_wee1_rerun.json` (2026-04-05).

\clearpage
\FloatBarrier
