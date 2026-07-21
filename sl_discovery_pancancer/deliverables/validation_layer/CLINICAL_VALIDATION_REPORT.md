# Trust, Benchmark & Clinical-Validation Layer — Pan-Cancer Synthetic-Lethality Portfolio

**Scope:** all 58 battery-surviving synthetic-lethality (SL) axes from the prior 5-layer confound-validated discovery phase.
**Question answered:** *how do we trust these findings, what is the benchmark, what is the proof, and how far can they be clinically validated?*
**Data vintage:** DepMap 24Q4 (integrated Chronos + per-screen `ScreenGeneEffect`); Broad Drug Repurposing Hub; clinpgx annotations (CREATED_2026-01-05). Engine commit `90a56b7`.

---

## 0. The honest ceiling (read this first)

CRISPR-derived SL candidates **cannot be clinically validated in this environment.** Clinical validation requires prospective patient or trial data that we do not have and did not fabricate. **Nothing in this portfolio is labeled "clinically proven."**

What this layer delivers instead is a **graded trust ladder** built entirely from real orthogonal datasets. Its top rung is *"reproducible, mechanistically coherent SL whose dependency node is a clinically actionable target with trial precedent"* — **not** *"clinically validated finding."* Novel candidates remain discovery-only nominations. Every axis carries an explicit tier and the reason it earned that tier. This distinction is enforced in every file and figure.

Three independent yardsticks were applied, each answering a different trust question:

| Rung | Question | Yardstick | Headline result |
|---|---|---|---|
| **A** | Does the pipeline recover SL we already believe? | Recall vs a 27-pair clinically-anchored truth set | **45.5% recall** on CRISPR-visible pairs; every miss mechanistically explained |
| **B** | Is the signal an artifact of the Chronos integration? | Re-test on the raw per-screen matrix | **53/53 replicate**, Spearman rho **0.973** — an *integration-robustness* check, not independent replication |
| **C** | Is the dependency node actionable, with what precedent? | Repurposing Hub + clinpgx + targeted literature | **13/58** MOA-consistent druggable; **17/58** druggable or with clinical-trial precedent |

---

## 1. Rung A — Known-SL recovery benchmark (precision/recall vs a clinical truth set)

**Truth set:** 27 SL pairs hand-curated from the primary citations already in `known_sl_panel.py`, each tagged with a real clinical anchor and an anchor tier (`approved` = 4, `pivotal_trial` = 13, `preclinical_validated` = 10). No new unsourced pairs were invented. A pair is scored "recovered" if **any** valid battery route survived (handles copy-number duplicate rows).

### Primary metric
- **CRISPR-visible recall = 10/22 = 45.5%** (the 5 PARP-trapping pairs are excluded here — see caveat below).
- **By anchor tier:** approved 1/2 (50%), pivotal_trial 6/11 (55%), preclinical_validated 3/9 (33%).

**Recovered (10):** MTAP→PRMT5 (FDR 1.7e-8), PTEN→PIK3CB (2.0e-5), ARID1A→ARID1B (8.8e-5), MTAP→MAT2A (5.6e-4), VHL→EPAS1 (4.4e-3), ARID1A→EZH2 (0.014), SMARCA4→SMARCA2 (0.014), KRAS→STK33 (0.025), CCNE1→PKMYT1 (0.037), CCNE1→WEE1 (0.037).

### Every miss is mechanistically explained (not a black-box failure)
The 12 missed pairs fall into exactly two known boundaries of CRISPR-KO viability screening:

1. **Broadly-essential partners that the pan-essential guard correctly rejects.** RB1→SKP2 is the instructive case: it reached FDR 3.6e-9 but was **guard-rejected because SKP2 is pan-essential** — i.e., the battery deliberately suppressed what would otherwise be a false headline. Also RB1→CDK4, KRAS→PTPN11, STK11→MTOR, ATM→ATR, CCNE1→ATR (ATR/MTOR broadly essential).
2. **Drug-stress / checkpoint-dependent SL that unperturbed CRISPR KO does not model.** TP53→WEE1, TP53→PKMYT1, TP53→CHEK1 are synthetic-lethal *under replication stress / chemotherapy*, not in baseline viability screens; NF1→MAP2K1 and SMARCA4→EZH2 are lineage-diluted in a pan-cancer analysis.

**What this bounds:** the pipeline is a **high-precision detector of CRISPR-tractable SL** — paralog buffering (ARID1B, SMARCA2), metabolic (PRMT5/MAT2A), and lineage/pathway addiction — and is **not** a detector of drug-mechanism-specific or broadly-essential SL. That is a feature (it refuses to over-call pan-essential genes), stated as a limitation of coverage.

### Two honest caveats on this benchmark
- **PARP/BRCA is modality-limited and reported separately, not as failure.** CRISPR knockout removes PARP1 protein; it does **not** reproduce the *trapping* mechanism that makes PARP inhibitors lethal in HR-deficient cells. Of 5 PARP-trapping pairs, only BRCA1→PARP1 reached the survival gate (FDR 0.25), and only under a relaxed threshold. These pairs are quarantined from the CRISPR recall number by design.
- **Precision against the truth set is not meaningfully computable.** Most pipeline positives are *novel* pairs absent from any curated truth set, so a naive precision calculation would penalize genuine discovery. Only **recall** and **reclassification-correctness** (does the battery correctly demote known false headlines like RB1→SKP2?) are interpretable here.

---

## 2. Rung B — Cross-screen reproducibility (integration-robustness, NOT independent replication)

> **Framing gate (critical):** the "independent" per-screen matrix (`ScreenGeneEffect`, 1186 models after ScreenID→ModelID collapse) shares **1177/1186 cell lines with the integrated Chronos discovery matrix — only 9 are truly held out.** Therefore this is an **integration-robustness check** (raw per-screen effect vs integrated Chronos model), **not** independent-cohort replication. A survivor replicating here confirms its signal is **not an artifact of the Chronos integration pipeline**; it does **not** prove replication in an independent patient or cell-line cohort. This is stated everywhere the number appears.

The L1 differential-dependency test (one-sided Mann–Whitney, Cohen's d, BH-FDR) was re-run on the per-screen matrix for all 58 survivors.

- **Testable = 53** (4 copy-number axes — MTAP→PRMT5/MAT2A, CCNE1→WEE1/PKMYT1 — are `na_cn`: CN calls are absent from the per-screen matrix, so they retain their expression-proxy CN detection from discovery; 1 axis, NRAS→N6AMT1, is `partner_absent` because N6AMT1 is not in the per-screen matrix).
- **Replication rate among testable = 53/53 = 100%; sign concordance = 100%.**
- **Effect-size concordance: Spearman rho = 0.973 (p = 4.8e-34), Pearson r = 0.990.**
- **Positive controls all replicate:** PTEN→PIK3CB (p=7.4e-6, d=−0.79), BRAF→MAPK1 (p=6.6e-17, d=−1.79), BRAF→MAP2K1 (p=2.3e-12, d=−1.49), STAG2→STAG1 (p=4.0e-7, d=−3.08), VHL→EPAS1 (p=1.1e-3, d=−2.07).

**Interpretation:** 100% replication is *expected* given the cell-line overlap; its value is as a **negative-control that did not fail** — none of the 58 survivors is a Chronos-integration artifact. A drop here would have been a red flag. It is a **necessary robustness condition, not sufficient clinical proof.**

---

## 3. Rung C — Clinical actionability + precedent (real datasets, with citations)

Actionability is keyed on the **partner** (the CRISPR dependency = the druggable node). Three sources were joined:

**Broad Drug Repurposing Hub** (6,798 drugs; 2,183 target-genes indexed). Raw hub matching is **promiscuous/noisy** — e.g., BRAF→MAPK1's top raw drug is "arsenic-trioxide," KRAS→CTNNB1 maps to "urea." To avoid over-claiming, a **MOA-consistency filter** was applied: a drug counts only if its MOA string contains an inhibitor/antagonist/degrader/blocker/agonist/modulator term **and** its annotated target set is ≤8 genes (selective). Promiscuous-only nodes are flagged, not counted.

- **13/58 clinically actionable** (MOA-consistent selective drug, any phase).
- **9/58 actionable at late phase** (Phase 3+/approved).
- MOA-correct launched examples: BRAF→MAP2K1 (trametinib/cobimetinib/binimetinib, MEK), KRAS/NRAS→RAF1 (dabrafenib/vemurafenib), PTEN→PIK3CB (alpelisib/copanlisib), CDKN2A→EGFR (erlotinib/osimertinib), ARID1A→EZH2 (tazemetostat), STK11→HDAC4 (panobinostat). Clinical-stage: CCNE1→WEE1 (adavosertib/MK-1775), VHL→EPAS1 (PT-2385).

**clinpgx** (127,601 gene-drug-phenotype relationships; drugLabels.byGene). **23/58** survivor partners carry a clinpgx annotation; **5/58** have an FDA/EMA label (e.g., EGFR with 148 annotations + label; PIK3CB; WRN with PMIDs 28791697/28796378).

**Targeted literature precedent** (12 headline axes searched via LiteratureSearch; the remaining 46 axes are marked `not_searched`, which means *unexamined*, **not** "no precedent"). This closed the gap for nodes the hub misses because the drug is too new to be indexed:
- **WRN family — HRO761 (Novartis) Phase 1 NCT05838768 ± irinotecan/tislelizumab [5]; VVD-214/RO7589831 (Vividion/Roche) clinical-stage [1,4]**.
- **STAG2→STAG1** — established cohesin paralog SL (van der Lelij 2017 [13]; Benedetti 2017 [15]); STAG1 degraders in development [14].

Combining hub + literature, **17/58** axes have a druggable node **or** a clinical-trial precedent (`external_clinical_target` captures WRN/STAG1 via literature).

---

## 4. Trust-tier assignment (added to the master matrix)

Each survivor's tier combines all three rungs. CN axes (`na_cn`) are treated as replicated-effective for tiering (flagged `repro_basis = discovery_only_cn_not_retestable`); NRAS→N6AMT1 is `partner_absent_in_screen`.

| Tier | Definition | n | Members |
|---|---|---|---|
| **T1 — clinically anchored** | reproduced + late-phase drug (rank ≥3) or trial precedent + clinical annotation/precedent | **12** | BRAF→MAPK1, BRAF→MAP2K1, KRAS→RAF1, NRAS→RAF1, STAG2→STAG1, ARID1A→WRN, CDKN2A→EGFR, PTEN→PIK3CB, KMT2D→WRN, STK11→HDAC4, EP300→WRN, ARID1A→EZH2 |
| **T2 — reproduced + druggable** | reproduced + selective drug (any phase), no strong clinical precedent | **5** | CDKN2A→TUBB4B, PTEN→BRD2, VHL→EPAS1, CCNE1→WEE1, BRCA1→PARP1 |
| **T3 — reproduced, not yet druggable** | reproduced, no MOA-consistent clinical-stage drug for the node | **40** | e.g. KRAS→CTNNB1, KRAS→TCF7L2, MTAP→PRMT5, RB1→E2F3, NRAS→SHOC2, EP300→CREBBP … |
| **T4 — discovery-only** | survives battery but partner absent from the reproducibility screen | **1** | NRAS→N6AMT1 |

`drug_source` breakdown: none = 33, repurposing_hub = 20, literature = 5. Every tier records which rungs passed/failed; no axis was deleted, consistent with the prior phase. **The 58 survivors and 44 reclassified axes from the prior matrix are preserved** — this phase only *added* columns and tiers (integrity asserted in the receipts).

---

## 5. Emergent hypothesis surfaced by the trust layer: chromatin-remodeler loss → WRN dependency (MSI-independent angle)

Three of the twelve T1 axes converge on a single dependency node — **WRN** — driven by loss of distinct chromatin regulators: **ARID1A→WRN** (FDR 6.6e-6), **KMT2D→WRN** (2.7e-5), **EP300→WRN** (5.3e-5). All three replicate in the robustness screen.

- The canonical WRN SL is **MSI-H-specific** (Chan 2019 [11]; Picco 2021 [12]), and clinical WRN inhibitors show diminished effect in MSS tumors [2].
- Our axes survived the prior phase's **MSI purge**, i.e., they persist in an MSI-adjusted analysis — pointing to a **chromatin-loss / MSI-independent** route to WRN dependency.
- Independent groups now report exactly this: **ARID1A-mutant WRN dependency shown in vivo via a distinct DDR mechanism** (Science Advances 2026 [10]; Cancer Research 2025 [7]).
- The node is **clinically targeted today** (HRO761 Phase 1; VVD-214) [4,5].

This is strong external corroboration of a *mechanistic direction*, not clinical proof of our specific axes. The KMT2D→WRN and EP300→WRN links are less established in the literature than ARID1A→WRN and remain discovery-grade nominations with a clinically targeted node.

Note on **STAG2→STAG1** (our largest effect, d=−3.08): real and clinically interesting, but **literature-known cohesin paralog SL** — it is novel-to-this-pipeline, **not** novel-to-the-field, and is labeled as such.

---

## 6. Limitations (non-negotiable)

1. **Orthogonal corroboration ≠ clinical proof.** No axis is clinically validated. The top tier means "reproducible + clinically actionable node with precedent," nothing more.
2. **Benchmark B is integration-robustness, not independent replication** (1177/1186 cell-line overlap). Truly independent replication needs a non-DepMap CRISPR cohort or patient data.
3. **46/58 axes were not literature-searched** (`not_searched` ≠ "no precedent"). Absence of a precedent annotation is absence of evidence, not evidence of absence.
4. **4 copy-number axes are not re-testable** in the per-screen matrix and retain expression-proxy CN detection.
5. **PARP/BRCA remains modality-limited** — CRISPR KO cannot model PARP trapping.
6. **Hub target annotations are promiscuous;** the MOA-consistency filter mitigates but does not eliminate mis-annotation. Drug lists are actionability leads, not treatment recommendations.
7. **No efficacy is claimed anywhere.** T1/T2 labels trace to a repurposing-hub phase and/or clinpgx PMID; none asserts patient efficacy.
8. **Novel axes stay discovery-only** until externally validated in an independent cohort.

---

## 7. Deliverables

| File | Contents |
|---|---|
| `CLINICAL_VALIDATION_REPORT.md` | this report |
| `sl_master_matrix.csv` / `.json` | prior matrix extended (24→50 cols); 58 survivors + 44 reclassified preserved |
| `clinical_validation.csv` | per-axis full validation ladder (all three rungs + tier) |
| `benchmark_known_sl_recovery.json` | Rung A: recall vs clinical truth set, misses with reasons, PARP caveat |
| `benchmark_cross_screen.csv` / `_summary.json` | Rung B: per-axis robustness replication |
| `clinical_actionability.csv` | Rung C: hub + clinpgx actionability per axis |
| `fig_benchmark_scorecard.png/.svg` | recall by anchor tier + replication counts |
| `fig_cross_screen_concordance.png/.svg` | discovery vs per-screen effect-size scatter (rho 0.973) |
| `fig_trust_tier_landscape.png/.svg` | tiers by stream + T1/T2 significance-with-drug |
| `fig_clinical_actionability_map.png/.svg` | 17 druggable/precedent axes by clinical phase |
| `SL_VALIDATION_RECEIPTS.json` | provenance + SHA256 + dataset versions + engine commit |

*Citations `[N]` refer to the literature sources retrieved during this analysis.*
