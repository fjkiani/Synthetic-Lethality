# Clinical-Grade Validation of the CrisPRO Synthetic-Lethality Layer

**Schema:** `v2.0-clinical-grade` · **Scope:** all 58 SL survivor axes · **Discovery cohort:** DepMap 24Q4 (Broad Avana + Chronos, 1,178 lines) · **Independent cohort:** Sanger Project Score CERES (Figshare 9116732; Behan 2019, *Nature*)

---

## 0. Why this report does **not** say "clinically proven"

The request was to make the validation layer "clinical proof" and "clinical grade." One thing is stated up front, in writing, because it governs every claim below:

> **We decline to label any axis "clinically proven," "clinically validated," or "validated in patients."** That evidentiary bar requires **prospective human outcome data for the specific driver→partner relationship** — a randomized or single-arm trial in which disrupting the partner produces clinical benefit in patients selected by the driver genotype. **No such data exist for these novel axes,** and fabricating or implying them would be scientific misconduct. For the handful of axes where a matched drug is already approved (e.g. PARP inhibitors in *BRCA1*-mutant disease), the *clinical proof* belongs to that drug programme, not to our computational nomination.

What we *can* deliver — and what "clinical grade" honestly means here — is a **graded, citation-traceable evidence ladder** plus **genuine independent replication on an orthogonal experimental platform**. That is the strongest defensible position, and it is what follows. Where a caveat is only *partially* closed, we say so explicitly rather than round up.

### Status of the four original caveats

| # | Original caveat | Status after this phase | What closed it / what remains |
|---|---|---|---|
| 1 | "Orthogonal corroboration is not clinical proof." | **Reframed, not overturned.** | We keep the honest label. Corroboration is now formalised as an evidence ladder (E1–E5); the ceiling is precedent/actionability, never efficacy. |
| 2 | "Prior Benchmark B is robustness, not independent replication." | **Partially closed.** | Replaced the same-pipeline robustness check (ρ=0.973, reused Chronos) with **cross-platform replication** on the Sanger KY/CERES platform (ρ=0.791). This is orthogonal library + wet-lab + algorithm — strictly stronger — but the cell panels overlap (315/318), so it is **not** out-of-cohort replication. Stated everywhere. |
| 3 | "46/58 axes were not literature-searched." | **Closed.** | All **58/58** axes literature-searched; every axis carries `lit_searched=True`, a `precedent_class`, and citation indices. "None found" now means *searched and nothing found*, a real result. |
| 4 | "Novel axes remain discovery-only nominations." | **Closed for grading; honestly retained where true.** | Every axis now has an explicit evidence grade. 28/58 have external biological precedent (E1–E4); 30/58 remain E5 = computational nominations — and are **labeled as such**, not upgraded. |

---

## 1. Independent replication (cross-platform)

### 1.1 What was done

Every survivor axis was re-tested with the **identical** L1 differential-dependency statistic used in discovery — one-sided Mann–Whitney (partner more depleted in driver-mutant lines), Cohen's *d*, Benjamini–Hochberg FDR across the tested set — but on a **fully independent experimental platform**:

| | Discovery | Independent replication |
|---|---|---|
| Institute | Broad | **Sanger** |
| Library | Avana | **KY1.0/1.1** (Tzelepis 2016) |
| Wet-lab / screen | Broad | **Sanger Project Score** (Behan 2019) |
| Dependency algorithm | Chronos | **CERES** |

This differs from the prior Benchmark B, which re-used the Chronos pipeline on integrated data and therefore measured *integration robustness*, not replication.

### 1.2 The honest scope limitation (stated plainly)

The verifiable Figshare CERES release (`10.6084/m9.figshare.9116732`) was **pre-filtered by its authors** to Sanger lines that also have Broad copy-number data, and cross-dataset normalised. Consequently:

- **315 / 318** Sanger lines **overlap** the discovery panel (99.1%).
- **3** lines are held-out (Sanger-only) — and **none** carry mutation calls.
- **→ Out-of-cohort (held-out) replication rate = N/A** (reported, not hidden).
- **311** shared lines carry mutation calls and form the testable universe.

A bounded best-effort hunt for an unfiltered, non-overlapping Sanger matrix (Score portal, Cell Model Passports S3) **failed** — the portal is a JavaScript app returning its HTML shell for guessed data paths, and the object store rejects unauthenticated keys. The Figshare CERES file is the only receipt-backed Sanger source, so this is **cross-platform** replication (orthogonal technology on shared lines), **not cross-cohort**. See `fig_replication_overlap.png`.

### 1.3 Results

- **51 / 58** axes testable on the Sanger platform (4 are copy-number-driven and not re-derivable in this matrix; 3 have < 5 driver-mutant lines in Sanger → `insufficient_n`; 0 partners absent from the KY library).
- **47 / 51** sign-concordant with discovery.
- **41 / 51 replicated** = sign-concordant **and** FDR < 0.05 on the orthogonal platform.
- **Effect-size concordance (n = 51 tested):** **Spearman ρ = 0.791** (p = 5.1×10⁻¹²), **Pearson r = 0.813** (p = 4.5×10⁻¹³).

This ρ = 0.791 is the **honest orthogonal-technology number**. It is lower than the prior integration-robustness ρ = 0.973 *by design* — a different library, wet-lab, and algorithm introduce real technical variance, and the effect survives it. See `fig_independent_replication.png`.

### 1.4 Positive-control gate (must pass)

| Positive control | Sanger status | Cohen's *d* | FDR | Call |
|---|---|---|---|---|
| PTEN → PIK3CB | tested (n=27) | −0.75 | 2×10⁻³ | **replicated** |
| BRAF → MAPK1 | tested (n=18) | −1.34 | 2.6×10⁻⁵ | **replicated** |
| BRAF → MAP2K1 | tested (n=18) | −1.84 | 2.6×10⁻⁵ | **replicated** |
| STAG2 → STAG1 | tested (n=9) | −1.93 | 4.5×10⁻² | **replicated** |
| VHL → EPAS1 | **insufficient_n** (n=3) | — | — | not counted (direction concordant) |

Four of five canonical SL controls replicate on the orthogonal platform with sign preserved and FDR < 0.05. VHL → EPAS1 has only 3 VHL-mutant lines in the Sanger panel; it is correctly flagged `insufficient_n` — **not** a failure, and its effect direction remains concordant. Gate **PASSED**.

### 1.5 Informative non-replications (the metric discriminates)

- **KRAS → STK33** — sign **flips** on the orthogonal platform (Sanger *d* = +0.02 vs discovery −0.29). This recapitulates the historical record: Scholl 2009 claimed STK33 synthetic lethality with mutant KRAS; Babij 2011 could not reproduce it. The cross-platform test independently flags it. E5.
- **TP53 → KRTAP4-11** and **ARID1A → KRTAP4-11** — keratin-associated-protein partners, biologically implausible; both **flip sign**. Exactly the discovery-only artifacts an orthogonal test should reject. E5.
- **SMARCA4 → CIT** flips sign (E5).
- **SMARCA4 → SMARCA2** narrowly misses (FDR = 0.053) and **ARID1A → EZH2** (FDR = 0.083) are *concordant but not significant* — real paralog/polycomb biology attenuated by cross-platform variance, retained as concordant-not-sig rather than over-claimed.

A metric that let implausible axes "replicate" would be worthless; this one does not.

### 1.6 Emergent hypothesis holds on the orthogonal platform

All three **chromatin-loss → WRN** axes replicate cross-platform: ARID1A → WRN (FDR = 7.3×10⁻⁸), KMT2D → WRN (FDR = 1.7×10⁻⁶), EP300 → WRN (FDR = 3.6×10⁻³). This strengthens the emergent "chromatin-remodeler loss confers WRN dependence" hypothesis — **with the caveat that the WRN precedent is MSI-driven** (Behan/Picco/van Wietmarschen), so the *specific* driver→WRN link is a pipeline nomination that is MSI-confounded, not a driver-specific clinical claim. Graded E3 on the strength of the WRN-in-vivo/MSI literature, not on a driver-specific trial.

---

## 2. Clinical-evidence ladder

### 2.1 Rubric (ordered, conservative)

| Grade | Meaning | Evidence required |
|---|---|---|
| **E1** Approved-in-context | approved drug for the node in the matched genomic context | FDA/EMA label + context match |
| **E2** Clinical-trial precedent | node in an active trial (ideally matched context) | NCT ID / trial citation |
| **E3** Preclinical in vivo | SL shown in animal / xenograft / PDX | primary-literature PMID/DOI |
| **E4** Preclinical in vitro / paralog-established | SL shown in cells, or well-established paralog buffering | primary-literature PMID/DOI |
| **E5** Computational-only | our pipeline + reproducibility only; no external biology yet | none (explicit nomination) |

Grades are assigned **conservatively**: E1/E2 require the drug to be MOA-consistent for the node **and** the context to match (or be explicitly, plausibly generalisable). **Promiscuous repurposing-hub hits cannot raise a grade.** Evidence grade (external precedent) and replication status (internal, orthogonal) are kept as **two orthogonal dimensions** — an axis can replicate strongly yet be E5, or be E1 yet untestable in Sanger.

### 2.2 Distribution

| Grade | n | | Replication status | n |
|---|---|---|---|---|
| E1 Approved-in-context | 4 | | replicated cross-platform | 41 |
| E2 Trial precedent | 8 | | concordant, not significant | 6 |
| E3 Preclinical in vivo | 6 | | discordant (sign flip) | 4 |
| E4 Preclinical in vitro / paralog | 10 | | insufficient n | 3 |
| E5 Computational-only | 30 | | not testable (CN axis) | 4 |
| **Total** | **58** | | **Total** | **58** |

`precedent_class`: none_found = 23, mechanistic_no_drug = 13, clinical_target_in_trials = 12, preclinical_in_vivo = 5, preclinical_in_vitro = 3, established_SL_drug_dev = 2. See `fig_evidence_ladder.png`.

### 2.3 Two-dimensional trust tier

Re-mapping the prior 1-D tiers onto (evidence grade × replication):

| Trust tier v2 | n | Meaning |
|---|---|---|
| **T1 anchored-and-replicated** | 5 | E1/E2 **and** replicated cross-platform — strongest honest position |
| **T2 anchored-precedent** | 21 | strong external precedent (E1/E2, or E3/E4 + replicated) |
| **T3 replicated, no precedent** | 22 | replicates on the orthogonal platform but no external biology yet |
| **T3 precedent-only** | 2 | precedent exists but did not (or could not) replicate here |
| **T4 discovery-only** | 8 | neither anchored nor replicated — explicit nominations |

See `fig_trust_tier_landscape.png` and `fig_clinical_actionability_map.png`.

---

## 3. Node narratives (re-cast with grades and caveats)

**BRAF → MAPK1 / MAP2K1 (E1, replicated).** Approved BRAF/MEK combinations; tissue-agnostic precedent [50, 48]. Both replicate cross-platform (FDR ≈ 2.6×10⁻⁵). The clearest anchored-and-replicated pair.

**BRCA1 → PARP1 (E1, replicated).** Approved PARP inhibitors in *BRCA*-mutant disease [67]. **Modality-limited caveat:** CRISPR knockout cannot model PARP *trapping*, the dominant clinical mechanism — so the *genetic* dependency replicates, but the CRISPR assay under-represents the pharmacology. Stated in the row citation.

**PTEN → PIK3CB (E2, replicated).** PI3Kβ dependence in PTEN-null cells [56]; GSK2636771 in PTEN-deficient trials [55]. Replicates (d = −0.75).

**STK11 → HDAC4 (E2, replicated).** HDAC-class agents in LKB1/STK11-mutant lung cancer and TNG260 (CoREST/HDAC1) in trials [60, 61]. **Caveat:** the *specific* partner is HDAC4; the precedent is **HDAC-class, not HDAC4-specific.** Graded E2 on class precedent, flagged in the row.

**VHL → EPAS1 (E1, insufficient_n).** Belzutifan (HIF-2α inhibitor) is FDA-approved for VHL disease [53] — a textbook matched-context E1. But only 3 VHL-mutant lines exist in the Sanger panel, so replication is **N/A (insufficient_n)**, not a failure. A case where external evidence is strong and the internal orthogonal test is simply under-powered.

**STAG2 → STAG1 (E4→ paralog-established, replicated).** Canonical paralog buffering [19, 20]; replicates (d = −1.93, FDR = 0.045).

**Chromatin → WRN (ARID1A / KMT2D / EP300 → WRN) (E3, all replicated).** WRN helicase dependence is an active clinical target in MSI/dMMR cancers [33, 43, 42]. All three replicate cross-platform. **Caveat: the precedent is MSI-driven, not driver-specific** — the driver→WRN edge is an MSI-confounded pipeline nomination; the E3 grade rests on WRN-in-vivo/MSI evidence.

**KRAS pathway (→ SHOC2 / DOCK5 / RAF1 / RAB10).** Pathway-level support is strong [82, 46, 47]; the *specific* partner edges are largely nominations (mostly E4). RAF1 has in-vivo regression data [46]. KRAS → STK33 is the cautionary counter-example (§1.5).

**MTAP → PRMT5 / MAT2A and CCNE1 → PKMYT1 / WEE1 (E2–E3, not testable here).** Strong trial-stage precedent (MRTX1719, AG-270, lunresertib, adavosertib) [36, 34, 23, 73], but these are **copy-number-driven** axes not re-derivable in the CERES matrix → `not_testable_cn`. Evidence grade high, replication N/A by assay design.

---

## 4. Production hardening (software sense)

- **Unit tests:** `tests/test_sl_clinical_validation.py`, **25/25 passing**, covering the replication statistic, status taxonomy (tested / insufficient_n / partner_absent / na_cn), evidence-grade mapping, 2-D tier logic, phase rank, and the ACH→ModelID bridge — each with a positive and a boundary case.
- **Single re-runnable entrypoint:** `run_clinical_validation.py` — download-if-missing (MD5-gated) → cross-platform replication → evidence/tier merge → matrix + receipts. Idempotent, parameterised paths, logged provenance. **Verified to reproduce the canonical numbers exactly** (ρ = 0.7908, r = 0.8125, 41 replicated) from a clean kernel.
- **Pinned data versions:** `data_versions.json` records source, file, and MD5 for every input (DepMap 24Q4, Sanger Figshare 9116732 [MD5 9a9e92d9…], mutations, model table).
- **Versioned schema:** `SCHEMA_VERSION = "v2.0-clinical-grade"`; the master matrix carries all prior columns plus the additive v2 columns (replication, evidence grade, precedent class, tier v2), preserving the 58-survivor / 44-reclassified split.

"Production grade" here means **reproducible software** (tests, pinned versions, deterministic entrypoint, versioned schema) — **not** a regulatory or GxP claim.

---

## 5. Limitations that remain (nothing swept under the rug)

1. **Not out-of-cohort.** Cell panels overlap 315/318; replication is cross-platform, not cross-cohort. Held-out rate = N/A.
2. **Not clinical proof.** No prospective human outcome data for any specific axis. Evidence grades are precedent/actionability, never efficacy.
3. **CN-driven axes untested here** (4 axes) — high external evidence, but not re-derivable in the CERES matrix.
4. **Under-powered controls** — VHL → EPAS1 and other axes with < 5 Sanger mutant lines cannot be judged internally.
5. **Modality mismatch for some E1s** — CRISPR KO does not model PARP trapping (BRCA1 → PARP1) or drug-specific pharmacology generally; genetic dependency ≠ drug response (Vermeulen 2025 [92]: pharmacological inhibitors often fail to phenocopy genetic knockouts).
6. **Precedent-vs-node gaps** — STK11 → HDAC4 (HDAC-class not HDAC4-specific), chromatin → WRN (MSI-driven not driver-specific), CDKN2A → FOSL1 (KRAS-context not CDKN2A-specific) are graded on the nearest defensible precedent, with the gap flagged in each row's citation.
7. **E5 axes are nominations** (30/58) — computational + reproducibility only; explicitly not differentiated biology yet.

---

## 6. Deliverables

| File | Contents |
|---|---|
| `clinical_evidence_ladder.csv` | per-axis evidence grade + replication status + citations (headline table) |
| `benchmark_independent_sanger.csv` + `_summary.json` | per-axis cross-platform replication, concordance stats, cell-line accounting |
| `sl_master_matrix.csv` / `.json` | extended master matrix (v2.0 schema, 58 survivors + 44 reclassified) |
| `data_versions.json` | pinned inputs + MD5s |
| `fig_independent_replication.png/.svg` | discovery vs Sanger effect-size scatter (ρ = 0.791) |
| `fig_evidence_ladder.png/.svg` | 58 axes across E1–E5 × replication status |
| `fig_replication_overlap.png/.svg` | cell-line accounting + per-axis testability |
| `fig_trust_tier_landscape.png/.svg` | 2-D trust landscape (grade × replication) |
| `fig_clinical_actionability_map.png/.svg` | druggability × evidence grade × replication |
| `run_clinical_validation.py` | single re-runnable entrypoint (in engine repo) |
| `tests/test_sl_clinical_validation.py` | 25 unit tests (in engine repo) |
| `SL_VALIDATION_RECEIPTS.json` | extended receipts (datasets, MD5s, test summary, schema, commit) |

*Citations `[N]` refer to the literature records compiled during this analysis; DOIs are given inline in each node's row of `clinical_evidence_ladder.csv`.*
