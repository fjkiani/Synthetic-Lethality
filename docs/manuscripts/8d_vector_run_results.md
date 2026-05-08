# 8D Mechanism Fit Run — GBM Trial Classes vs NSCLC BrM Patient Profiles

**Generated:** 2026-05-08  
**Vector version:** 8D.v1 `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux, RSS]`  
**Formula:** `combined = (0.7 × eligibility) + (0.3 × fit)` where `fit = (patient·trial) / ‖trial‖₂`  
**Eligibility:** Uniform 0.75 for all trial classes — mechanism alignment is the variable.  
**RUO:** Research Use Only. Not cleared for clinical decision-making.

## What this run shows

The manuscript's 6 GBM trial failure classes are encoded as 8D mechanism vectors representing
what each drug was designed to hit. These are scored against 3 NSCLC BrM patient profiles.
High mechanism fit = the drug's mechanism aligns with the patient's dominant pathway state.
**High fit in a profile that failed in GBM = the SL engine identifies where the mechanism
was right but the context was wrong.** Low fit = the mechanism was mismatched to the biology.

---

## Profile: ZEB1-high NSCLC BrM (EGFR-wt, KRAS-mut)

**Canonical SL exploit profile. KRAS-mut → MAPK dominant. ZEB1-high EMT state → RSS burden. Moderate IO (PD-L1 via ZEB1→CD274). ITGAV dependency via ZEB1→ITGAV SL axis.**

Patient 8D vector: `{'ddr': 0.1, 'mapk': 0.7, 'pi3k': 0.4, 'vegf': 0.2, 'her2': 0.0, 'io': 0.4, 'efflux': 0.0, 'rss': 0.6}`

| Rank | Trial Class | Mech Fit | Combined | Dominant Axes |
|------|-------------|----------|----------|---------------|
| 1 | EGFR TKIs (erlotinib/gefitinib) | 0.6325 | 0.7148 | mapk=0.490, pi3k=0.200 |
| 2 | Checkpoint inhibitors (CheckMate 143/498/548) | 0.5692 | 0.6958 | io=0.360, rss=0.180 |
| 3 | Cilengitide (CENTRIC) | 0.5207 | 0.6812 | pi3k=0.240, vegf=0.100, mapk=0.070 |
| 4 | mTOR/PI3K inhibitors (everolimus/BKM120) | 0.4749 | 0.6675 | pi3k=0.360, mapk=0.070 |
| 5 | Rindopepimut (ACT IV) | 0.4000 | 0.6450 | io=0.280 |
| 6 | Bevacizumab (AVAglio/RTOG 0825) | 0.2000 | 0.5850 | vegf=0.180 |

**Interpretation — EGFR TKIs rank #1 — but the SL engine routes around them**

EGFR TKIs score highest (mfit=0.63) because KRAS-mut drives MAPK, which overlaps with the EGFR→MAPK axis. But EGFR TKIs fail in KRAS-mut NSCLC (KRAS is downstream of EGFR — the target is wrong). The SL engine's ZEB1→ITGAV route is not in this trial class list because it is a non-antigenic, non-pathway-inhibition mechanism — it operates in a resistance class the 8D ranker does not score against any of these 6 failed approaches. Cilengitide ranks #3 (mfit=0.52) — the PI3K/VEGF overlap is real, but the SL engine adds ZEB1 patient selection that CENTRIC lacked. Rindopepimut ranks #5 (mfit=0.40) — IO alignment is moderate because ZEB1-high tumors are PD-L1-high, but the antigen escape problem remains. The SL engine is designed to exploit the gap between ranks 3 and 1.

---

## Profile: EGFR-mutant NSCLC BrM

**SOC anchor profile. EGFR-mut → MAPK dominant, PI3K activated. IO-cold (low TMB, immune exclusion in EGFR-mut). Low RSS. Osimertinib-eligible.**

Patient 8D vector: `{'ddr': 0.1, 'mapk': 0.8, 'pi3k': 0.5, 'vegf': 0.1, 'her2': 0.2, 'io': 0.1, 'efflux': 0.3, 'rss': 0.2}`

| Rank | Trial Class | Mech Fit | Combined | Dominant Axes |
|------|-------------|----------|----------|---------------|
| 1 | EGFR TKIs (erlotinib/gefitinib) | 0.9625 | 0.8138 | mapk=0.560, pi3k=0.250, efflux=0.180 |
| 2 | mTOR/PI3K inhibitors (everolimus/BKM120) | 0.5853 | 0.7006 | pi3k=0.450, mapk=0.080 |
| 3 | Cilengitide (CENTRIC) | 0.5461 | 0.6888 | pi3k=0.300, mapk=0.080, vegf=0.050 |
| 4 | Checkpoint inhibitors (CheckMate 143/498/548) | 0.1581 | 0.5724 | io=0.090, rss=0.060 |
| 5 | Bevacizumab (AVAglio/RTOG 0825) | 0.1000 | 0.5550 | vegf=0.090 |
| 6 | Rindopepimut (ACT IV) | 0.1000 | 0.5550 | io=0.070 |

**Interpretation — EGFR TKIs dominate (mfit=0.96) — osimertinib solves the delivery problem**

EGFR-mut profile aligns almost perfectly with EGFR TKI mechanism (mfit=0.96). This is the SOC anchor — osimertinib is already approved for this profile and solves the P-gp exclusion problem that killed erlotinib/gefitinib in GBM. Checkpoint inhibitors rank #4 (mfit=0.16) — IO-cold profile, consistent with known poor checkpoint response in EGFR-mut NSCLC. Rindopepimut ranks last (mfit=0.10) — antigen-targeted vaccine has near-zero mechanism fit in an EGFR-mut, IO-cold profile. The 8D run confirms the SOC anchor assignment in brm.json is mechanistically correct.

---

## Profile: PTEN-null ZEB1-high NSCLC BrM

**Candidate convergence signal profile. PTEN-null → PI3K hyperactivated. ZEB1-high → RSS burden. IO-suppressed (PTEN-loss immune exclusion). Candidate enrichment for ITGAV dependency (hypothesis, not confirmed).**

Patient 8D vector: `{'ddr': 0.1, 'mapk': 0.5, 'pi3k': 0.9, 'vegf': 0.2, 'her2': 0.0, 'io': 0.2, 'efflux': 0.0, 'rss': 0.4}`

| Rank | Trial Class | Mech Fit | Combined | Dominant Axes |
|------|-------------|----------|----------|---------------|
| 1 | mTOR/PI3K inhibitors (everolimus/BKM120) | 0.9497 | 0.8099 | pi3k=0.810, mapk=0.050 |
| 2 | Cilengitide (CENTRIC) | 0.8763 | 0.7879 | pi3k=0.540, vegf=0.100, mapk=0.050 |
| 3 | EGFR TKIs (erlotinib/gefitinib) | 0.7334 | 0.7450 | pi3k=0.450, mapk=0.350 |
| 4 | Checkpoint inhibitors (CheckMate 143/498/548) | 0.3162 | 0.6199 | io=0.180, rss=0.120 |
| 5 | Bevacizumab (AVAglio/RTOG 0825) | 0.2000 | 0.5850 | vegf=0.180 |
| 6 | Rindopepimut (ACT IV) | 0.2000 | 0.5850 | io=0.140 |

**Interpretation — mTOR/PI3K ranks #1 (mfit=0.95) — but this is exactly the feedback loop trap**

PTEN-null profile aligns maximally with mTOR/PI3K inhibitors (mfit=0.95) — PI3K is hyperactivated, so the mechanism fits perfectly. But this is the AKT feedback loop trap: mTOR inhibition reactivates AKT in PTEN-null tumors. The 8D ranker scores mechanism alignment, not resistance architecture. Cilengitide ranks #2 (mfit=0.88) — PI3K/VEGF overlap is strong in PTEN-null context. This is the candidate convergence signal: PTEN-null tumors that resist mTOR inhibitors (rank 1) may be enriched for ITGAV dependency (cilengitide, rank 2). Rindopepimut ranks last (mfit=0.20) — IO-suppressed profile, antigen-targeted vaccine has minimal mechanism fit. The 8D run supports the manuscript's PTEN-as-candidate-convergence-signal argument without overclaiming.

---

## Cross-Profile Summary

| Trial Class | ZEB1-high KRAS-mut | EGFR-mut | PTEN-null ZEB1-high | Key Insight |
|-------------|-------------------|----------|---------------------|-------------|
| Cilengitide (CENTRIC) | #3 (0.52) | #3 (0.55) | #2 (0.88) | PI3K/VEGF overlap real; ZEB1 selection adds what CENTRIC lacked |
| Bevacizumab (AVAglio/RTOG 0825) | #6 (0.20) | #5 (0.10) | #5 (0.20) | VEGF-only; low fit across all profiles — invasion escape confirmed |
| EGFR TKIs (erlotinib/gefitinib) | #1 (0.63) | #1 (0.96) | #3 (0.73) | MAPK/PI3K dominant; osimertinib solves delivery; KRAS-mut mismatch |
| Rindopepimut (ACT IV) | #5 (0.40) | #6 (0.10) | #6 (0.20) | IO-only; lowest fit in IO-cold/suppressed profiles — antigen escape irrelevant to SL |
| Checkpoint inhibitors (CheckMate 143/498/548) | #2 (0.57) | #4 (0.16) | #4 (0.32) | IO+RSS; moderate in ZEB1-high (PD-L1 link); cold in EGFR-mut/PTEN-null |
| mTOR/PI3K inhibitors (everolimus/BKM120) | #4 (0.47) | #2 (0.59) | #1 (0.95) | PI3K dominant; maximal in PTEN-null but AKT feedback trap applies |

---

## Critical Limitations

**1. Mechanism alignment ≠ clinical efficacy.**
The ranker scores whether a drug's *intended mechanism* overlaps with a patient's *dominant pathway state*. It does not score whether the drug would work. High mechanism fit in a profile that failed in GBM means the mechanism was contextually relevant — not that the drug would succeed in NSCLC BrM.

**2. Resistance architecture is not encoded.**
The 8D vector captures pathway activation state, not resistance escape routes. mTOR/PI3K inhibitors rank #1 for PTEN-null because PI3K is hyperactivated — but the AKT feedback loop (mTORC1 inhibition → IRS-1 → PI3K reactivation) is a resistance mechanism that operates *within* the PI3K axis and is invisible to the ranker. High mechanism fit in this case identifies the correct pathway but not the correct drug.

**3. EGFR TKI ranking for KRAS-mut is a known false positive.**
EGFR TKIs rank #1 for ZEB1-high/KRAS-mut because KRAS-mut drives MAPK, which overlaps with the EGFR→MAPK axis. But KRAS is downstream of EGFR — EGFR TKIs do not suppress KRAS-driven MAPK signaling. This is a well-documented clinical failure (KRAS-mut NSCLC does not respond to EGFR TKIs). The ranker correctly identifies the pathway overlap but cannot detect the target-pathway mismatch.

**4. Trial class vectors are curated, not empirical.**
The 8D vectors for each GBM trial class encode the mechanism the drug was *designed* to hit, based on published pharmacology. They are not derived from patient response data, transcriptomic profiling, or CRISPR screens. They represent a mechanistic prior, not a measured signal.

**5. Patient profile vectors are biologically grounded but not patient-derived.**
The three patient profiles are constructed from known biology (KRAS-mut → MAPK, PTEN-null → PI3K, etc.) and RSS encoding from Konstantinopoulos 2021. They are not derived from individual patient data. They represent canonical molecular archetypes for NSCLC BrM subgroups.

**Correct use of this output:** Use the ranked tables to identify which GBM failure modes are mechanistically proximal to each NSCLC BrM profile — and therefore which manuscript arguments are most relevant to each patient subgroup. Do not use combined scores as clinical predictions.

---

*RUO — Research Use Only. Not cleared for clinical decision-making.*  
*Provenance: brm_targetability_matrix_v3.json (commit 64258b6) | cilengitide entry (commit f64758a) | DepMap 24Q4 | 95 NSCLC lines*