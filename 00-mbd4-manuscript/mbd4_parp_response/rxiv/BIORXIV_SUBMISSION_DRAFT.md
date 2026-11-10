# bioRxiv gateway field pack — MBD4 PARP response preprint

**Status:** `SUBMISSION_READY = YES` (all required gateway fields filled; **do not submit** until Alpha clicks Submit)  
**Gateway URL:** https://submit.biorxiv.org/submit_gateway.dtl?newms=true  
**Bundle root:** `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/`  
**Source manuscript:** `rxiv/manuscript.md`  
**Upload PDF (exact path):** `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/mbd4_parp_response_biorxiv_submission.pdf` (~3.1 MB; verified 2026-11-10)  
**Last updated:** 2026-11-10 (complete field pack per Alpha; no LLC/RUO author metadata)

**Research framing:** Computational pharmacogenomics + machine learning on public, de-identified cancer cell lines. No new patient-identifiable data.

---

## Quick reference — paste order

| # | Field | Section |
|---|--------|---------|
| 1 | Article type | §1 |
| 2 | Title | §2 |
| 3 | Running title | §3 |
| 4 | Abstract | §4 |
| 5 | Authors + affiliations + corresponding | §5 |
| 6 | Subject category | §6 |
| 7 | Keywords | §7 |
| 8 | Competing interests | §8 |
| 9 | Funding | §9 |
| 10 | Data / code availability | §10 |
| 11 | Cover letter | §11 |
| 12 | Files to upload | §12 |
| 13 | License | §13 |
| 14 | Dual submission / prior preprint / related MS | §14 |
| 15 | Ethics (human / animal / clinical trial) | §15 |
| 16 | Author contributions (CRediT) | §16 |
| 17 | Acknowledgments | §17 |
| 18 | Co-author confirmation emails | §18 |
| 19 | ORCID | §19 |
| 20 | Pre-submit verification | §20 |

---

## 1. Article type / New manuscript

**Selection:** **New Results** (new manuscript; not a review, not a commentary).

**Gateway action:** Start at `submit_gateway.dtl?newms=true` → choose **New Results** (or equivalent “New manuscript / New Results” option).

---

## 2. Title (exact — paste verbatim)

```
A Three-Axis Framework for Therapy Prioritization in MBD4-Deficient Tumors: De-confounded Synthetic Lethality Modeling, Replication-Stress Pharmacology, and PARP1 Trapping Decoupling
```

---

## 3. Running title (short title)

**Character count:** 48 (within typical bioRxiv short-title limits; no truncation needed).

**Paste:**

```
MBD4-LOF: ATR vulnerability without PARPi sensitivity
```

---

## 4. Abstract (full text — paste verbatim)

```
Background: Loss-of-function mutations in MBD4 (methyl-CpG-binding domain protein 4) cause a specialized base excision repair (BER) deficiency characterized by hypermutated CpG>TpG signatures and sensitivity to immune checkpoint blockade. However, post-immunotherapy treatment prioritization remains empirically unmapped. In clinical oncology, DNA damage response (DDR) defects are frequently extrapolated to poly(ADP-ribose) polymerase inhibitor (PARPi) sensitivity by analogy to homologous recombination deficiency (HRD), without testing whether upstream glycosylase loss generates the nicked single-strand break (SSB) intermediates required for toxic PARP1-DNA trapping.

Methods: We developed a three-axis framework integrating mechanistic curation, large-scale pharmacogenomics (DepMap 24Q2 CRISPR/RNA-seq, n=1,517; GDSC2 dose-response profiling, n=928; PRISM Repurposing), and a 5-Layer De-confounded Synthetic Lethality Machine Learning (SL-ML) Engine (SL-XGBRanker). To prevent functional dilution from variants of uncertain significance, MBD4 loss-of-function (MBD4-LOF) was strictly restricted to truncating and canonical splice-site alleles (21 of 36 mutated lines), independently validated by nonsense-mediated mRNA decay (NMD; two-sided p2s = 0.0367), while excluding 15 lines harboring 19 non-truncating variants (18 missense, 1 inframe deletion). All primary pharmacogenomic comparisons were evaluated using two-sided Mann-Whitney U tests (p2s), bootstrap 95% confidence intervals (CIs), Benjamini-Hochberg false discovery rate (BH-FDR) multiplicity correction across six pre-specified candidate axes, and Two One-Sided Tests (TOST) for expression equivalence.

Results: Our analysis delineates three distinct therapeutic axes in MBD4-LOF tumors: (1) Cytidine Analog Vulnerability (Axis 1): Mechanistic structural evidence demonstrates that MBD4 loss impairs excision of 5-substituted pyrimidine analogs and deaminated bases, creating a direct biochemical bottleneck to gemcitabine and 5-fluorouracil incorporation. (2) ATR-Check1 Replication-Stress Vulnerability (Axis 2): In GDSC2 profiling, MBD4-LOF cell lines exhibit significant sensitivity to the ATR inhibitor ceralasertib (AZD6738). In the clinically relevant TP53-mutant background—which disables the G1/S checkpoint and forces reliance on ATR-mediated S/G2 arrest—MBD4-LOF lines (n=11 vs n=619 WT) demonstrate pronounced ceralasertib hypersensitivity (ΔLN_IC50 = -1.07, 95% CI [-1.62, -0.41], Cohen's d = -0.74, p2s = 0.0060, BH-FDR q = 0.036; ΔAUC = -0.087, d = -0.89, p2s = 0.0017), with 11/11 leave-one-out (LOO) jackknife iterations retaining two-sided significance (p2s ≤ 0.014). This genotype-associated dependency remains robust in the unstratified pan-cancer cohort (n=14 vs n=914 WT; ΔLN_IC50 = -0.73, 95% CI [-1.31, -0.11], d = -0.50, p2s = 0.0430; ΔAUC = -0.056, d = -0.56, p2s = 0.0259) and after purging microsatellite instability-high (MSI-H) lines (n=10 vs n=906 WT MSS; ΔLN_IC50 = -0.91, 95% CI [-1.50, -0.22], d = -0.62, p2s = 0.0307). By contrast, Wee1 inhibition (adavosertib; n=15 vs n=920 WT; ΔLN_IC50 = -0.35, p2s = 0.1667) does not reach statistical significance. (3) Absence of Single-Agent PARPi Sensitivity (Axis 3): While PARP1 transcript abundance strongly governs pan-cancer PARPi cytotoxicity across n=481 cell lines (ρ = -0.416, p < 10^-21 across five clinical PARPi), MBD4-LOF tumors (n=19 vs n=1,498 WT) show zero compensatory PARP1 upregulation (median 6.77 vs 6.66 log2[TPM+1], Δmedian = +0.11, 95% CI [-0.39, +0.71], p2s = 0.605; TOST equivalence p_TOST = 0.0199 within ± 0.50 log2[TPM+1]). Direct GDSC2 evaluation of MBD4-LOF (n=8) versus WT (n=478) confirms complete absence of single-agent PARPi sensitization across talazoparib (p2s = 0.973), olaparib (p2s = 0.522), niraparib (p2s = 0.877), rucaparib (p2s = 0.804), veliparib (p2s = 0.538), and the composite 5-drug mean PARPi Z-score (Δmedian Z = -0.024, 95% CI [-0.84, +0.60], p2s = 0.841). Under Leave-One-Anchor-Out cross-validation across 12 canonical DDR genes (370 gene-drug pairs), our 5-Layer De-confounded SL-ML Engine achieves an out-of-fold AUROC of 0.9716, AUPRC of 0.8007, and Top-1 accuracy of 82.86%, zero-shot ranking gemcitabine (#1) and ceralasertib (#2) as the top therapeutic dependencies for held-out MBD4 while rejecting all PARP inhibitors (ranks #13, #27, #28, #31).

Conclusion: MBD4-LOF tumors exhibit a selective, TP53-potentiated pharmacological dependency on ATR kinase inhibition alongside mechanistic vulnerability to cytidine analogs, but lack both the nicked BER repair intermediates and the PARP1 target expression required for single-agent PARP inhibitor efficacy.
```

*Note: LaTeX math in the manuscript PDF uses proper notation; gateway abstract uses plain-text equivalents above for paste compatibility. If the gateway accepts rich text, you may paste from `manuscript.md` Abstract section (lines 57–63) instead for identical math formatting.*

---

## 5. Author list (submission order)

| Order | Display name | Email | Affiliation (paste) | Corresponding? | ORCID (§19) |
|-------|--------------|-------|---------------------|----------------|-------------|
| 1 | Fahad Kiani | fahad@crispro.io | CrisPRO.org, New York, NY | **Yes** | Not provided — leave blank / N/A |
| 2 | Ayesha Kiani | ayesha.kiani@gmail.com | New York University, New York, NY | No | Not provided — leave blank / N/A |
| 3 | Rahima Nayeem, MD | nayeem.rahima@gmail.com | CrisPRO.org, New York, NY | No | Not provided — leave blank / N/A |

**Corresponding author (gateway):** Fahad Kiani — fahad@crispro.io

**Affiliation numbering (if gateway uses numbered affiliations):**

- **Affiliation 1:** CrisPRO.org, New York, NY → Fahad Kiani; Rahima Nayeem, MD  
- **Affiliation 2:** New York University, New York, NY → Ayesha Kiani  

**Author name format (if gateway wants surname first):** Kiani, Fahad; Kiani, Ayesha; Nayeem, Rahima (MD)

---

## 6. Subject category / taxonomy

**Primary subject (single choice — paste/select):** **Cancer Biology**

*(Optional secondary, if the gateway allows a second category: **Bioinformatics** — not required.)*

---

## 7. Keywords (comma-separated; 9 terms)

```
MBD4, Base Excision Repair, PARP Inhibitor, ATR Inhibitor, Ceralasertib, Gemcitabine, Synthetic Lethality, De-confounded Machine Learning, Replication Stress
```

---

## 8. Competing interest statement (paste)

```
Fahad Kiani is the founder of CrisPRO.org (https://crispro.org). Ayesha Kiani declares no competing interests. Rahima Nayeem, MD declares no competing interests. All analyses were conducted using publicly available, de-identified cell-line pharmacogenomic datasets.
```

---

## 9. Funding statement (paste)

```
This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.
```

---

## 10. Data and code availability (paste)

```
All pharmacogenomic datasets analyzed in this study are publicly available from the Broad Institute DepMap Portal (DepMap Public 24Q2: https://depmap.org/portal/) and the Sanger GDSC2 repository (https://www.cancerrxgene.org/). Analysis scripts, frozen JSON statistical receipts (including canonical_atr_wee1_rerun.json, gdsc2_parpi_mbd4_summary.json, and s2_evaluation_metrics.json), and figure-generation code are available at https://github.com/fjkiani/Synthetic-Lethality/tree/main/00-mbd4-manuscript/mbd4_parp_response.
```

---

## 11. Cover letter (final — paste if requested)

```
Dear Editors,

We submit for posting as a preprint on bioRxiv:

Title: A Three-Axis Framework for Therapy Prioritization in MBD4-Deficient Tumors: De-confounded Synthetic Lethality Modeling, Replication-Stress Pharmacology, and PARP1 Trapping Decoupling

This work integrates mechanistic BER literature, pan-cancer pharmacogenomic screens (DepMap 24Q2, GDSC2), and a de-confounded synthetic-lethality ranking engine to distinguish three therapeutic axes in strict MBD4 loss-of-function models: cytidine-analog vulnerability, TP53-potentiated ATR inhibitor sensitivity (ceralasertib), and absence of single-agent PARP inhibitor sensitivity despite intact PARP1 expression (TOST equivalence). Findings are derived from de-identified cancer cell lines and are intended as hypothesis-generating and translational prioritization support, not as clinical recommendations.

We believe this manuscript will interest the cancer biology, DDR, and pharmacogenomics communities, particularly groups studying BER versus HRD and PARP inhibitor biomarker assumptions.

Corresponding author: Fahad Kiani (fahad@crispro.io).

Sincerely,
Fahad Kiani (on behalf of all authors)
```

---

## 12. File checklist (absolute paths)

| Item | Absolute path | Notes |
|------|----------------|-------|
| **Main manuscript PDF (REQUIRED)** | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/mbd4_parp_response_biorxiv_submission.pdf` | 25 pages; references rendered; no raw `[@` citations in PDF (PyPDF2 check 2026-11-10) |
| **Source markdown** | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/manuscript.md` | Not uploaded unless gateway asks for source |
| **Bibliography** | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/references.bib` | Baked into PDF via citeproc |
| **Build instructions** | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/BUILD.md` | Rebuild only if manuscript/bib changes |

**Figures (embedded in PDF; separate PNGs on disk for reference only):**

| Figure | Absolute path |
|--------|----------------|
| Fig 1 | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES/fig1_evidence_matrix.png` |
| Fig 2 | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES/fig2_ceralasertib_volcano.png` |
| Fig 3 | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES/fig3_stress_tests.png` |
| Fig 4 | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES/fig4_parp1_expression.png` |
| Fig 5 | `/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES/fig5_mechanism_model.png` |

**Supplementary files:** **None required** for bioRxiv (receipt JSONs under `artifacts/` optional if Alpha wants a supplemental zip; not in main MS).

**Rebuild PDF (only if `manuscript.md` or `references.bib` change):**

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/00-mbd4-manuscript/mbd4_parp_response/rxiv
pandoc manuscript.md \
  -o mbd4_parp_response_biorxiv_submission.pdf \
  --pdf-engine=tectonic \
  --filter pandoc-crossref \
  --citeproc
```

---

## 13. License / copyright

**Recommendation (filled for Alpha to confirm at submit):** **CC-BY 4.0** (Creative Commons Attribution 4.0 International)

**Rationale:** bioRxiv default/recommended license for broad reuse and citation; allows commercial reuse with attribution. **Alternative if Alpha prefers restricted reuse:** CC-BY-NC 4.0 — use only if Alpha explicitly chooses NC at the license step.

**Copyright holder:** Authors (standard author-retained copyright on bioRxiv preprint).

---

## 14. Dual submission / prior preprint / related manuscripts

| Question | Answer |
|----------|--------|
| **Previously posted as preprint?** | **No** — first submission to bioRxiv (no medRxiv/SSRN/arXiv DOI for this manuscript in repo) |
| **Concurrent journal submission?** | **No** — not concurrently under review at a journal (per current package; update only if status changes) |
| **Dual submission to another preprint server at same time?** | **No** |
| **Related manuscripts / companion papers?** | **No** separate manuscript; GitHub repo is the code/data companion only (`https://github.com/fjkiani/Synthetic-Lethality/tree/main/00-mbd4-manuscript/mbd4_parp_response`) |
| **Cross-posting / transfer from another server?** | **No** |

---

## 15. Ethics and regulatory declarations

| Gateway field | Answer |
|---------------|--------|
| **Clinical trial registration** | **No** / Not applicable |
| **Human subjects research** | **No** — analysis of publicly available, de-identified cancer cell-line databases (DepMap, GDSC2) only; no new human subjects data |
| **Animal research** | **No** |
| **IRB / ethics approval** | **Not applicable** (no human subjects data collection) |

---

## 16. Author contributions (CRediT — paste if gateway asks)

```
Fahad Kiani: Conceptualization, Methodology, Software, Formal analysis, Visualization, Writing – original draft, Writing – review & editing, Supervision, Project administration, Corresponding author.
Ayesha Kiani: Writing – review & editing.
Rahima Nayeem: Writing – review & editing, Clinical contextualization.
```

---

## 17. Acknowledgments

**Paste (if gateway has a separate acknowledgments box):**

```
No funding acknowledgments. We thank the Broad Institute DepMap and Sanger GDSC teams for public data resources.
```

---

## 18. Co-author confirmation emails (bioRxiv author approval workflow)

Co-authors who must confirm authorship (all except corresponding, unless gateway emails everyone):

| Co-author | Email for confirmation |
|-----------|-------------------------|
| Ayesha Kiani | ayesha.kiani@gmail.com |
| Rahima Nayeem, MD | nayeem.rahima@gmail.com |

**Corresponding author account email:** fahad@crispro.io (submitter / corresponding).

---

## 19. ORCID

No ORCID iDs are present in `manuscript.md`, `references.bib`, or the `rxiv/` bundle.

| Author | ORCID |
|--------|-------|
| Fahad Kiani | **Not applicable** — leave gateway ORCID field blank or skip |
| Ayesha Kiani | **Not applicable** — leave blank |
| Rahima Nayeem, MD | **Not applicable** — leave blank |

*Alpha may optionally add ORCIDs at submit time without changing the manuscript.*

---

## 20. Pre-upload verification (completed 2026-11-10)

| Check | Result |
|-------|--------|
| Required-field TODOs in this pack | **0** |
| PDF exists at path in §12 | **Yes** (~3.1 MB; 3,256,733 bytes; rebuild 2026-11-10) |
| Broken `[@cite]` in PDF text layer | **No** |
| References section in PDF | **Yes** |
| All three authors on PDF title page (PyPDF2 p.1) | **Yes** — Fahad Kiani, Ayesha Kiani, Rahima Nayeem, MD |
| Author metadata locked (no LLC/RUO in author block) | **Yes** — CrisPRO.org, NYU, CrisPRO.org |
| Submit button clicked by agent | **No** (forbidden) |

---

*End of bioRxiv submission field pack. Alpha: open gateway → paste fields §1–§19 → upload PDF from §12 → confirm **CC-BY 4.0** at §13 → complete co-author emails §18 → Submit manually.*
