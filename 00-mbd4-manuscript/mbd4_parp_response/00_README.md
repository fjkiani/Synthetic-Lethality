# MBD4 Deficiency: Dual Therapeutic Vulnerability & Mechanistic Characterization

**MBD4-LOF Manuscript Project** | bioRxiv-ready preprint package | *Research Use Only (RUO)*

**Public GitHub repository:** [github.com/crispro-ai/MBD4-LOF-Dual-Therapeutic-Vulnerability](https://github.com/crispro-ai/MBD4-LOF-Dual-Therapeutic-Vulnerability)

---

## Overview

This repository contains the manuscript and supporting data artifacts for our study on the therapeutic landscape of **MBD4 loss-of-function (LOF)** tumors. We define MBD4-LOF as a base excision repair (BER) defective state with distinct synthetic lethalities at the replication fork.

### Key Findings
1.  **Axis A: Cytidine Analog Synthetic Lethality** – Confirmed as the "gold-standard" vulnerability (gemcitabine/cytarabine).
2.  **Axis B: Hypermutator Sentinel** – MBD4-LOF drives a CpG>TpG hypermutator phenotype predictive of immune checkpoint inhibitor (ICI) response.
3.  **Axis C: ATR Checkpoint Dependency** – Novel discovery of ATRi (ceralasertib) sensitivity in MBD4-LOF models, validated across four orthogonal confound stress tests.
4.  **Mechanistic Falsification** – We falsify the hypothesis that PARP1 transcriptional upregulation drives PARPi sensitivity in MBD4-LOF.

---

## Repository Structure

```tree
.
├── README.md               # GitHub landing (links into this tree)
├── 00_README.md            # Detailed project readme
├── SUPPLEMENTAL_DATA_PROVENANCE.md  # Supplemental tables / receipts audit
├── MISSING_RECEIPTS_PLAN.md # Strategic plan (if present)
├── artifacts/              # Data receipts & analysis outputs (JSON/CSV)
│   ├── preclinical/        # GDSC2/DepMap stratification results
│   └── clinical/           # cBioPortal/clinical outcomes queries
└── rxiv/                   # Manuscript source
    ├── manuscript.md       # Full text (Markdown/Pandoc)
    ├── references.bib      # Bibliography
    ├── FIGURES/            # Generated publication-quality figures
    └── SUPPLEMENTARY/      # Supplemental tables and data
```

---

## Quick Links

-   **[Manuscript Draft](rxiv/manuscript.md)**: The current version of the preprint.
-   **[Supplemental data provenance](SUPPLEMENTAL_DATA_PROVENANCE.md)**: Audit trail for supplemental tables and receipts.
-   **[Figures](rxiv/FIGURES/)**: Visual summaries of the findings.

---

## Data & Reproducibility

This project utilizes the **CrisPRO Multimodal Evidence Matrix** framework as the organizing environment for evidence integration. **Public reproducibility for this paper**—manuscript source, figure scripts, canonical reruns, and frozen artifacts—is provided only via the GitHub repository linked at the top of this file.

All claims categorized as "Verified" are backed by explicit receipts located in the `artifacts/` directory.

-   **DepMap 24Q2**: Somatic mutation and expression data.
-   **GDSC2**: Pharmacological screening results.
-   **cBioPortal**: Clinical outcome data.
-   **SL Agent**: Automated hypothesis fuser results.

---

## Current Status

-   **Axis A**: ✅ **Verified** (Gold standard calibrated).
-   **Axis B**: 🟡 **Partial** (Case-level receipts; cohort expansion required).
-   **Axis C (ATRi)**: 🟢 **Strong** (Multimodally validated; confound-purged).
-   **Axis C (PARPi)**: ❌ **Falsified** (Transcriptional model rejected).

For contributors: See **MISSING_RECEIPTS_PLAN.md** (if present) for current blocking items.


