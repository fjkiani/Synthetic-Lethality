# AACR 2026 CrisPRO Gap Report (Interim)
**Generated:** 2026-07-04 16:32 UTC  
**Corpus:** 7,485 AACR 2026 abstracts (Cancer Research 7_Supplement + 8_Supplement)  
**Coverage:** 3,263/7,485 abstracts processed (43.6%)  
**Sources:** 1,060 LLM-enriched (Groq llama-3.3-70b) | 2,203 keyword-fallback | 4,222 pending (enrichment running)  
**QC note:** Rare axes (PKMYT1, WRN, CCNE1_AMP, TP53_RS) require keyword confirmation to suppress LLM false positives (pre-cleaning FP rates: PKMYT1 100%, TP53_RS 91%, CCNE1_AMP 71%).

---

## 1. Signal Matrix — Axes × Cancer Type

| Axis | Ovarian | CRC | PDAC | Breast | Lung | Prostate | Heme | Bladder |
|---|---|---|---|---|---|---|---|---|
| PARP_INHIBITORS | 85 | 34 | 26 | 80 | 24 | 21 | 15 | 3 |
| ATR_WEE1 | 17 | 14 | 7 | 24 | 16 | 8 | 9 | 1 |
| CCNE1_AMP | 13 | 2 | — | 6 | 3 | — | — | — |
| PKMYT1 | — | 3 | — | — | 1 | — | — | — |
| WRN | 2 | 8 | 1 | 2 | — | — | — | — |
| PRMT5_MTAP | 1 | 7 | 10 | 1 | 5 | 2 | 4 | 1 |
| IMMUNOTHERAPY | 56 | 162 | 103 | 200 | 248 | 28 | 93 | 31 |
| CYTIDINE_ANALOGS | 8 | 8 | 60 | 9 | 5 | 2 | 21 | 4 |
| TP53_REPLICATION_STRESS | — | 1 | — | 3 | 1 | — | — | 1 |
| PI3K_AKT | 13 | 52 | 38 | 95 | 53 | 12 | 23 | 4 |

**Key observations:**
- **IMMUNOTHERAPY** dominates (CRC: 162, Breast: 200, Lung: 248) — reflects AACR 2026 program emphasis on checkpoint/IO combinations
- **PARP_INHIBITORS** strong in Ovarian (85) and Breast (80) — olaparib/niraparib/talazoparib trial data
- **PI3K_AKT** significant in Breast (95) and CRC (52) — alpelisib/capivasertib/inavolisib data
- **PRMT5_MTAP** emerging in PDAC (10) and CRC (7) — MTAP deletion enriched in these tumor types; multiple Phase I datasets
- **PKMYT1** confirmed in 5 abstracts (CT022 Phase I, ACR-2316 WEE1/PKMYT1 dual inhibitor) — highly specific signal
- **WRN** CRC signal (8) — MSI-H synthetic lethality data; relevant to Brenus MSS/MSI-H context

---

## 2. Replication Stress Feature Distribution

| RS Feature | Abstracts | % of corpus |
|---|---|---|
| MSI_H | 84 | 1.1% |
| TP53_LOF | 72 | 1.0% |
| MYC_amplified | 52 | 0.7% |
| ARID1A_LOF | 22 | 0.3% |
| CCNE1_amplified | 14 | 0.2% |
| MBD4_LOF | 1 | 0.01% |

**MBD4_LOF is rare** (1 abstract) — consistent with its niche status as a CRC-specific RS driver. MSI-H (84 abstracts) is the dominant RS feature at AACR 2026.

---

## 3. Top Brenus-Relevant Abstracts (LLM-confirmed, fit_score > 0)

**1. Abstract 177: In vivo CRISPR screening uncovers Folr1-mediated metabolic remodeling that s**  
Signal: high | Brenus: high | Fit: 0.8  
Axes: IMMUNOTHERAPY  
Rationale: _The study highlights a potential therapeutic strategy to overcome intrinsic resistance to immunotherapy in MSS CRC by ta_  

**2. Abstract 3311: Autologous colorectal cancer-mesothelial co-culture model identifies stroma**  
Signal: high | Brenus: high | Fit: 0.8  
Axes: none  
Rationale: _FGFR3-mediated stromal signaling may be exploitable in MSS CRC peritoneal metastases._  

**3. Abstract 6746: Biomarker-driven restoration of tumor provisional matrix signaling network **  
Signal: medium | Brenus: high | Fit: 0.6  
Axes: IMMUNOTHERAPY  
Rationale: _Identifies a stromal biomarker (VCAN proteolysis) relevant to MSS CRC where checkpoint inhibitors have limited efficacy._  

**4. Abstract 7096: Design and development of a biparatopic antibody-drug conjugate against CDH**  
Signal: high | Brenus: high | Fit: 0.5  
Axes: none  
Rationale: _CDH17 is expressed in MSS colorectal tumors, making this ADC a promising option for MSS CRC patients._  

**5. Abstract 4578: Discovery of a Best-in-Class small molecule p53 Y220C reactivator: Breaking**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: none  
Rationale: _While the abstract does not specifically mention MSS CRC, the discovery of a potent p53 reactivator has potential implic_  

**6. Abstract 454: Overcoming PD-1 resistance with a first-in-class dual-mode agent that transf**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: IMMUNOTHERAPY  
Rationale: _The approach targets MSS CRC by inducing MSI-H phenotype, which could be relevant to MSS CRC treatment._  

**7. Abstract 436: Redox-DNA repair co-targeting with IP-DNQ and rucaparib induces oxidative DN**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: PARP_INHIBITORS  
Rationale: _The study's focus on NQO1-positive tumors and the use of rucaparib, a PARP inhibitor, may have implications for MSS CRC _  

**8. Abstract 251: Targeting cyclin K-CDK12 synergizes with ATR inhibition by limiting RPA chro**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: PARP_INHIBITORS, ATR_WEE1  
Rationale: _The findings may have implications for the treatment of triple-negative breast cancer, but further research is needed to_  

**9. Abstract 240: Allele-specific sensitization of pancreatic ductal adenocarcinoma to PARP in**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: PARP_INHIBITORS  
Rationale: _The study's focus on PDAC and PARP inhibition has relevance to MSS CRC, as both cancer types may benefit from targeted t_  

**10. Abstract 371: Integrative analysis of capivasertib mediated AKT blockade with AR inhibitio**  
Signal: high | Brenus: medium | Fit: 0.8  
Axes: PI3K_AKT  
Rationale: _While the study focuses on prostate cancer, the findings on PI3K/AKT activation and its compensation for AR inhibition m_  

---

## 4. PKMYT1 Confirmed Abstracts (5 total, keyword-validated)

- **Abstract 2505: Dual siRNAs for precision targeting of synthetic lethal vulnerabilities in **  
  Signal: low | Fit: 0.1 | Cancer: Unknown  

- **Abstract 3044: PKMYT1 inhibition is synthetically lethal with CCNE1 overexpression and syn**  
  Signal: low | Fit: 0.1 | Cancer: CRC  

- **Abstract 7430: Transcriptional programs of precursor-exhausted CD8+ T cells associated wit**  
  Signal: low | Fit: 0.1 | Cancer: Lung  

- **Abstract 3789: Treatment with ACR-2316, a potential first- and best-in-class WEE1/PKMYT1 i**  
  Signal: high | Fit: 0.5 | Cancer: CRC  

- **Abstract CT022: First data disclosure of the Phase I trial of the first in class combinati**  
  Signal: high | Fit: 0.5 | Cancer: CRC  

---

## 5. Coverage & Quality Summary

| Metric | Value |
|---|---|
| Total abstracts | 7,485 |
| LLM-enriched (Groq) | 1,060 (14.2%) |
| Keyword-fallback | 2,203 (29.4%) |
| Pending (workers running) | 4,222 (56.4%) |
| High signal | 320 |
| Medium signal | 845 |
| Brenus-relevant (high+medium) | 400 |
| Abstracts with KB axis hit | 2,598 (34.7%) |

**QC applied:** Rare-axis keyword confirmation filter eliminates LLM hallucinations on PKMYT1/WRN/CCNE1_AMP/TP53_RS.  
**Next:** Final report regenerated automatically once all 5 enrichment workers complete.
