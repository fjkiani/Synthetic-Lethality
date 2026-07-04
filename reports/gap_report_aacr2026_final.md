# AACR 2026 Gap Report — CrisPRO Fit-Gap Intelligence
**Generated:** 2026-07-04 17:25 UTC  
**Corpus:** 7,485 AACR 2026 abstracts  
**Enrichment:** 1,527 LLM (20.4%) + 2,041 keyword fallback (27.3%) = 3,568 with axis data (47.7%)  
**QC:** Rare-axis keyword confirmation applied (PKMYT1, WRN, CCNE1_AMP, TP53_RS)

---

## 1. Axis Signal Counts (Post-QC)

| Axis | Total | CRC | Ovarian | Breast | PDAC | Lung |
|------|------:|----:|--------:|-------:|-----:|-----:|
| PARP_INHIBITORS | 366 | 31 | 82 | 77 | 23 | 19 |
| ATR_WEE1 | 137 | 13 | 16 | 21 | 8 | 14 |
| CCNE1_AMP | 31 | 2 | 13 | 6 | 0 | 3 |
| PKMYT1 | 5 | 3 | 0 | 0 | 0 | 1 |
| WRN | 28 | 8 | 2 | 2 | 1 | 0 |
| PRMT5_MTAP | 50 | 7 | 1 | 1 | 9 | 5 |
| IMMUNOTHERAPY | 1515 | 148 | 54 | 193 | 90 | 219 |
| CYTIDINE_ANALOGS | 122 | 8 | 8 | 9 | 53 | 4 |
| TP53_REPLICATION_STRESS | 8 | 1 | 0 | 3 | 0 | 1 |
| PI3K_AKT | 382 | 47 | 11 | 83 | 34 | 43 |

---

## 2. Replication Stress Feature Distribution

| RS Feature | Count |
|------------|------:|
| TP53_LOF | 138 |
| MSI_H | 124 |
| MYC_amplified | 103 |
| ARID1A_LOF | 38 |
| CCNE1_amplified | 14 |

---

## 3. CrisPRO White Space Analysis

### 3.1 MSS CRC Positioning
- **MSS CRC + vaccine abstracts:** 1 (STC-1010 only — Brenus uncontested)
- **MSS CRC + replication stress:** 0 (zero — white space confirmed)
- **STC-1010 mentions:** 1 (CT051 — Brenus's own abstract)

### 3.2 STC-1010 at AACR 2026 (CT051)
**Title:** Abstract CT051: From preclinical models to first-in-human evaluation of STC-1010 immunotherapy in unresectable advanced colorectal cancer  
**DOI:** https://doi.org/10.1158/1538-7445.am2026-ct051  
**Key finding:** First-in-human BreAK CRC-001 Phase I/IIa data presented. mSTC-1010 + FOLFOX reduced tumor volumes vs control in syngeneic CRC models; increased CD8⁺ T cell infiltration. MSS CRC (85-95% of CRC) remains the target population.

### 3.3 MBD4_LOF White Space
- **MBD4_LOF RS feature:** 0 abstracts (near-zero — white space confirmed)
- **ATR_WEE1 total:** 137 abstracts — active competitive space
- **MBD4 + ATR_WEE1 intersection:** Brenus's unique angle (MBD4_LOF → ATR sensitivity) not represented

---

## 4. Competitive Landscape by Axis

### PKMYT1 (n=5)
- **Abstract 2505: Dual siRNAs for precision targeting of synthetic lethal vulnerabilities in **  
  Cancer:  | Signal: medium | Fit: 0.00  
  Drugs: []  

- **Abstract 3044: PKMYT1 inhibition is synthetically lethal with CCNE1 overexpression and syn**  
  Cancer: CRC | Signal: medium | Fit: 0.00  
  Drugs: []  

- **Abstract 7430: Transcriptional programs of precursor-exhausted CD8+ T cells associated wit**  
  Cancer: Lung | Signal: medium | Fit: 0.00  
  Drugs: []  

- **Abstract 3789: Treatment with ACR-2316, a potential first- and best-in-class WEE1/PKMYT1 i**  
  Cancer: CRC | Signal: medium | Fit: 0.00  
  Drugs: []  

- **Abstract CT022: First data disclosure of the Phase I trial of the first in class combinati**  
  Cancer: CRC | Signal: medium | Fit: 0.00  
  Drugs: []  


### WRN (n=28)
- **Abstract 4436: Sutro’s site-specific dual-payload ADCs combining TOPO1i and DNA damage res**  
  Cancer: string | Signal: medium  

- **Abstract 423: ETX-880, a potential best-in-class, oral, highly potent and selective covale**  
  Cancer: MSI-H cancers | Signal: high  

- **Abstract 244: DHX9 inhibition is synthetically lethal with homologous recombination defici**  
  Cancer:  | Signal: medium  

- **Abstract 233: Preclinical development of EIK1005, a potent and selective inhibitor of Wern**  
  Cancer:  | Signal: medium  

- **Abstract 1685: The HER2-targeting dual-payload antibody-drug conjugate combining a topoiso**  
  Cancer: Breast | Signal: medium  


### PRMT5_MTAP (n=50)
- **Abstract 741: Illuminate your targets: Endogenous HiBiT knock-in monoclones—inventory + ra**  
  Cancer:  | Signal: none | Drugs: []  

- **Abstract 7022: Investigating PRMT5 as a therapeutic target in EGFR TKI resistant NSCLC and**  
  Cancer: Lung | Signal: low | Drugs: []  

- **Abstract 4504: Synergistic antitumor activity of the MTA-cooperative PRMT5 inhibitor ABSK1**  
  Cancer: None | Signal: high | Drugs: ['ABSK131', 'KRAS inhibitors', 'EGFR inhibitors', 'chemotherapy agents', 'MAT2A inhibitors']  

- **Abstract 7077: ISM1745, an MTA-cooperative PRMT5 inhibitor for the treatment of MTAP-delet**  
  Cancer: CRC | Signal: high | Drugs: ['ISM1745', 'ISM3412']  

- **Abstract 5941: Large-scale analysis reveals distinct molecular subtypes in real-world gast**  
  Cancer: Gastric | Signal: medium | Drugs: []  


---

## 5. Brenus-Relevant Abstracts (High Priority)

| # | Title | Cancer | Axes | Signal | Fit |
|---|-------|--------|------|--------|-----|
| 1 | Abstract 177: In vivo CRISPR screening uncovers Folr1-mediat | CRC | IMMUNOTHERAPY | high | 0.80 |
| 2 | Abstract 4351: A sex-specific role of CD36 targeting therapy | CRC |  | high | 0.80 |
| 3 | Abstract 6779: Molecular glue degraders of the RNA binding p | CRC |  | high | 0.80 |
| 4 | Abstract 1303: Real-world overall survival with trastuzumab  | CRC |  | high | 0.80 |
| 5 | Abstract 3311: Autologous colorectal cancer-mesothelial co-c | CRC |  | high | 0.80 |
| 6 | Abstract 7096: Design and development of a biparatopic antib | Colorectal |  | high | 0.50 |
| 7 | Abstract 4089: Immunosuppressive cellular topography and gen | Colorectal |  | high | 0.50 |
| 8 | Abstract 3824: Cross-cohort robust detection of colorectal c | CRC |  | high | 0.20 |
| 9 | Abstract 3021: ONM-421, a pH-responsive polymer-drug conjuga | CRC |  | high | 0.00 |
| 10 | Abstract 2210: Ouabain exerts both anti-tumor and senolytic  | Glioblastoma |  | high | 0.00 |
| 11 | Abstract 6453: Serum pharmacodynamic biomarkers of IPI/NIVO/ | melanoma |  | high | 0.00 |
| 12 | Abstract 1256: A novel therapeutic approach to overcome meta | Lung |  | high | 0.00 |
| 13 | Abstract 7172: Combination of B7-H4-TOP1i ADC, PD-1/TIGIT bi | CRC | IMMUNOTHERAPY|PARP_INHIBITORS | high | 0.00 |
| 14 | Abstract 6076: Preclinical models of cancer cachexia: Bridgi | None |  | high | 0.00 |
| 15 | Abstract 1136: Analytical evaluation of a whole genome tumor | Bladder, Breast, Melanoma, NSCLC, CRC |  | high | 0.00 |
| 16 | Abstract 4504: Synergistic antitumor activity of the MTA-coo | None | PRMT5_MTAP | high | 0.00 |
| 17 | Abstract 7077: ISM1745, an MTA-cooperative PRMT5 inhibitor f | CRC | PRMT5_MTAP | high | 0.00 |
| 18 | Abstract 5825: Discovery of FL-261 as a theranostic RDC vect | Lung |  | high | 0.00 |
| 19 | Abstract 5668: USP7 and ITCH/UBE4B ubiquitin ligase regulate | neuroblastoma |  | high | 0.00 |
| 20 | Abstract 7083: ABSK211, a highly potent and orally available | Pancreatic, Lung, Colorectal |  | high | 0.00 |

---

## 6. Key Intelligence for CrisPRO Positioning

### Unoccupied White Space (Brenus Advantage)
1. **MSS CRC + whole-cell vaccine:** 0 competitor abstracts. STC-1010 (CT051) is the only entry.
2. **MBD4_LOF + ATR_WEE1 SL axis:** 0 MBD4_LOF abstracts — mechanistic angle not on AACR radar.
3. **MSS CRC + replication stress:** 0 abstracts — positioning angle uncontested.

### Active Competitive Axes (Monitor)
- **PARP_INHIBITORS:** 366 abstracts — crowded, Ovarian/Breast dominant
- **IMMUNOTHERAPY:** 1515 abstracts — very crowded, Lung/Breast dominant  
- **PI3K_AKT:** 382 abstracts — active, Breast/CRC/Lung
- **PRMT5_MTAP:** 50 abstracts — emerging Phase I data (AMG 193, MRTX1719)

### CRC-Specific Competitive Intelligence
- **WRN (CRC n=8):** MSI-H SL angle; VX-803/HRO-761 data — distinct from MSS
- **PRMT5_MTAP (CRC n=7):** MTAP-deleted CRC subpopulation
- **ATR_WEE1 (CRC n=13):** Replication stress angle — closest to Brenus mechanism

---

*Generated by CrisPRO AACR 2026 Enrichment Pipeline v2 | 2026-07-04 17:25 UTC*
