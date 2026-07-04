"""
Literature Receipt Layer — frozen evidence from validated publications.

The cytidine-analog axis for MBD4 is the calibration gold standard.
Its receipts are frozen here as authoritative ground-truth.
All other axes are evaluated against this bar.

Frozen receipts come from:
  - npj Precis Oncol 2022, PMID 35428381 (MBD4-KO + cytidine SL)
  - Additional curated records that can be extended at runtime.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from .models import (
    CandidateAxis,
    ModalityEvidence,
    ModalityStatus,
    Modality,
)


# ─────────────────────────────────────────────────────────────────────────────
# FROZEN RECEIPT STORE
# Each entry is keyed by (gene, axis) and provides pre-populated
# ModalityEvidence objects for in_vitro, in_vivo, and clinical modalities.
# These come from published, peer-reviewed literature and are not recalculated.
# ─────────────────────────────────────────────────────────────────────────────

# Type alias
_ReceiptKey = tuple  # (gene_upper, CandidateAxis)
_ReceiptStore = Dict[_ReceiptKey, Dict[str, ModalityEvidence]]

_FROZEN_RECEIPTS: _ReceiptStore = {

    # ── MBD4 + Cytidine Analogs ───────────────────────────────────────────────
    # Gold standard calibration axis.
    # Source: Ohashi et al., npj Precis Oncol 2022 (PMID 35428381)
    ("MBD4", CandidateAxis.CYTIDINE_ANALOGS): {

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            delta_ic50_log2=-2.3,    # ~5× more sensitive (MBD4-KO vs WT)
            n_mut=3,
            n_wt=3,
            summary=(
                "MBD4-KO isogenic cell lines show ~5× increased sensitivity to gemcitabine "
                "and cytarabine vs WT in rescue experiments. Effect is MBD4-specific: "
                "re-expression of MBD4 (but not catalytic-dead mutant) restores resistance. "
                "Two independent isogenic systems validated."
            ),
            pmids=["35428381"],
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MBD4-deficient uveal melanoma PDX derived from the index patient "
                "showed dramatic, durable response to gemcitabine monotherapy. "
                "PDX response mirrors patient response; tumor regression observed. "
                "WT-MBD4 PDX controls did not respond."
            ),
            pmids=["35428381"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "Index patient (uveal melanoma, MBD4 germline LOF): dramatic, durable "
                "response to gemcitabine — exceptional for uveal melanoma, which is "
                "normally chemo-refractory. Response mirrors PDX and isogenic model data. "
                "Additional MBD4-deficient colorectal patients from MSI cohort showed "
                "preferential cytidine-analog response vs MBD4-WT."
            ),
            pmids=["35428381"],
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MBD4 LOF causes BER substrate accumulation (dU:dG mispairs). "
                "Gemcitabine and cytarabine exploit this: dFdCTP/araC incorporation "
                "into dU-rich repair patches stalls replication. Expression of BER "
                "components (POLB, PCNA) is elevated in MBD4-KO cells, consistent "
                "with chronic BER stress."
            ),
            pmids=["35428381"],
        ),
    },

    # ── MBD4 + PARP Inhibitors ────────────────────────────────────────────────
    # Mechanistic signal via RNF144A–PARP1 axis; CRISPR negative.
    # Sources: Ohashi 2022 (PMID 35428381) + RNF144A-PARP1 paper
    ("MBD4", CandidateAxis.PARP_INHIBITORS): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MBD4 LOF → BER strand-break accumulation → compensatory PARP1 "
                "upregulation observed in some MBD4-deficient models. RNF144A "
                "ubiquitinates PARP1 for degradation; low RNF144A (seen in subset "
                "of MBD4-LOF tumors) → elevated PARP1 → theoretical trapping substrate. "
                "This creates an expression-level hook, not a dependency-level hook."
            ),
            pmids=["35428381"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MIXED,
            summary=(
                "No robust isogenic KO/WT SL validated for PARP inhibitors in MBD4-deficient "
                "cells as of current literature. RNF144A knockdown in MDA-MB-231 (BRCA-WT) "
                "sensitizes to olaparib, but this is PARP1 upregulation effect, not MBD4 "
                "direct SL. MBD4-specific PARP inhibitor in vitro data absent."
            ),
            pmids=["35428381"],
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.MISSING,
            summary="No MBD4-deficient PDX data for PARP inhibitor response.",
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No published patient-level response data for PARP inhibitors specifically "
                "in MBD4-deficient tumors. PARP inhibitor trials in MSI/hypermutated contexts "
                "have mixed results and do not control for MBD4 status."
            ),
        ),
    },

    # ── MBD4 + ATR/WEE1 ──────────────────────────────────────────────────────
    # BER stress → replication stress → ATR/WEE1 checkpoint dependency
    # TIER UPGRADE: MECHANISTIC → STRONG (2026-07-03)
    # Source: canonical GDSC2 ceralasertib rerun (fjkiani/crispro,
    #   publications/00-mbd4-manuscript/mbd4_parp_response/
    #   artifacts/canonical_atr_wee1_rerun_20260405/)
    # Evidence: n_LOF=14 vs n_WT=914, Δ LN_IC50=−0.732, p=0.021, d=−0.503
    #   AUC: Δ=−0.056, p=0.013, d=−0.557
    #   4 confound stress tests all p<0.05 (MSI purge, TP53 co-strat, non-Bowel, LOO)
    # Tier logic: pharma_pos=True (gdsc=POSITIVE) + expression=POSITIVE → n_positive=2
    #   → "Strong candidate dependency axis" per _assign_recommendation_tier()
    ("MBD4", CandidateAxis.ATR_WEE1): {

        "gdsc": ModalityEvidence(
            modality=Modality.GDSC_PHARMACOLOGIC,
            status=ModalityStatus.POSITIVE,
            delta_ic50_log2=None,          # LN scale: Δ LN_IC50=−0.732
            delta_auc=-0.056,
            p_value=0.021,
            effect_size=-0.503,            # Cohen's d (primary LN_IC50)
            n_mut=14,
            n_wt=914,
            drug_screen_dataset="GDSC2",
            stratifier="MBD4_LOF_vs_WT",
            summary=(
                "Ceralasertib (AZD6738, ATRi): MBD4-LOF cell lines (n=14) vs WT (n=914). "
                "Primary: Δ LN_IC50=−0.732, p=0.021, Cohen's d=−0.503. "
                "AUC: Δ=−0.056, p=0.013, d=−0.557. "
                "4 confound stress tests all p<0.05: MSI purge (p=0.015, d=−0.623), "
                "TP53 co-stratification (p=0.003, d=−0.740), non-Bowel (p=0.025, d=−0.599), "
                "LOO worst-case (p=0.045). "
                "Adavosertib (WEE1i): Δ LN_IC50=−0.508, p=0.074 (trend only, not significant). "
                "Receipt locked: canonical_atr_wee1_rerun_20260405."
            ),
            notes=(
                "Upgrade criterion met: BH-adjusted p<0.05 AND Cohen's d<−0.5 for ceralasertib. "
                "Adavosertib is trend-only and does NOT independently meet the upgrade bar. "
                "PARPi axis FALSIFIED: RNF144A not downregulated (p=0.476), PARP1 not upregulated (p=0.605)."
            ),
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "Mechanistic chain: MBD4-LOF → chronic BER substrate accumulation (dU:dG mispairs) "
                "→ elevated ssDNA at stalled replication forks → constitutive ATR activation "
                "→ ATR/WEE1 checkpoint dependency. "
                "PARP1 pan-cancer expression correlate: ρ=−0.416, n=481, p=1.36×10⁻²¹ "
                "(publishable signal; PARP1 as pan-cancer PARPi correlate, not MBD4-direct SL)."
            ),
            pmids=["35428381"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No published MBD4-specific isogenic ATR/WEE1 inhibitor in vitro data. "
                "Path to upgrade: HAP1-MBD4-KO isogenic ATRi/WEE1i IC50 assay "
                "(ceralasertib + berzosertib vs parental) + γH2AX/pCHK1 western blot."
            ),
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.MISSING,
            summary=(
                "No MBD4-specific ATR/WEE1 PDX data. "
                "Path to upgrade: MBD4-KO organoid ATRi dose-response (patient-derived, UCEC preferred)."
            ),
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No MBD4-specific clinical ATR/WEE1 inhibitor data. "
                "Path to upgrade: PATRIOT data access request; prospective UCEC trial with MBD4-LOF enrichment arm."
            ),
        ),
    },

    # ── MBD4 + WRN ───────────────────────────────────────────────────────────
    # WRN is an MSI-H SL target; MBD4 LOF is associated with hypermutation
    # but the SL is MSI-context-dependent, not MBD4-direct.
    ("MBD4", CandidateAxis.WRN): {

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "WRN SL is established for MSI-H context (Chan et al. Nature 2019). "
                "MBD4 LOF → hypermutator phenotype → subset of tumors become MSI-H. "
                "WRN SL with MBD4-direct LOF (MSI-independent) not established in vitro."
            ),
            pmids=["30971823"],  # Chan 2019 Nature
            notes="Signal is MSI-H specific; MBD4 LOF alone insufficient.",
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary="No clinical WRN inhibitor data for MBD4-deficient tumors.",
        ),
    },

    # ── CCNE1 + PKMYT1 ──────────────────────────────────────────────────────
    # Sanger/DepMap: CCNE1-amp → PKMYT1 dependency; RP-6306 Phase I ongoing
    ("CCNE1", CandidateAxis.PKMYT1): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.POSITIVE,
            summary=(
                "Sanger/DepMap CCNE1-amp → PKMYT1 dependency."
            ),
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "CCNE1 amplification associated with PKMYT1 dependency in expression data."
            ),
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "CCNE1-amplified cell lines show sensitivity to PKMYT1 inhibition in vitro."
            ),
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary="RP-6306 Phase I ongoing — no clinical data yet.",
        ),
    },

    # ── MBD4 + PKMYT1 ────────────────────────────────────────────────────────
    # MBD4 has no direct PKMYT1 evidence — all MISSING
    ("MBD4", CandidateAxis.PKMYT1): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific PKMYT1 CRISPR dependency data.",
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific PKMYT1 expression association data.",
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific PKMYT1 in vitro data.",
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific PKMYT1 clinical data.",
        ),
    },

    # ── MBD4 + Immunotherapy ─────────────────────────────────────────────────
    # Hypermutator phenotype → high TMB → IO response
    ("MBD4", CandidateAxis.IMMUNOTHERAPY): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MBD4 LOF → CpG→TpG hypermutator clock → highest TMB in many tumor types. "
                "High TMB is associated with IO response (pan-cancer FDA approval). "
                "TMB > 10 mut/Mb threshold frequently exceeded in MBD4-deficient tumors."
            ),
            pmids=["35428381"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MIXED,
            summary=(
                "Exceptional IO responses reported in MBD4-germline-LOF uveal melanoma "
                "and colorectal cancer cases. Pan-cancer TMB-high FDA approval supports "
                "IO use. Not a direct SL axis; relies on hypermutator → TMB mechanism. "
                "Small n — case reports, not controlled trial data."
            ),
            pmids=["35428381"],
        ),
    },
    # ══════════════════════════════════════════════════════════════════════════
    # BREAST CANCER EXPANSION — 2026-07-03
    # Tier hierarchy: VALIDATED > STRONG > MECHANISTIC
    # Approved: Fahad Kiani 2026-07-03
    # ══════════════════════════════════════════════════════════════════════════

    # ── BRCA1 + PARP Inhibitors ───────────────────────────────────────────────
    # VALIDATED: Phase III clinical data (OlympiAD), isogenic in vitro, CRISPR
    # Source: Robson et al. NEJM 2017 (OlympiAD, PMID 28578601)
    ("BRCA1", CandidateAxis.PARP_INHIBITORS): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA1-deficient cells show strong PARP1 dependency in DepMap CRISPR screens. "
                "BRCA1-LOF → HR deficiency → reliance on PARP-mediated SSB repair. "
                "Multiple independent CRISPR datasets confirm dependency."
            ),
            pmids=["28578601"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            delta_ic50_log2=-3.2,
            summary=(
                "BRCA1-deficient cell lines (HCC1937, SUM149PT) show >10× sensitivity "
                "to olaparib vs BRCA1-WT controls. Isogenic BRCA1-KO rescue experiments "
                "confirm BRCA1-specific effect. Validated across multiple PARP inhibitors "
                "(olaparib, niraparib, talazoparib, rucaparib)."
            ),
            pmids=["28578601", "18772396"],
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA1-mutant breast cancer PDX models show durable tumor regression "
                "with olaparib monotherapy. Response mirrors clinical outcomes. "
                "BRCA1-WT PDX controls do not respond."
            ),
            pmids=["28578601"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "OlympiAD (Phase III, n=302): olaparib vs chemotherapy in BRCA1/2-mutant "
                "HER2-negative metastatic breast cancer. PFS: 7.0 vs 4.2 months (HR=0.58, p<0.001). "
                "FDA approved olaparib for BRCA1/2-mutant HER2-negative metastatic breast cancer. "
                "EMBRACA (Phase III, n=431): talazoparib PFS 8.6 vs 5.6 months (HR=0.54, p<0.001)."
            ),
            pmids=["28578601", "29863979"],
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA1 LOF → HR deficiency → BRCAness signature. "
                "PARP1 trapping mechanism: BRCA1-deficient cells cannot resolve PARP1-DNA complexes. "
                "RAD51 foci formation impaired in BRCA1-LOF cells (HR deficiency marker)."
            ),
            pmids=["28578601"],
        ),
    },

    # ── BRCA2 + PARP Inhibitors ───────────────────────────────────────────────
    # VALIDATED: Phase III clinical data (OlympiAD includes BRCA2), isogenic in vitro
    ("BRCA2", CandidateAxis.PARP_INHIBITORS): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA2-deficient cells show strong PARP1 dependency in DepMap CRISPR screens. "
                "BRCA2-LOF → HR deficiency → PARP dependency. "
                "Confirmed across multiple cancer lineages."
            ),
            pmids=["28578601"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            delta_ic50_log2=-3.5,
            summary=(
                "BRCA2-deficient cell lines (CAPAN-1, PE01) show >10× sensitivity to olaparib. "
                "Isogenic BRCA2-KO rescue experiments confirm BRCA2-specific effect. "
                "Validated across olaparib, niraparib, talazoparib, rucaparib."
            ),
            pmids=["28578601", "18772396"],
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA2-mutant breast and ovarian cancer PDX models show durable response "
                "to PARP inhibitors. Validated in multiple independent PDX cohorts."
            ),
            pmids=["28578601"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "OlympiAD (Phase III): BRCA2-mutant subgroup shows PFS benefit with olaparib. "
                "EMBRACA (Phase III): talazoparib PFS benefit in BRCA2-mutant breast cancer. "
                "FDA approved olaparib and talazoparib for BRCA1/2-mutant HER2-negative metastatic breast cancer."
            ),
            pmids=["28578601", "29863979"],
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "BRCA2 LOF → HR deficiency → BRCAness signature. "
                "RAD51 foci formation impaired (HR deficiency marker). "
                "PARP1 trapping mechanism confirmed in BRCA2-deficient models."
            ),
            pmids=["28578601"],
        ),
    },

    # ── PIK3CA + PI3K/AKT Inhibitors ─────────────────────────────────────────
    # STRONG: Phase III clinical data (SOLAR-1), in vitro validated
    # Source: André et al. NEJM 2019 (SOLAR-1, PMID 31091374)
    ("PIK3CA", CandidateAxis.PI3K_AKT): {

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "PIK3CA-mutant breast cancer cell lines show selective sensitivity to "
                "alpelisib (BYL719) vs PIK3CA-WT controls. "
                "PIK3CA-mutant lines: MCF7, T47D, ZR-75-1 — all show >5× sensitivity. "
                "Validated across alpelisib, capivasertib, inavolisib."
            ),
            pmids=["31091374"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "SOLAR-1 (Phase III, n=572): alpelisib + fulvestrant vs placebo + fulvestrant "
                "in PIK3CA-mutant HR+/HER2- advanced breast cancer. "
                "PFS: 11.0 vs 5.7 months (HR=0.65, p<0.001). FDA approved alpelisib (Piqray) "
                "for PIK3CA-mutant HR+/HER2- advanced breast cancer. "
                "CAPItello-291 (Phase III): capivasertib + fulvestrant PFS benefit in PIK3CA/AKT1/PTEN-altered BC."
            ),
            pmids=["31091374"],
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "PIK3CA hotspot mutations (H1047R, E545K, E542K) activate PI3K/AKT/mTOR pathway. "
                "pAKT and pS6K1 elevated in PIK3CA-mutant tumors. "
                "PIK3CA mutation frequency: ~40% of HR+/HER2- breast cancers."
            ),
            pmids=["31091374"],
        ),

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.MIXED,
            summary=(
                "PIK3CA is a gain-of-function oncogene — CRISPR KO of PIK3CA itself reduces "
                "pathway activity and is not a dependency readout. "
                "Downstream effectors AKT1 and mTOR show elevated dependency scores in "
                "PIK3CA-mutant lines in DepMap, but this is an indirect/context dependency, "
                "not a direct PIK3CA SL. Marked MIXED: pathway dependency exists but "
                "PIK3CA-direct CRISPR SL is not established."
            ),
            pmids=["31091374"],
            notes="Tier driven by clinical=POSITIVE (SOLAR-1) + in_vitro=POSITIVE, not CRISPR.",
        ),
    },

    # ── PALB2 + PARP Inhibitors ───────────────────────────────────────────────
    # STRONG: Phase II clinical data (Tung 2020 NEJM), in vitro validated
    # Source: Tung et al. NEJM 2020 (PMID 32023372)
    ("PALB2", CandidateAxis.PARP_INHIBITORS): {

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "PALB2-deficient cells show sensitivity to PARP inhibitors comparable to BRCA2-LOF. "
                "PALB2 is a BRCA2 binding partner; PALB2-LOF → HR deficiency → PARP dependency. "
                "Isogenic PALB2-KO cell lines confirm PARP inhibitor sensitivity."
            ),
            pmids=["33119476"],  # Tung et al. JCO 2020 (DOI: 10.1200/jco.20.02151)
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "TBCRC 048 (Phase II, n=54): olaparib in PALB2-mutant breast cancer. "
                "ORR: 75% in germline PALB2-mutant cohort. "
                "Tung et al. NEJM 2020 (PMID 32023372). "
                "Supports PALB2 as a PARP inhibitor biomarker beyond BRCA1/2."
            ),
            pmids=["33119476"],  # Tung et al. JCO 2020 (DOI: 10.1200/jco.20.02151)
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "PALB2 LOF → HR deficiency → BRCAness signature. "
                "PALB2 bridges BRCA1 and BRCA2 in HR repair complex. "
                "RAD51 foci formation impaired in PALB2-LOF cells."
            ),
            pmids=["33119476"],  # Tung et al. JCO 2020 (DOI: 10.1200/jco.20.02151)
        ),
    },

    # ── ARID1A + ATR/WEE1 ────────────────────────────────────────────────────
    # STRONG: CRISPR dependency data (Williamson 2016), in vitro validated
    # Source: Williamson et al. Cancer Cell 2016 (PMID 27070702)
    ("ARID1A", CandidateAxis.ATR_WEE1): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.POSITIVE,
            summary=(
                "ARID1A-LOF → ATR dependency confirmed in DepMap CRISPR screens. "
                "ARID1A-deficient cells show elevated ATR dependency (Chronos delta < -0.3). "
                "Williamson et al. 2016: ARID1A-LOF → chromatin remodeling defect → "
                "fork instability → ATR checkpoint dependency."
            ),
            pmids=["27958275"],  # Williamson et al. Nat Commun 2016
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "ARID1A-deficient cell lines (OVISE, TOV21G) show sensitivity to "
                "ATR inhibitors (berzosertib, ceralasertib) vs ARID1A-WT controls. "
                "Williamson 2016: ARID1A-LOF → ATRi sensitivity in ovarian and breast cancer lines. "
                "Effect is ARID1A-specific: re-expression of ARID1A restores resistance."
            ),
            pmids=["27958275"],  # Williamson et al. Nat Commun 2016 (DOI: 10.1038/ncomms13837)
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "ARID1A LOF → SWI/SNF chromatin remodeling defect → fork stalling at "
                "nucleosome-dense regions → elevated γH2AX and RPA32 phosphorylation. "
                "RS markers elevated in ARID1A-LOF cells. "
                "ARID1A mutation frequency: ~10% of breast cancers (TCGA)."
            ),
            pmids=["27958275"],  # Williamson et al. Nat Commun 2016
        ),
    },

    # ── RB1 + ATR/WEE1 ───────────────────────────────────────────────────────
    # MECHANISTIC: expression/pathway only — no CRISPR or in vitro data yet
    ("RB1", CandidateAxis.ATR_WEE1): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "RB1 LOF → E2F transcription factor deregulation → premature S-phase entry "
                "→ replication stress. E2F targets include RRM2, TYMS (nucleotide synthesis) "
                "and CDC6, MCM2-7 (origin licensing) — all elevated in RB1-LOF cells. "
                "Elevated RS markers (γH2AX, pRPA32) in RB1-deficient breast cancer lines. "
                "RB1 mutation frequency: ~20% of triple-negative breast cancer (TCGA)."
            ),
            pmids=["25490446"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No published isogenic RB1-KO ATRi sensitivity data in breast cancer. "
                "Path to upgrade: RB1-KO MCF10A or MDA-MB-231 ATRi dose-response."
            ),
        ),
    },

    # ── FBXW7 + PKMYT1 ───────────────────────────────────────────────────────
    # MECHANISTIC: FBXW7-LOF → cyclin E stabilization → CCNE1-like RS → PKMYT1 dependency
    ("FBXW7", CandidateAxis.PKMYT1): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "FBXW7 is the E3 ubiquitin ligase that degrades cyclin E (CCNE1). "
                "FBXW7 LOF → cyclin E stabilization → CCNE1-like replication stress. "
                "CCNE1 overexpression → PKMYT1 dependency (established). "
                "FBXW7-LOF tumors phenocopy CCNE1-amplified tumors in RS signature. "
                "FBXW7 mutation frequency: ~5% of breast cancers (TCGA)."
            ),
            pmids=["19360079"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No published FBXW7-specific PKMYT1 inhibitor data. "
                "Hypothesis: FBXW7-LOF → CCNE1 stabilization → RP-6306 sensitivity. "
                "Path to upgrade: FBXW7-KO cell line RP-6306 dose-response."
            ),
        ),
    },

    # ── CDH1 + ATR/WEE1 ──────────────────────────────────────────────────────
    # MECHANISTIC: CDH1-LOF → EMT → RS — novel hypothesis
    ("CDH1", CandidateAxis.ATR_WEE1): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "CDH1 (E-cadherin) LOF → epithelial-to-mesenchymal transition (EMT). "
                "EMT → cytoskeletal reorganization → altered nuclear architecture → "
                "replication fork instability. EMT cells show elevated RS markers. "
                "CDH1 mutation frequency: ~15% of lobular breast cancer (TCGA). "
                "Hypothesis: CDH1-LOF → EMT-driven RS → ATR/WEE1 checkpoint dependency."
            ),
            notes="Novel hypothesis — no direct experimental validation yet.",
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "No published CDH1-specific ATRi sensitivity data. "
                "Path to upgrade: CDH1-KO MCF7 or T47D ATRi dose-response + γH2AX quantification."
            ),
        ),
    },

    # ── PIK3R1 + PI3K/AKT ────────────────────────────────────────────────────
    # MECHANISTIC: PIK3R1-LOF → PI3K hyperactivation → PI3K/AKT dependency
    ("PIK3R1", CandidateAxis.PI3K_AKT): {

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "PIK3R1 (p85α) is the regulatory subunit of PI3K. "
                "PIK3R1 LOF → loss of PI3K negative regulation → PI3K/AKT hyperactivation. "
                "PIK3R1-LOF tumors phenocopy PIK3CA-mutant tumors in PI3K pathway activation. "
                "pAKT elevated in PIK3R1-LOF breast cancer lines. "
                "PIK3R1 mutation frequency: ~8% of breast cancers (TCGA)."
            ),
            pmids=["31091374"],
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "Limited PIK3R1-specific PI3Ki sensitivity data. "
                "Some evidence that PIK3R1-LOF lines respond to alpelisib similarly to PIK3CA-mutant. "
                "Path to upgrade: PIK3R1-KO isogenic alpelisib/capivasertib dose-response."
            ),
        ),
    },


    # ── MTAP-deleted tumors + PRMT5 inhibitors ───────────────────────────────
    # MTAP deletion (9p21.3) → PRMT5 synthetic lethality
    # Drugs: BMS-986504, MRTX1719, AMG 193, IDE397
    # Source: AACR 2026 abstract 1558 + published PRMT5/MTAP SL literature
    ("MTAP", CandidateAxis.PRMT5_MTAP): {

        "crispr": ModalityEvidence(
            modality=Modality.CRISPR_DEPENDENCY,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MTAP deletion creates PRMT5 dependency via MTA accumulation. "
                "CRISPR screens confirm PRMT5 as essential in MTAP-deleted lines. "
                "Multiple independent DepMap/Sanger datasets confirm dependency."
            ),
            pmids=["28783718", "28783719"],  # Mavrakis 2016, Kryukov 2016
        ),

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MTAP-deleted cell lines show selective sensitivity to PRMT5 inhibitors "
                "(BMS-986504, MRTX1719, AMG 193, IDE397) vs MTAP-WT controls. "
                "MTA accumulation in MTAP-deleted cells partially inhibits PRMT5, "
                "creating hypersensitivity to further PRMT5 inhibition."
            ),
            pmids=["28783718"],
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MIXED,
            summary=(
                "Multiple Phase I/II trials ongoing: BMS-986504, MRTX1719 (NCT04089449), "
                "AMG 193 (NCT04888312), IDE397 (NCT04794699). "
                "Early clinical signals in MTAP-deleted solid tumors. "
                "No Phase III data yet — MIXED pending mature readouts."
            ),
        ),

        "expression": ModalityEvidence(
            modality=Modality.EXPRESSION_ASSOC,
            status=ModalityStatus.POSITIVE,
            summary=(
                "MTAP deletion (9p21.3) co-occurs with CDKN2A deletion in ~15% of solid tumors. "
                "MTA accumulation measurable by metabolomics in MTAP-deleted tumors. "
                "PRMT5 protein levels elevated in MTAP-deleted vs WT tumors (compensatory upregulation)."
            ),
        ),
    },


}


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def get_literature_receipts(
    gene: str,
    axis: CandidateAxis,
) -> Dict[str, ModalityEvidence]:
    """
    Return frozen literature-backed ModalityEvidence cells for (gene, axis).
    Keys are modality names: "in_vitro", "in_vivo", "clinical", "expression".
    Returns empty dict if no frozen data exists.
    """
    return dict(_FROZEN_RECEIPTS.get((gene.upper(), axis), {}))


def get_calibration_narrative(gene: str = "MBD4") -> str:
    """
    Return a structured narrative for the gold-standard calibration axis
    (cytidine analogs for MBD4). Used to define what "high evidence SL" looks like.
    """
    receipts = get_literature_receipts(gene, CandidateAxis.CYTIDINE_ANALOGS)
    if not receipts:
        return (
            f"No frozen calibration receipts available for {gene} / cytidine analogs."
        )

    lines = [
        f"=== CALIBRATION GOLD STANDARD: {gene} + Cytidine Analogs ===",
        "",
        "This is the bar that any SL axis must approach to be classified as 'Validated.'",
        "",
        "Evidence across modalities:",
    ]
    for mod_name, ev in receipts.items():
        pmid_str = " [PMID: " + ", ".join(ev.pmids) + "]" if ev.pmids else ""
        lines.append(f"  [{mod_name.upper()}] {ev.status.value.upper()}: {ev.summary}{pmid_str}")

    lines += [
        "",
        "Pattern that defines 'High Evidence / Validated SL lever':",
        "  • In vitro isogenic KO/WT with rescue validation",
        "  • In vivo PDX matching the in vitro sensitivity",
        "  • At least one patient-level clinical receipt",
        "  • Mechanistically coherent (not MSI confound alone)",
        "",
        "PARP inhibitors for MBD4 do NOT yet meet this bar.",
        "ATR/WEE1 inhibitors for MBD4 are mechanistically plausible but lack in vitro/in vivo receipts.",
        "WRN inhibitors are MSI-H context-dependent and should NOT be declared MBD4-direct SL.",
    ]
    return "\n".join(lines)


def list_receipts_for_gene(gene: str) -> Dict[str, List[str]]:
    """Return a summary of which axes have frozen receipts for a gene."""
    gene_upper = gene.upper()
    result = {}
    for (g, axis), cells in _FROZEN_RECEIPTS.items():
        if g == gene_upper:
            positive = [k for k, v in cells.items() if v.status == ModalityStatus.POSITIVE]
            result[axis.value] = positive
    return result
