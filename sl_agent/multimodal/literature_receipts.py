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
    ("MBD4", CandidateAxis.ATR_WEE1): {

        "in_vitro": ModalityEvidence(
            modality=Modality.IN_VITRO_FUNCTIONAL,
            status=ModalityStatus.MISSING,
            summary=(
                "Mechanistically plausible: chronic BER stress in MBD4-LOF → elevated "
                "ssDNA at stalled forks → ATR activation. No published MBD4-specific "
                "isogenic ATR/WEE1 inhibitor data."
            ),
        ),

        "in_vivo": ModalityEvidence(
            modality=Modality.IN_VIVO_PDX,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific ATR/WEE1 PDX data.",
        ),

        "clinical": ModalityEvidence(
            modality=Modality.CLINICAL,
            status=ModalityStatus.MISSING,
            summary="No MBD4-specific clinical ATR/WEE1 inhibitor data.",
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
