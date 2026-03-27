"""
Matrix Builder — assembles the full multi-modal evidence matrix for a gene.

Combines:
  - CRISPR SL engine results (from sl_engine.py)
  - Pharmacologic screen analysis (PRISM + GDSC via pharmacologic_analyzer.py)
  - Literature receipts (frozen validated data via literature_receipts.py)
  - Open KB evidence (CIViC + CGI + JAX via kb_engine.py)

Then invokes the modality fuser to derive evidence levels and tiers.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import pandas as pd

from .models import (
    CandidateAxis,
    EvidenceMatrix,
    EvidenceRow,
    ModalityEvidence,
    ModalityStatus,
    Modality,
    MultiModalQueryInput,
)
from .pharmacologic_analyzer import analyze_drug_screen, aggregate_by_axis
from .literature_receipts import get_literature_receipts, get_calibration_narrative
from .modality_fuser import fuse_matrix, build_agreement_report, build_updated_narrative
from .replication_stress import compute_rs_score, assess_combinatorial_vulnerability

logger = logging.getLogger(__name__)

# ── Axis metadata registry ────────────────────────────────────────────────────

_AXIS_META: Dict[CandidateAxis, Dict] = {
    CandidateAxis.CYTIDINE_ANALOGS: {
        "label": "Cytidine Analogs (gemcitabine, cytarabine)",
        "mechanism": (
            "BER substrate accumulation in MBD4-LOF (dU:dG mispairs). Gemcitabine / "
            "cytarabine triphosphate incorporation into BER repair patches stalls "
            "replication. MBD4-specific: rescue by WT-MBD4 but not catalytic-dead mutant."
        ),
    },
    CandidateAxis.PARP_INHIBITORS: {
        "label": "PARP Inhibitors (olaparib, niraparib, talazoparib, rucaparib)",
        "mechanism": (
            "BER strand-break accumulation → compensatory PARP1 upregulation. "
            "RNF144A–PARP1 axis: low RNF144A → elevated PARP1 → trapping substrate. "
            "Theoretical SSB → DSB conversion under PARP trapping. "
            "Context-dependent: only in MBD4-LOF tumors with high PARP1 / low RNF144A."
        ),
    },
    CandidateAxis.ATR_WEE1: {
        "label": "ATR / WEE1 Inhibitors (berzosertib, adavosertib, ceralasertib)",
        "mechanism": (
            "BER stress → ssDNA accumulation at stalled forks → ATR kinase activation. "
            "ATR/CHK1 and WEE1/CDK1 checkpoints are required to resolve replication stress. "
            "MBD4-LOF BER load may create constitutive low-level replication stress, "
            "sensitizing to checkpoint abrogation especially in TP53-loss context."
        ),
    },
    CandidateAxis.WRN: {
        "label": "WRN Helicase Inhibitors (VX-803, MRTX1719)",
        "mechanism": (
            "WRN SL is established for MSI-H context (reversed replication forks, "
            "aberrant AT-rich secondary structures). MBD4 LOF → hypermutator → subset "
            "of tumors become MSI-H. WRN SL with MBD4-direct LOF (MSI-independent) "
            "not yet validated. Signal may be MSI-H driven, not MBD4-specific."
        ),
    },
    CandidateAxis.IMMUNOTHERAPY: {
        "label": "Immunotherapy / Checkpoint Inhibitors (PD-1/PD-L1)",
        "mechanism": (
            "MBD4 LOF → CpG→TpG hypermutator clock → highest TMB in many tumor types. "
            "High TMB → increased neoantigen burden → IO response. "
            "Pan-cancer TMB-high FDA approval supports pembrolizumab use. "
            "Not a direct SL axis; relies on hypermutation → neoantigens mechanism."
        ),
    },
    CandidateAxis.PKMYT1: {
        "label": "PKMYT1 Kinase Inhibitors (RP-6306 class)",
        "mechanism": (
            "PKMYT1 kinase inhibition (RP-6306 class); replication stress axis. "
            "CCNE1-amplified tumors depend on PKMYT1 for S-phase checkpoint regulation. "
            "PKMYT1 inhibition causes premature mitotic entry and replication catastrophe "
            "in replication-stressed cells."
        ),
    },
}


def build_evidence_matrix(
    query: MultiModalQueryInput,
    sl_partners: Optional[List] = None,   # List[SLPartner] from sl_engine
    mutant_ids: Optional[List[str]] = None,
    wt_ids: Optional[List[str]] = None,
    prism_df: Optional[pd.DataFrame] = None,
    prism_meta: Optional[pd.DataFrame] = None,
    gdsc_df: Optional[pd.DataFrame] = None,
    sample_info: Optional[pd.DataFrame] = None,
    depmap_release: str = "unknown",
) -> EvidenceMatrix:
    """
    Build the complete multi-modal evidence matrix.

    Steps:
    1. For each candidate axis, instantiate an EvidenceRow.
    2. Fill CRISPR cells from sl_partners (if available).
    3. Fill pharmacologic cells from PRISM/GDSC analysis.
    4. Fill in_vitro / in_vivo / clinical cells from literature receipts.
    5. Fuse to derive evidence levels and recommendation tiers.
    """
    gene = query.gene.upper()
    warnings: List[str] = []

    # ── Step 1: Instantiate rows ──────────────────────────────────────────────
    axes_to_build = query.axes if query.axes else list(CandidateAxis)
    rows: List[EvidenceRow] = []
    for axis in axes_to_build:
        if axis == CandidateAxis.CUSTOM:
            continue
        meta = _AXIS_META.get(axis, {})
        rows.append(EvidenceRow(
            axis=axis,
            axis_label=meta.get("label", axis.value),
            mechanism=meta.get("mechanism", "Unknown"),
        ))

    # ── Step 2: Fill CRISPR cells ─────────────────────────────────────────────
    if sl_partners:
        crispr_gene_map = _build_crispr_map(sl_partners)
        for row in rows:
            axis_crispr = _axis_to_target_genes(row.axis)
            best = _best_crispr_result(crispr_gene_map, axis_crispr)
            if best is not None:
                delta, pval, fdr, n_mut, n_wt, d = best
                if fdr < 0.25 and delta < -0.05:
                    status = ModalityStatus.POSITIVE
                elif delta > 0.0:
                    status = ModalityStatus.NEGATIVE
                else:
                    status = ModalityStatus.MIXED
                row.crispr = ModalityEvidence(
                    modality=Modality.CRISPR_DEPENDENCY,
                    status=status,
                    delta_dep=round(delta, 4),
                    p_value=round(pval, 6),
                    fdr=round(fdr, 4),
                    effect_size=round(d, 3) if d else None,
                    n_mut=n_mut,
                    n_wt=n_wt,
                    summary=(
                        f"Δdep={delta:+.3f}, FDR={fdr:.3f}, n_mut={n_mut}, n_wt={n_wt}"
                    ),
                )
            else:
                # Axis-level targets not found in CRISPR results
                row.crispr = ModalityEvidence(
                    modality=Modality.CRISPR_DEPENDENCY,
                    status=ModalityStatus.MISSING,
                    summary="No CRISPR SL partners identified for this axis target genes.",
                )
    else:
        for row in rows:
            row.crispr = ModalityEvidence(
                modality=Modality.CRISPR_DEPENDENCY,
                status=ModalityStatus.MISSING,
                summary="CRISPR analysis not run or no data available.",
            )
        warnings.append("CRISPR analysis not provided; CRISPR cells are all MISSING.")

    # ── Step 3: Fill pharmacologic cells ─────────────────────────────────────
    pharma_results = []
    if query.include_pharmacologic_stratification and mutant_ids and wt_ids:
        target_axes_str = [row.axis.value for row in rows]
        try:
            pharma_results = analyze_drug_screen(
                gene=gene,
                mutant_ids=mutant_ids,
                wt_ids=wt_ids,
                prism_df=prism_df,
                prism_meta=prism_meta,
                gdsc_df=gdsc_df,
                sample_info=sample_info,
                target_axes=target_axes_str,
                stratify_by_msi=query.stratify_by_msi,
            )
        except Exception as e:
            logger.warning("Pharmacologic analysis failed: %s", e)
            warnings.append(f"Pharmacologic analysis error: {e}")

    prism_by_axis = aggregate_by_axis(pharma_results, stratifier="MBD4_LOF_vs_WT")
    gdsc_by_axis: Dict[str, ModalityEvidence] = {}  # separate GDSC pass if needed

    for row in rows:
        axis_key = row.axis.value
        if axis_key in prism_by_axis:
            row.prism = prism_by_axis[axis_key]
        elif prism_df is None:
            row.prism = ModalityEvidence(
                modality=Modality.PRISM_PHARMACOLOGIC,
                status=ModalityStatus.MISSING,
                summary="PRISM data not loaded (data not available in this deployment).",
            )
        else:
            row.prism = ModalityEvidence(
                modality=Modality.PRISM_PHARMACOLOGIC,
                status=ModalityStatus.NEGATIVE,
                summary=f"No significant PRISM signal for {row.axis_label} in MBD4-LOF vs WT.",
            )

        if axis_key in gdsc_by_axis:
            row.gdsc = gdsc_by_axis[axis_key]
        else:
            row.gdsc = ModalityEvidence(
                modality=Modality.GDSC_PHARMACOLOGIC,
                status=ModalityStatus.MISSING,
                summary="GDSC not loaded (optional enrichment — see orchestrator.py).",
            )

    # ── Step 4: Fill literature receipts ─────────────────────────────────────
    if query.include_literature_receipts:
        for row in rows:
            receipts = get_literature_receipts(gene, row.axis)
            if "in_vitro" in receipts:
                row.in_vitro = receipts["in_vitro"]
            if "in_vivo" in receipts:
                row.in_vivo = receipts["in_vivo"]
            if "clinical" in receipts:
                row.clinical = receipts["clinical"]
            if "expression" in receipts and row.expression.status == ModalityStatus.MISSING:
                row.expression = receipts["expression"]

    # ── Step 5: Build matrix object ───────────────────────────────────────────
    cal_narrative = None
    if query.include_calibration_narrative:
        cal_narrative = get_calibration_narrative(gene)

    matrix = EvidenceMatrix(
        query_gene=gene,
        cancer_type=query.cancer_type,
        depmap_release=depmap_release,
        rows=rows,
        gold_standard_summary=cal_narrative,
        notes=warnings,
    )

    # ── Step 6: Compute RS score if features provided ───────────────────────
    rs_score = compute_rs_score(query.rs_features) if query.rs_features else None

    # ── Step 7: Fuse — compute levels, tiers, agreements ─────────────────────
    matrix = fuse_matrix(matrix, rs_score=rs_score)
    matrix.rs_score = rs_score

    # ── Step 8: Assess combinatorial vulnerabilities ──────────────────────────
    if rs_score is not None and rs_score.score > 0.0:
        combo = assess_combinatorial_vulnerability(
            gene=gene,
            rs_score=rs_score,
            evidence_matrix=matrix,
        )
        if combo is not None:
            matrix.combination_vulnerabilities = [combo]

    return matrix


def run_multimodal_analysis(
    query: MultiModalQueryInput,
    sl_partners=None,
    mutant_ids=None,
    wt_ids=None,
    prism_df=None,
    prism_meta=None,
    gdsc_df=None,
    sample_info=None,
    depmap_release="unknown",
):
    """
    Full multi-modal analysis entry point.
    Returns (EvidenceMatrix, agreement_report, narrative, recommendation_map, guardrails).
    """
    from .models import MultiModalResult

    matrix = build_evidence_matrix(
        query=query,
        sl_partners=sl_partners,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
        prism_df=prism_df,
        prism_meta=prism_meta,
        gdsc_df=gdsc_df,
        sample_info=sample_info,
        depmap_release=depmap_release,
    )

    agreement_report = build_agreement_report(matrix)
    guardrails = _collect_guardrails(matrix)
    narrative = build_updated_narrative(
        gene=query.gene,
        cancer_type=query.cancer_type,
        matrix=matrix,
        guardrails=guardrails,
    )
    rec_map = matrix.recommendation_summary()

    return MultiModalResult(
        query_gene=query.gene,
        cancer_type=query.cancer_type,
        evidence_matrix=matrix,
        agreement_report=agreement_report,
        updated_narrative=narrative,
        recommendation_map=rec_map,
        guardrails_applied=guardrails,
        warnings=matrix.notes,
    )


# ── Helper functions ──────────────────────────────────────────────────────────

def _build_crispr_map(sl_partners) -> Dict[str, tuple]:
    """Build a gene → (delta, p, fdr, n_mut, n_wt, d) map from SLPartner list."""
    result = {}
    for p in sl_partners:
        result[p.gene.upper()] = (
            p.delta_dependency,
            p.p_value,
            p.fdr,
            p.n_mut,
            p.n_wt,
            p.effect_size_cohend,
        )
    return result


def _axis_to_target_genes(axis: CandidateAxis) -> List[str]:
    """Map a candidate axis to target gene symbols for CRISPR lookup."""
    return {
        CandidateAxis.CYTIDINE_ANALOGS: ["DCTD", "DCK", "CDA", "RRM1", "RRM2"],
        CandidateAxis.PARP_INHIBITORS:  ["PARP1", "PARP2", "TNKS", "TNKS2"],
        CandidateAxis.ATR_WEE1:        ["ATR", "WEE1", "CHEK1", "CDC25A"],
        CandidateAxis.WRN:             ["WRN"],
        CandidateAxis.IMMUNOTHERAPY:   [],   # No CRISPR targets for IO
        CandidateAxis.PKMYT1:          ["PKMYT1"],
        CandidateAxis.CUSTOM:          [],
    }.get(axis, [])


def _best_crispr_result(
    crispr_map: Dict[str, tuple],
    target_genes: List[str],
) -> Optional[tuple]:
    """
    Find the most significant CRISPR result among target genes.
    Returns the tuple with lowest FDR, or None if no targets found.
    """
    candidates = [crispr_map[g] for g in target_genes if g in crispr_map]
    if not candidates:
        return None
    return min(candidates, key=lambda x: x[2])  # sort by FDR (index 2)


def _collect_guardrails(matrix: EvidenceMatrix) -> List[str]:
    """
    Check all rows for guardrail triggers and return descriptions.
    These are surfaced in the narrative and API response.
    """
    guardrails = []
    for row in matrix.rows:
        cells = row.cells()
        crispr_neg = cells["crispr"].status == ModalityStatus.NEGATIVE
        pharma_missing = (
            cells["prism"].status == ModalityStatus.MISSING
            and cells["gdsc"].status == ModalityStatus.MISSING
        )
        if crispr_neg and pharma_missing:
            guardrails.append(
                f"[{row.axis_label}] CRISPR is negative but pharmacologic screens have NOT "
                "been examined. Cannot declare 'no vulnerability' — axis remains 'Mechanistic "
                "candidate only' until PRISM/GDSC data are interrogated."
            )

        # Confound flag
        for mod_name, cell in cells.items():
            if cell.is_confound_flagged:
                guardrails.append(
                    f"[{row.axis_label}] Pharmacologic signal in {mod_name} may be "
                    "MSI/MMR confounded, not MBD4-specific. Requires MBD4-LOF & MMR-intact "
                    "sub-stratification to confirm specificity."
                )

    return guardrails
