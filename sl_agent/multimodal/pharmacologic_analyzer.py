"""
Pharmacologic Analyzer — stratified drug-screen analysis for multi-modal integration.

Interrogates PRISM and GDSC data with three stratification axes:
  1. MBD4 LOF vs WT
  2. MSI-H vs MSS
  3. Combinations (MBD4-LOF & MMR-intact vs MSI-H without MBD4)

Principle: CRISPR is one axis. A robust, mechanistically coherent pharmacologic
signal in the absence of CRISPR dependency = "compound-level vulnerability
without genetic dependency" — not "no effect."
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .models import ModalityEvidence, ModalityStatus, Modality

logger = logging.getLogger(__name__)

# ── Drug name → canonical axis mapping ───────────────────────────────────────
# Used to bucket individual compounds into axes for the matrix.

_DRUG_TO_AXIS: Dict[str, str] = {
    # Cytidine analogs
    "gemcitabine": "cytidine_analogs",
    "cytarabine": "cytidine_analogs",
    "decitabine": "cytidine_analogs",
    "azacitidine": "cytidine_analogs",
    "5-azacytidine": "cytidine_analogs",
    "ara-c": "cytidine_analogs",
    # PARP inhibitors
    "olaparib": "parp_inhibitors",
    "niraparib": "parp_inhibitors",
    "talazoparib": "parp_inhibitors",
    "rucaparib": "parp_inhibitors",
    "veliparib": "parp_inhibitors",
    "pamiparib": "parp_inhibitors",
    "fluzoparib": "parp_inhibitors",
    # ATR/WEE1
    "vx-970": "atr_wee1",
    "vx970": "atr_wee1",
    "berzosertib": "atr_wee1",
    "elimusertib": "atr_wee1",
    "adavosertib": "atr_wee1",
    "azd1775": "atr_wee1",
    "azd-1775": "atr_wee1",
    "mk-1775": "atr_wee1",
    "ceralasertib": "atr_wee1",
    # PKMYT1
    "rp-6306":    "pkmyt1",
    "rp6306":     "pkmyt1",
    "pkmyt1i":    "pkmyt1",
    "pkmyt1 inhibitor": "pkmyt1",
    # WRN
    "vx-803": "wrn",
    "mrtx1719": "wrn",
    # Immunotherapy (IO — handled separately, not PRISM/GDSC primary)
    "pembrolizumab": "immunotherapy",
    "nivolumab": "immunotherapy",
    "atezolizumab": "immunotherapy",
}

# MSI/MMR marker columns we look for in sample_info
_MSI_COLS = ("msi_status", "microsatellite_instability_status", "mss_msi")
_MMR_COLS = ("mmr_status", "mlh1_methylated", "mmr_deficient")

MIN_N_PER_GROUP = 5   # minimum cell lines per group for any test


@dataclass
class StratifiedDrugResult:
    """Pharmacologic screen result for one drug, one stratification."""
    drug_name: str
    axis: str                         # e.g. "cytidine_analogs"
    dataset: str                      # "PRISM_repurposing", "GDSC1"
    stratifier: str                   # e.g. "MBD4_LOF_vs_WT"
    n_mut: int
    n_wt: int
    delta_response: float             # mean(mutant) - mean(WT), lower = more sensitive
    effect_size: float                # Cohen's d
    p_value: float
    fdr: float = 1.0
    is_significant: bool = False
    direction_consistent_with_mechanism: Optional[bool] = None
    notes: str = ""


def _cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled Cohen's d (a vs b)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0.0
    pooled_std = np.sqrt(((na - 1) * a.std(ddof=1) ** 2 + (nb - 1) * b.std(ddof=1) ** 2)
                         / (na + nb - 2))
    return (a.mean() - b.mean()) / pooled_std if pooled_std > 0 else 0.0


def _run_wilcoxon(mut_vals: np.ndarray, wt_vals: np.ndarray) -> Tuple[float, float]:
    """One-sided Wilcoxon rank-sum (mutant more sensitive = more negative AUC/LFC)."""
    if len(mut_vals) < MIN_N_PER_GROUP or len(wt_vals) < MIN_N_PER_GROUP:
        return 1.0, 0.0
    try:
        stat, p = stats.mannwhitneyu(mut_vals, wt_vals, alternative="less")
        return float(p), float(stat)
    except Exception:
        return 1.0, 0.0


def _get_msi_ids(sample_info: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """
    Return (msi_h_ids, mss_ids) from sample_info if any MSI column exists.
    Returns ([], []) if no MSI data available.
    """
    si = sample_info
    for col in _MSI_COLS + _MMR_COLS:
        if col in si.columns:
            col_vals = si[col].astype(str).str.upper()
            msi_h = si[col_vals.isin({"MSI-H", "MSI_H", "HIGH", "TRUE", "1", "DEFICIENT"})].index.tolist()
            mss = si[col_vals.isin({"MSS", "MSI-L", "MSI_S", "FALSE", "0", "PROFICIENT", "STABLE"})].index.tolist()
            if len(msi_h) + len(mss) >= 10:
                return msi_h, mss
    return [], []


def analyze_drug_screen(
    gene: str,
    mutant_ids: List[str],
    wt_ids: List[str],
    prism_df: Optional[pd.DataFrame],
    prism_meta: Optional[pd.DataFrame],
    gdsc_df: Optional[pd.DataFrame],
    sample_info: Optional[pd.DataFrame],
    target_axes: Optional[List[str]] = None,
    stratify_by_msi: bool = True,
) -> List[StratifiedDrugResult]:
    """
    Main entry point for pharmacologic analysis.

    Returns one StratifiedDrugResult per (drug, stratifier) combination tested.
    """
    results: List[StratifiedDrugResult] = []
    target_axes_set = set(target_axes) if target_axes else None

    # ── PRISM ─────────────────────────────────────────────────────────────────
    if prism_df is not None and prism_meta is not None:
        prism_results = _analyze_prism(
            gene=gene,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
            prism_df=prism_df,
            prism_meta=prism_meta,
            sample_info=sample_info,
            target_axes=target_axes_set,
            stratify_by_msi=stratify_by_msi,
        )
        results.extend(prism_results)

    # ── GDSC ──────────────────────────────────────────────────────────────────
    if gdsc_df is not None:
        gdsc_results = _analyze_gdsc(
            gene=gene,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
            gdsc_df=gdsc_df,
            sample_info=sample_info,
            target_axes=target_axes_set,
            stratify_by_msi=stratify_by_msi,
        )
        results.extend(gdsc_results)

    # BH FDR correction across all PRISM results together
    if results:
        pvals = [r.p_value for r in results]
        try:
            _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")
            for r, fdr in zip(results, fdrs):
                r.fdr = float(fdr)
                r.is_significant = fdr < 0.25 and abs(r.delta_response) > 0.1
        except Exception as e:
            logger.warning("FDR correction failed: %s", e)

    return results


def _analyze_prism(
    gene: str,
    mutant_ids: List[str],
    wt_ids: List[str],
    prism_df: pd.DataFrame,
    prism_meta: pd.DataFrame,
    sample_info: Optional[pd.DataFrame],
    target_axes: Optional[set],
    stratify_by_msi: bool,
) -> List[StratifiedDrugResult]:
    results = []

    # PRISM: rows = cell lines, cols = compound_ids
    # Find drug columns that map to our target axes
    meta_lower = prism_meta.copy()
    name_col = next((c for c in prism_meta.columns if "name" in c.lower()), None)
    if name_col is None:
        logger.warning("PRISM meta has no name column; skipping PRISM pharmacologic analysis")
        return results

    # Build drug_id → axis mapping
    drug_axis_map: Dict[str, str] = {}
    for idx, row in prism_meta.iterrows():
        dname = str(row.get(name_col, "")).lower().strip()
        axis = _match_drug_to_axis(dname)
        if axis and (target_axes is None or axis in target_axes):
            drug_axis_map[str(idx)] = axis

    if not drug_axis_map:
        return results

    # Restrict to relevant drug columns
    drug_cols = [c for c in prism_df.columns if str(c) in drug_axis_map]
    if not drug_cols:
        return results

    # Filter to cell lines present in PRISM
    all_prism_lines = set(prism_df.index.tolist())
    mut_in_prism = [m for m in mutant_ids if m in all_prism_lines]
    wt_in_prism  = [w for w in wt_ids      if w in all_prism_lines]

    if len(mut_in_prism) < MIN_N_PER_GROUP or len(wt_in_prism) < MIN_N_PER_GROUP:
        logger.info(
            "PRISM: insufficient lines for gene=%s (mut=%d, wt=%d)",
            gene, len(mut_in_prism), len(wt_in_prism),
        )
        return results

    # MSI stratification
    msi_h_ids, mss_ids = ([], [])
    if stratify_by_msi and sample_info is not None:
        msi_h_ids, mss_ids = _get_msi_ids(sample_info)
        msi_h_in_prism = [m for m in msi_h_ids if m in all_prism_lines]
        mss_in_prism   = [m for m in mss_ids   if m in all_prism_lines]
    else:
        msi_h_in_prism = mss_in_prism = []

    for drug_col in drug_cols:
        axis = drug_axis_map[str(drug_col)]
        drug_name_raw = str(prism_meta.loc[drug_col, name_col]) if drug_col in prism_meta.index else str(drug_col)

        drug_series = prism_df[drug_col].dropna()

        # Primary stratification: MBD4 LOF vs WT
        mut_vals = drug_series.reindex(mut_in_prism).dropna().values
        wt_vals  = drug_series.reindex(wt_in_prism).dropna().values

        if len(mut_vals) < MIN_N_PER_GROUP or len(wt_vals) < MIN_N_PER_GROUP:
            continue

        p_val, _ = _run_wilcoxon(mut_vals, wt_vals)
        delta = float(np.mean(mut_vals) - np.mean(wt_vals))
        d = _cohens_d(mut_vals, wt_vals)

        results.append(StratifiedDrugResult(
            drug_name=drug_name_raw,
            axis=axis,
            dataset="PRISM_repurposing",
            stratifier="MBD4_LOF_vs_WT",
            n_mut=len(mut_vals),
            n_wt=len(wt_vals),
            delta_response=delta,
            effect_size=d,
            p_value=p_val,
            direction_consistent_with_mechanism=(delta < 0),  # lower = more sensitive
        ))

        # MSI stratification: test if signal is MSI-driven vs MBD4-specific
        if len(msi_h_in_prism) >= MIN_N_PER_GROUP and len(mss_in_prism) >= MIN_N_PER_GROUP:
            msi_vals = drug_series.reindex(msi_h_in_prism).dropna().values
            mss_vals = drug_series.reindex(mss_in_prism).dropna().values
            if len(msi_vals) >= MIN_N_PER_GROUP and len(mss_vals) >= MIN_N_PER_GROUP:
                p_msi, _ = _run_wilcoxon(msi_vals, mss_vals)
                delta_msi = float(np.mean(msi_vals) - np.mean(mss_vals))
                d_msi = _cohens_d(msi_vals, mss_vals)
                results.append(StratifiedDrugResult(
                    drug_name=drug_name_raw,
                    axis=axis,
                    dataset="PRISM_repurposing",
                    stratifier="MSI_H_vs_MSS",
                    n_mut=len(msi_vals),
                    n_wt=len(mss_vals),
                    delta_response=delta_msi,
                    effect_size=d_msi,
                    p_value=p_msi,
                    notes="MSI confound check — compare with MBD4_LOF_vs_WT",
                ))

    return results


def _analyze_gdsc(
    gene: str,
    mutant_ids: List[str],
    wt_ids: List[str],
    gdsc_df: pd.DataFrame,
    sample_info: Optional[pd.DataFrame],
    target_axes: Optional[set],
    stratify_by_msi: bool,
) -> List[StratifiedDrugResult]:
    """Analyze GDSC AUC/IC50 data with same stratification logic as PRISM."""
    results = []
    drug_name_col = next(
        (c for c in gdsc_df.columns if c.lower() in ("drug_name", "compound", "drug")), None
    )
    cell_id_col = next(
        (c for c in gdsc_df.columns if c.lower() in ("cell_line_name", "model_id", "ccle_name")), None
    )
    response_col = next(
        (c for c in gdsc_df.columns if c.lower() in ("ln_ic50", "auc", "z_score")), None
    )
    if not all([drug_name_col, cell_id_col, response_col]):
        logger.warning("GDSC DataFrame missing required columns; skipping GDSC analysis")
        return results

    all_gdsc_ids = set(gdsc_df[cell_id_col].dropna().tolist())
    mut_in_gdsc = [m for m in mutant_ids if m in all_gdsc_ids]
    wt_in_gdsc  = [w for w in wt_ids      if w in all_gdsc_ids]

    if len(mut_in_gdsc) < MIN_N_PER_GROUP or len(wt_in_gdsc) < MIN_N_PER_GROUP:
        return results

    for drug_name, group in gdsc_df.groupby(drug_name_col):
        axis = _match_drug_to_axis(str(drug_name).lower())
        if axis is None:
            continue
        if target_axes and axis not in target_axes:
            continue

        sub = group.set_index(cell_id_col)[response_col].dropna()
        mut_vals = sub.reindex(mut_in_gdsc).dropna().values
        wt_vals  = sub.reindex(wt_in_gdsc).dropna().values

        if len(mut_vals) < MIN_N_PER_GROUP or len(wt_vals) < MIN_N_PER_GROUP:
            continue

        p_val, _ = _run_wilcoxon(mut_vals, wt_vals)
        delta = float(np.mean(mut_vals) - np.mean(wt_vals))
        d = _cohens_d(mut_vals, wt_vals)

        results.append(StratifiedDrugResult(
            drug_name=str(drug_name),
            axis=axis,
            dataset="GDSC",
            stratifier="MBD4_LOF_vs_WT",
            n_mut=len(mut_vals),
            n_wt=len(wt_vals),
            delta_response=delta,
            effect_size=d,
            p_value=p_val,
            direction_consistent_with_mechanism=(delta < 0),
        ))

    return results


def _match_drug_to_axis(drug_name_lower: str) -> Optional[str]:
    """Map a drug name string to a CandidateAxis string, or None if not recognized."""
    for kw, axis in _DRUG_TO_AXIS.items():
        if kw in drug_name_lower:
            return axis
    return None


def aggregate_by_axis(
    results: List[StratifiedDrugResult],
    stratifier: str = "MBD4_LOF_vs_WT",
) -> Dict[str, ModalityEvidence]:
    """
    Collapse individual drug results into one ModalityEvidence per axis/dataset.

    Returns dict: axis_name → ModalityEvidence for the pharmacologic modality.
    """
    # Group by axis+dataset+stratifier
    by_axis: Dict[str, List[StratifiedDrugResult]] = {}
    for r in results:
        if r.stratifier != stratifier:
            continue
        by_axis.setdefault(r.axis, []).append(r)

    aggregated: Dict[str, ModalityEvidence] = {}
    for axis, drug_results in by_axis.items():
        # Summarise: use most significant result per axis
        sig_results = [r for r in drug_results if r.is_significant]
        best = min(drug_results, key=lambda r: r.p_value)
        all_consistent = all(
            r.direction_consistent_with_mechanism is True for r in drug_results
        )
        any_significant = any(r.is_significant for r in drug_results)

        if any_significant and all_consistent:
            status = ModalityStatus.POSITIVE
        elif any_significant and not all_consistent:
            status = ModalityStatus.MIXED
        elif best.p_value < 0.25 and best.direction_consistent_with_mechanism:
            status = ModalityStatus.MIXED
        else:
            status = ModalityStatus.NEGATIVE

        # Check if MSI confound test was run
        msi_results = [r for r in results if r.axis == axis and r.stratifier == "MSI_H_vs_MSS"]
        is_confounded = False
        if msi_results:
            best_msi = min(msi_results, key=lambda r: r.p_value)
            # Confound if MSI effect is at least as large and significant
            if best_msi.is_significant and abs(best_msi.delta_response) >= abs(best.delta_response):
                is_confounded = True

        summaries = []
        for r in sig_results[:3]:  # top 3 significant
            summaries.append(
                f"{r.drug_name} ({r.dataset}): Δ={r.delta_response:+.3f}, "
                f"d={r.effect_size:.2f}, p={r.p_value:.3f}"
            )

        modality = (
            Modality.PRISM_PHARMACOLOGIC
            if "PRISM" in (drug_results[0].dataset if drug_results else "")
            else Modality.GDSC_PHARMACOLOGIC
        )

        aggregated[axis] = ModalityEvidence(
            modality=modality,
            status=ModalityStatus.CONFOUNDED if is_confounded else status,
            delta_auc=best.delta_response,
            effect_size=best.effect_size,
            p_value=best.p_value,
            fdr=best.fdr,
            n_mut=best.n_mut,
            n_wt=best.n_wt,
            drug_screen_dataset=best.dataset,
            stratifier=stratifier,
            summary=(
                f"{len(sig_results)}/{len(drug_results)} drugs significant. "
                + "; ".join(summaries)
            ) if summaries else f"No significant hits. Best: {best.drug_name} p={best.p_value:.3f}",
            is_confound_flagged=is_confounded,
        )

    return aggregated
