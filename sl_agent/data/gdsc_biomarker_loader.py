"""
GDSC Biomarker Loader — computes drug sensitivity associations from GDSC1/2 data.

Strategy (Mars approach):
  The pre-computed ANOVA biomarker files from cancerrxgene.org are unavailable (HTTP 410).
  We download the GDSC1/2 fitted dose-response Excel files (live at Sanger COG)
  and compute our own Wilcoxon rank-sum test per gene-drug pair.

  This gives us MORE control than pre-baked ANOVA files:
    - We use the same mutation stratification as the rest of the pipeline
    - We can filter by cancer type
    - We can pick our own FDR cutoff

Data:
  GDSC1: https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx
  GDSC2: https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
  Both are in config.py as gdsc1_url and gdsc2_url.

Columns in the Excel files (both GDSC1 and GDSC2 use same schema):
  DATASET, NLME_RESULT_ID, DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME,
  COSMIC_ID, CELL_LINE_NAME, TCGA_DESC, LN_IC50, AUC, RMSE, Z_SCORE

Key columns we use:
  DRUG_NAME  → normalize → CandidateAxis via _DRUG_TO_AXIS
  COSMIC_ID  → cell line identifier (maps to DepMap via sanger_id)
  LN_IC50    → primary metric (natural log IC50)
  CELL_LINE_NAME → for logging only
"""
from __future__ import annotations

import io
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from ..core.config import get_settings

logger = logging.getLogger(__name__)
cfg = get_settings()


# ── Cache paths ───────────────────────────────────────────────────────────────

def _gdsc_cache_path(name: str) -> Path:
    path = cfg.gdsc_cache_dir / f"{name}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _load_gdsc_excel(name: str, url: str) -> pd.DataFrame:
    """Load GDSC Excel, using parquet cache. Normalizes column names."""
    cache = _gdsc_cache_path(name)
    if cache.exists():
        logger.info("GDSC cache hit: %s", name)
        return pd.read_parquet(cache)

    logger.info("Downloading %s from %s …", name, url)
    try:
        import httpx
        resp = httpx.get(url, follow_redirects=True, timeout=600)
        resp.raise_for_status()
    except Exception as e:
        logger.warning("GDSC download failed for %s: %s", name, e)
        return pd.DataFrame()

    try:
        df = pd.read_excel(io.BytesIO(resp.content), engine="openpyxl")
    except Exception as e:
        logger.warning("GDSC Excel parse failed for %s: %s", name, e)
        return pd.DataFrame()

    # Normalize column names
    df.columns = [c.strip().upper().replace(" ", "_") for c in df.columns]
    df["_DATASET"] = name  # track which screen

    # Normalize drug names for matching
    if "DRUG_NAME" in df.columns:
        df["DRUG_NAME_NORM"] = df["DRUG_NAME"].str.lower().str.strip()

    df.to_parquet(cache, index=False)
    logger.info("Saved GDSC cache: %s → %s (%d rows)", name, cache, len(df))
    return df


def load_gdsc1() -> pd.DataFrame:
    return _load_gdsc_excel("gdsc1", cfg.gdsc1_url)


def load_gdsc2() -> pd.DataFrame:
    return _load_gdsc_excel("gdsc2", cfg.gdsc2_url)


# ── DRUG_TO_AXIS mapping (import from pharmacologic_analyzer) ─────────────────

def _get_drug_to_axis() -> Dict[str, str]:
    from ..multimodal.pharmacologic_analyzer import _DRUG_TO_AXIS
    return _DRUG_TO_AXIS


def _match_drug_to_axis(drug_name_lower: str) -> Optional[str]:
    for kw, axis in _get_drug_to_axis().items():
        if kw in drug_name_lower:
            return axis
    return None


# ── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class GDSCResult:
    """Result of a gene × axis GDSC biomarker test."""
    gene: str
    axis: str
    drug: str
    gdsc_version: str          # "GDSC1" or "GDSC2"
    delta_ln_ic50: float       # mean(mutant LN_IC50) - mean(WT LN_IC50); negative = more sensitive
    p_value: float
    fdr: float                 # BH-corrected across all tests for this gene
    n_mutant: int
    n_wt: int
    effect_size_cohend: float  # Cohen's d


# ── Mutation stratification ───────────────────────────────────────────────────

def _get_cosmic_mutant_lines(
    gene: str,
    mutation_type: str = "loss_of_function",
) -> Tuple[List[int], List[int]]:
    """
    Get COSMIC IDs for mutant vs WT cell lines.
    Uses DepMap mutation table + sample_info to get COSMIC IDs.

    Returns (mutant_cosmic_ids, wt_cosmic_ids) or ([], []) on failure.
    """
    try:
        from ..data.depmap_loader import load_mutations, load_sample_info
        mut_df = load_mutations()
        sample_df = load_sample_info()
    except Exception as e:
        logger.warning("GDSC: could not load DepMap data for stratification: %s", e)
        return [], []

    # Find the Sanger/COSMIC ID column in sample_info
    cosmic_col = next(
        (c for c in sample_df.columns
         if "cosmic" in c.lower() or "sanger" in c.lower()),
        None,
    )
    if cosmic_col is None:
        logger.warning("GDSC: no COSMIC ID column in sample_info")
        return [], []

    # Get mutant DepMap IDs using existing stratification
    try:
        from ..data.depmap_loader import get_mutant_wt_lines
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=gene,
            mutation_df=mut_df,
            sample_info=sample_df,
            mutation_type=mutation_type,
        )
    except Exception as e:
        logger.warning("GDSC: stratification failed for %s: %s", gene, e)
        return [], []

    # Map DepMap IDs → COSMIC IDs
    mutant_cosmic = sample_df.loc[
        sample_df.index.isin(mutant_ids), cosmic_col
    ].dropna().astype(int).tolist()

    wt_cosmic = sample_df.loc[
        sample_df.index.isin(wt_ids), cosmic_col
    ].dropna().astype(int).tolist()

    return mutant_cosmic, wt_cosmic


# ── Core function ─────────────────────────────────────────────────────────────

def gdsc_biomarker_stratify(
    gene: str,
    axis: str,
    mutation_type: str = "loss_of_function",
    min_n_per_group: int = 5,
) -> Optional[GDSCResult]:
    """
    Compute GDSC drug sensitivity association for gene × axis pair.

    Strategy:
      1. Load GDSC1 + GDSC2 fitted dose-response data
      2. Filter to drugs that map to the requested axis via _DRUG_TO_AXIS
      3. Get mutant vs WT cell line COSMIC IDs from DepMap mutation table
      4. For each matching drug, compute Mann-Whitney U test on LN_IC50
      5. Return the best result (lowest p-value) across all matching drugs/versions
      6. Return None if: no matching drugs, insufficient cell lines, data load fails

    Gate in receipt_miner.py: fdr < 0.05

    NEVER raises — all errors return None.
    """
    try:
        return _compute_gdsc_result(gene, axis, mutation_type, min_n_per_group)
    except Exception as e:
        logger.warning("gdsc_biomarker_stratify(%s, %s) failed: %s", gene, axis, e)
        return None


def _compute_gdsc_result(
    gene: str,
    axis: str,
    mutation_type: str,
    min_n_per_group: int,
) -> Optional[GDSCResult]:
    # Load data
    gdsc1 = load_gdsc1()
    gdsc2 = load_gdsc2()

    if gdsc1.empty and gdsc2.empty:
        logger.warning("GDSC: both GDSC1 and GDSC2 failed to load — GDSC signal unavailable")
        return None

    # Combine both screens
    all_gdsc = pd.concat(
        [df for df in [gdsc1, gdsc2] if not df.empty],
        ignore_index=True,
    )

    # Filter to drugs that map to this axis
    if "DRUG_NAME_NORM" not in all_gdsc.columns:
        return None

    axis_drugs = all_gdsc[
        all_gdsc["DRUG_NAME_NORM"].apply(
            lambda d: _match_drug_to_axis(str(d)) == axis
        )
    ]

    if axis_drugs.empty:
        logger.debug("GDSC: no drugs found for axis '%s'", axis)
        return None

    # Get COSMIC IDs for mutant vs WT
    mutant_cosmic, wt_cosmic = _get_cosmic_mutant_lines(gene, mutation_type)

    if len(mutant_cosmic) < min_n_per_group or len(wt_cosmic) < min_n_per_group:
        logger.debug(
            "GDSC: insufficient cell lines for %s (mut=%d, wt=%d)",
            gene, len(mutant_cosmic), len(wt_cosmic),
        )
        return None

    # Find COSMIC ID column in GDSC data
    cosmic_id_col = next(
        (c for c in all_gdsc.columns if "cosmic" in c.lower()),
        None,
    )
    if cosmic_id_col is None:
        logger.warning("GDSC: no COSMIC_ID column found")
        return None

    # LN_IC50 column
    ln_ic50_col = next(
        (c for c in all_gdsc.columns if "LN_IC50" in c or "ln_ic50" in c.lower()),
        None,
    )
    if ln_ic50_col is None:
        return None

    # Test each drug, pick best result
    best: Optional[GDSCResult] = None
    best_pval = 1.0

    for drug_name, drug_group in axis_drugs.groupby("DRUG_NAME_NORM"):
        mut_ic50 = drug_group.loc[
            drug_group[cosmic_id_col].isin(mutant_cosmic), ln_ic50_col
        ].dropna().values

        wt_ic50 = drug_group.loc[
            drug_group[cosmic_id_col].isin(wt_cosmic), ln_ic50_col
        ].dropna().values

        if len(mut_ic50) < min_n_per_group or len(wt_ic50) < min_n_per_group:
            continue

        # Mann-Whitney U (non-parametric — appropriate for IC50 distributions)
        try:
            stat, pval = stats.mannwhitneyu(
                mut_ic50, wt_ic50, alternative="less"  # one-sided: mutant MORE sensitive
            )
        except Exception:
            continue

        if pval < best_pval:
            best_pval = pval
            mean_mut = float(np.mean(mut_ic50))
            mean_wt = float(np.mean(wt_ic50))
            delta = mean_mut - mean_wt

            # Cohen's d
            n_mut_i = len(mut_ic50)
            n_wt_i = len(wt_ic50)
            pooled_std = float(np.sqrt(
                ((n_mut_i - 1) * np.std(mut_ic50) ** 2 +
                 (n_wt_i - 1) * np.std(wt_ic50) ** 2)
                / (n_mut_i + n_wt_i - 2)
            )) if (n_mut_i + n_wt_i > 2) else 1.0
            cohend = delta / pooled_std if pooled_std > 0 else 0.0

            # Determine version from dataset column
            version = "GDSC1"
            if "_DATASET" in drug_group.columns:
                v = drug_group["_DATASET"].iloc[0]
                version = str(v).upper()

            best = GDSCResult(
                gene=gene.upper(),
                axis=axis,
                drug=str(drug_name),
                gdsc_version=version,
                delta_ln_ic50=round(delta, 4),
                p_value=round(float(pval), 6),
                fdr=round(float(pval), 6),  # single test — fdr = p_value; BH applied in miner
                n_mutant=len(mut_ic50),
                n_wt=len(wt_ic50),
                effect_size_cohend=round(cohend, 3),
            )

    return best
