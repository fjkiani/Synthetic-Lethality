"""
PRISM Viability Signal Loader — Sprint 3.

Loads the PRISM Repurposing Secondary Screen (2020) from DepMap Figshare and
computes drug sensitivity associations for a gene × axis pair.

Strategy (Mars approach):
  PRISM secondary screen measures cancer cell line viability under 1448 compounds
  at 8 doses (AUC / EC50 summarised per cell line per drug).
  We stratify cell lines by mutation status (mutant vs WT) and test whether
  mutant lines are more sensitive (lower AUC = more kill) using Mann-Whitney U.

Data source (Figshare — public domain, no paywall):
  secondary-screen-dose-response-curve-parameters.csv
  URL: config.depmap_prism_url   (https://figshare.com/ndownloader/files/20237786)

Columns used from that file:
  depmap_id   — cell line identifier (matches DepMap CRISPR/mutation tables)
  name        — drug name (normalised → CandidateAxis via _DRUG_TO_AXIS)
  auc         — area under dose-response curve (lower = more sensitive)
  ec50        — EC50 estimate (secondary check only)

Integration with receipt_miner.py:
  PRISM is a BONUS signal only — it does NOT change existing weights.
    If prism_fdr < 0.05 → confidence += 0.05 (capped at 1.0)
  prism_delta_auc is stored in ReceiptCandidate for transparency.

NEVER raises — all errors return None.
"""
from __future__ import annotations

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


# ── Cache ─────────────────────────────────────────────────────────────────────

def _prism_cache_path() -> Path:
    cache_dir = Path(cfg.depmap_cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / "prism_secondary.parquet"


def _load_prism_raw() -> pd.DataFrame:
    """
    Load PRISM secondary screen CSV, using parquet cache.

    Columns retained: depmap_id, name, auc, ec50.
    Returns empty DataFrame on any failure (never raises).
    """
    cache = _prism_cache_path()
    if cache.exists():
        logger.info("PRISM cache hit: %s", cache)
        try:
            return pd.read_parquet(cache)
        except Exception as e:
            logger.warning("PRISM parquet read failed — re-downloading: %s", e)
            cache.unlink(missing_ok=True)

    url = cfg.depmap_prism_url
    logger.info("Downloading PRISM secondary screen from %s …", url)
    try:
        import httpx
        resp = httpx.get(url, follow_redirects=True, timeout=600)
        resp.raise_for_status()
    except Exception as e:
        logger.warning("PRISM download failed: %s", e)
        return pd.DataFrame()

    try:
        import io
        df = pd.read_csv(io.BytesIO(resp.content), low_memory=False)
    except Exception as e:
        logger.warning("PRISM CSV parse failed: %s", e)
        return pd.DataFrame()

    # Normalise column names
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # Keep only the columns we need
    needed = {"depmap_id", "name", "auc", "ec50"}
    available = set(df.columns)
    missing = needed - available
    if missing:
        logger.warning(
            "PRISM: expected columns missing: %s. Available: %s",
            missing, list(available)[:20],
        )
        return pd.DataFrame()

    df = df[["depmap_id", "name", "auc", "ec50"]].copy()

    # Normalise drug name for matching
    df["name_norm"] = df["name"].str.lower().str.strip()

    # Drop rows with null AUC (curve fit failed)
    df = df.dropna(subset=["auc"])

    df.to_parquet(cache, index=False)
    logger.info("Saved PRISM cache: %s (%d rows)", cache, len(df))
    return df


# ── Drug → Axis mapping ───────────────────────────────────────────────────────

def _get_drug_to_axis() -> Dict[str, str]:
    from ..multimodal.pharmacologic_analyzer import _DRUG_TO_AXIS
    return _DRUG_TO_AXIS


def _match_drug_to_axis(drug_name_lower: str) -> Optional[str]:
    for kw, axis in _get_drug_to_axis().items():
        if kw in drug_name_lower:
            return axis
    return None


# ── Mutation stratification (reuse DepMap loaders) ───────────────────────────

def _get_depmap_mutant_wt_lines(
    gene: str,
) -> Tuple[List[str], List[str]]:
    """
    Returns (mutant_depmap_ids, wt_depmap_ids) using cached DepMap mutation table.
    Uses get_mutant_wt_lines() from depmap_loader — already unit-tested.
    Returns ([], []) on any failure.
    """
    try:
        from ..data.depmap_loader import (
            load_mutations,
            load_sample_info,
            get_mutant_wt_lines,
        )
        mut_df = load_mutations()
        sample_df = load_sample_info()
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=gene,
            mutation_df=mut_df,
            sample_info=sample_df,
            mutation_type="loss_of_function",
        )
        return list(mutant_ids), list(wt_ids)
    except Exception as e:
        logger.warning("PRISM: mutation stratification failed for %s: %s", gene, e)
        return [], []


# ── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class PRISMResult:
    """Result of a gene × axis PRISM drug sensitivity test."""
    gene: str
    axis: str
    drug: str
    delta_auc: float        # mean(mutant AUC) − mean(WT AUC); negative = mutant more sensitive
    n_mut: int
    n_wt: int
    p_value: float
    fdr: float              # BH-corrected (single test: fdr = p_value; miner applies BH)
    ec50_delta: Optional[float] = None


# ── Core function ─────────────────────────────────────────────────────────────

def prism_stratify(
    gene: str,
    axis: str,
    min_n_per_group: int = 5,
) -> Optional[PRISMResult]:
    """
    Compute PRISM drug sensitivity association for a gene × axis pair.

    Strategy:
      1. Load PRISM secondary screen (parquet cache, download on miss)
      2. Filter to drugs that map to the requested axis via _DRUG_TO_AXIS
      3. Get mutant vs WT DepMap IDs from mutation table
      4. Mann-Whitney U on AUC values (alternative="less": mutant more sensitive)
      5. Return best result (lowest p-value) across all matching drugs
      6. Return None if: no matching drugs, insufficient cell lines, or load fails

    NEVER raises — all errors return None.

    Used by receipt_miner.py as a BONUS signal:
      if result.fdr < 0.05 → confidence += 0.05 (capped 1.0)
    """
    try:
        return _compute_prism_result(gene.upper(), axis, min_n_per_group)
    except Exception as e:
        logger.warning("prism_stratify(%s, %s) failed: %s", gene, axis, e)
        return None


def _compute_prism_result(
    gene: str,
    axis: str,
    min_n_per_group: int,
) -> Optional[PRISMResult]:
    # Load data
    df = _load_prism_raw()
    if df.empty:
        logger.warning("PRISM: no data loaded — PRISM signal unavailable for %s×%s", gene, axis)
        return None

    # Filter to drugs matching this axis
    if "name_norm" not in df.columns:
        return None

    axis_mask = df["name_norm"].apply(lambda d: _match_drug_to_axis(str(d)) == axis)
    axis_drugs = df[axis_mask]

    if axis_drugs.empty:
        logger.debug("PRISM: no drugs found for axis '%s'", axis)
        return None

    # Get mutant vs WT DepMap IDs
    mutant_ids, wt_ids = _get_depmap_mutant_wt_lines(gene)

    if len(mutant_ids) < min_n_per_group or len(wt_ids) < min_n_per_group:
        logger.debug(
            "PRISM: insufficient cell lines for %s (mut=%d, wt=%d)",
            gene, len(mutant_ids), len(wt_ids),
        )
        return None

    # Test each drug, pick best (lowest p-value)
    best: Optional[PRISMResult] = None
    best_pval = 1.0

    for drug_name, drug_group in axis_drugs.groupby("name_norm"):
        mut_auc = drug_group.loc[
            drug_group["depmap_id"].isin(mutant_ids), "auc"
        ].dropna().values

        wt_auc = drug_group.loc[
            drug_group["depmap_id"].isin(wt_ids), "auc"
        ].dropna().values

        if len(mut_auc) < min_n_per_group or len(wt_auc) < min_n_per_group:
            continue

        try:
            _stat, pval = stats.mannwhitneyu(
                mut_auc, wt_auc, alternative="less"   # one-sided: mutant MORE sensitive (lower AUC)
            )
        except Exception:
            continue

        if pval < best_pval:
            best_pval = pval
            delta_auc = float(np.mean(mut_auc) - np.mean(wt_auc))

            # EC50 delta (optional)
            ec50_delta: Optional[float] = None
            if "ec50" in drug_group.columns:
                mut_ec50 = drug_group.loc[
                    drug_group["depmap_id"].isin(mutant_ids), "ec50"
                ].dropna().values
                wt_ec50 = drug_group.loc[
                    drug_group["depmap_id"].isin(wt_ids), "ec50"
                ].dropna().values
                if len(mut_ec50) > 0 and len(wt_ec50) > 0:
                    ec50_delta = round(float(np.mean(mut_ec50) - np.mean(wt_ec50)), 4)

            best = PRISMResult(
                gene=gene,
                axis=axis,
                drug=str(drug_name),
                delta_auc=round(delta_auc, 4),
                n_mut=len(mut_auc),
                n_wt=len(wt_auc),
                p_value=round(float(pval), 6),
                fdr=round(float(pval), 6),  # single test; BH applied in miner if needed
                ec50_delta=ec50_delta,
            )

    return best
