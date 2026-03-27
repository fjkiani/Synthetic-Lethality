"""
DepMap data loader — downloads and caches all required matrices.

Strategy:
  1. Check disk cache first (parquet format for speed).
  2. Download from Figshare / DepMap portal on cache miss.
  3. Expose typed DataFrames to the analysis layer.

Cell-line ID convention: DepMap uses "ACH-XXXXXX" model IDs.
"""
from __future__ import annotations

import io
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import httpx
import pandas as pd

from ..core.config import get_settings

logger = logging.getLogger(__name__)
cfg = get_settings()

# ── Helpers ───────────────────────────────────────────────────────────────────

def _cache_path(name: str) -> Path:
    return cfg.depmap_cache_dir / f"{name}.parquet"


def _load_or_download(
    name: str,
    url: str,
    index_col: Optional[str] = None,
    transpose: bool = False,
) -> pd.DataFrame:
    """Return a DataFrame, using disk cache when available."""
    cache = _cache_path(name)
    if cache.exists():
        logger.info("Cache hit: %s", name)
        return pd.read_parquet(cache)

    logger.info("Downloading %s from %s …", name, url)
    response = httpx.get(url, follow_redirects=True, timeout=300)
    response.raise_for_status()

    content = response.content
    try:
        df = pd.read_csv(io.BytesIO(content), index_col=0, low_memory=False)
    except Exception:
        df = pd.read_csv(io.BytesIO(content), low_memory=False)

    if transpose:
        df = df.T
    if index_col and index_col in df.columns:
        df = df.set_index(index_col)

    df.to_parquet(cache, index=True)
    logger.info("Saved %s → %s", name, cache)
    return df


# ── Public loaders ────────────────────────────────────────────────────────────

def load_crispr_gene_effect() -> pd.DataFrame:
    """
    Returns: DataFrame shape (n_cell_lines, n_genes).
    Rows = ACH model IDs, Columns = 'GENE (ENTREZ)' format or plain symbols.
    Chronos scores: more negative = more essential.
    """
    df = _load_or_download("crispr_gene_effect", cfg.depmap_crispr_url)
    # Normalise column names: strip ENTREZ suffix → plain gene symbols
    df.columns = [c.split(" (")[0] if " (" in c else c for c in df.columns]
    return df


def load_mutations() -> pd.DataFrame:
    """
    Returns: long-format MAF-style DataFrame with columns including
    ModelID, HugoSymbol, VariantClassification, etc.
    """
    return _load_or_download("mutations", cfg.depmap_mutation_url)


def load_cna() -> pd.DataFrame:
    """
    Returns: DataFrame shape (n_cell_lines, n_genes).
    Values = relative copy number (log2 ratio centred at 0).
    """
    df = _load_or_download("cna", cfg.depmap_cna_url)
    df.columns = [c.split(" (")[0] if " (" in c else c for c in df.columns]
    return df


def load_expression() -> pd.DataFrame:
    """
    Returns: DataFrame shape (n_cell_lines, n_genes).
    Values = log2(TPM+1).
    """
    df = _load_or_download("expression", cfg.depmap_expression_url)
    df.columns = [c.split(" (")[0] if " (" in c else c for c in df.columns]
    return df


def load_sample_info() -> pd.DataFrame:
    """
    Returns: DataFrame indexed by ModelID with lineage, subtype, etc.
    Key columns: ModelID, OncotreeLineage, OncotreePrimaryDisease,
                 OncotreeSubtype, Age, Sex, PrimaryOrMetastasis.
    """
    df = _load_or_download("sample_info", cfg.depmap_sample_info_url)
    if "ModelID" in df.columns:
        df = df.set_index("ModelID")
    return df


def load_prism_viability() -> pd.DataFrame:
    """
    Returns: DataFrame shape (n_cell_lines, n_compounds).
    Values = log2-fold-change viability (PRISM LFC).
    """
    return _load_or_download("prism_viability", cfg.depmap_prism_url)


def load_prism_meta() -> pd.DataFrame:
    """
    Returns: compound metadata DataFrame.
    Key columns: column_name, name, target, moa, phase, smiles.
    """
    return _load_or_download("prism_meta", cfg.depmap_prism_meta_url)


# ── Stratification helpers ─────────────────────────────────────────────────────

def get_mutant_wt_lines(
    gene: str,
    mutation_df: pd.DataFrame,
    sample_info: pd.DataFrame,
    cancer_type: Optional[str] = None,
    mutation_type: str = "any_mutation",
    cna_df: Optional[pd.DataFrame] = None,
    homdel_threshold: float = -1.0,
    amp_threshold: float = 1.0,
) -> Tuple[list, list]:
    """
    Split cell lines into mutant vs WT groups for a query gene.

    mutation_type options:
        'loss_of_function'      – damaging SNVs/indels (splice, nonsense, frameshift, missense-damaging)
        'gain_of_function'      – missense activating
        'homozygous_deletion'   – CNA < homdel_threshold
        'amplification'         – CNA > amp_threshold
        'any_mutation'          – any non-silent mutation

    Returns: (mutant_ids, wt_ids)
    """
    # Filter by cancer type first
    if cancer_type:
        lineage_cols = [
            c for c in sample_info.columns
            if c in ("OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype",
                     "lineage", "primary_disease", "Subtype")
        ]
        mask = pd.Series(False, index=sample_info.index)
        for col in lineage_cols:
            mask |= sample_info[col].str.contains(cancer_type, case=False, na=False)
        subset_lines = sample_info.index[mask].tolist()
    else:
        subset_lines = sample_info.index.tolist()

    if not subset_lines:
        raise ValueError(f"No cell lines found for cancer_type='{cancer_type}'")

    # Determine mutant lines
    if mutation_type in ("any_mutation", "loss_of_function", "gain_of_function"):
        # Work from mutation long-format table
        gene_col = next(
            (c for c in mutation_df.columns if c in ("HugoSymbol", "gene", "Gene_Symbol")),
            None,
        )
        model_col = next(
            (c for c in mutation_df.columns if c in ("ModelID", "DepMap_ID", "model_id")),
            None,
        )
        if gene_col is None or model_col is None:
            raise ValueError("Mutation table missing expected gene/model columns")

        gene_muts = mutation_df[mutation_df[gene_col] == gene]

        if mutation_type == "loss_of_function":
            lof_classes = {
                "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                "Splice_Site", "Splice_Region", "In_Frame_Del", "In_Frame_Ins",
                "Translation_Start_Site", "Nonstop_Mutation",
            }
            vc_col = next(
                (c for c in mutation_df.columns
                 if c in ("VariantClassification", "variant_classification", "Variant_Classification")),
                None,
            )
            if vc_col:
                gene_muts = gene_muts[gene_muts[vc_col].isin(lof_classes)]
        elif mutation_type == "gain_of_function":
            gene_muts = gene_muts[
                gene_muts.get("VariantAnnotation", pd.Series(dtype=str)).str.contains(
                    "gain.of.function|activating", case=False, na=False
                )
            ]

        mutant_ids = gene_muts[model_col].unique().tolist()

    elif mutation_type == "homozygous_deletion":
        if cna_df is None:
            raise ValueError("CNA matrix required for homozygous_deletion stratification")
        if gene not in cna_df.columns:
            raise ValueError(f"Gene {gene} not found in CNA matrix")
        mutant_ids = cna_df.index[cna_df[gene] < homdel_threshold].tolist()

    elif mutation_type == "amplification":
        if cna_df is None:
            raise ValueError("CNA matrix required for amplification stratification")
        if gene not in cna_df.columns:
            raise ValueError(f"Gene {gene} not found in CNA matrix")
        mutant_ids = cna_df.index[cna_df[gene] > amp_threshold].tolist()

    else:
        raise ValueError(f"Unknown mutation_type: {mutation_type}")

    # Intersect with cancer-type subset
    mutant_ids = [m for m in mutant_ids if m in subset_lines]
    wt_ids = [l for l in subset_lines if l not in mutant_ids]

    return mutant_ids, wt_ids


# ── Stratified Result model ───────────────────────────────────────────────────

class StratifiedResult:
    """
    Result of stratified CRISPR dependency analysis for one checkpoint gene.
    Plain dataclass (not Pydantic) to avoid circular import with multimodal models.
    """
    __slots__ = [
        "checkpoint_gene", "context_gene",
        "mean_dep_mutant", "mean_dep_wt",
        "delta_dep", "p_value", "effect_size_cohend",
        "n_mutant", "n_wt",
        "is_more_essential_in_mutant",
        "notes",
    ]

    def __init__(
        self,
        checkpoint_gene: str,
        context_gene: str,
        mean_dep_mutant: float,
        mean_dep_wt: float,
        delta_dep: float,
        p_value: float,
        effect_size_cohend: float,
        n_mutant: int,
        n_wt: int,
        is_more_essential_in_mutant: bool,
        notes: str = "",
    ):
        self.checkpoint_gene = checkpoint_gene
        self.context_gene = context_gene
        self.mean_dep_mutant = mean_dep_mutant
        self.mean_dep_wt = mean_dep_wt
        self.delta_dep = delta_dep
        self.p_value = p_value
        self.effect_size_cohend = effect_size_cohend
        self.n_mutant = n_mutant
        self.n_wt = n_wt
        self.is_more_essential_in_mutant = is_more_essential_in_mutant
        self.notes = notes

    def __repr__(self) -> str:
        return (
            f"StratifiedResult({self.checkpoint_gene} in {self.context_gene}-LOF: "
            f"Δdep={self.delta_dep:+.3f}, p={self.p_value:.4f}, "
            f"d={self.effect_size_cohend:.3f}, n_mut={self.n_mutant})"
        )


def get_checkpoint_dependency_in_context(
    context_gene: str,
    checkpoint_genes: List[str],
    mutation_type: str = "loss_of_function",
    cancer_type: Optional[str] = None,
) -> Dict[str, "StratifiedResult"]:
    """
    For each checkpoint gene, compute Chronos dependency stratified by context_gene status.

    This answers: "Are {context_gene}-LOF cells more dependent on each checkpoint gene?"

    Args:
        context_gene:     Stratification gene (e.g. "MBD4")
        checkpoint_genes: List of checkpoint/target genes (e.g. ["ATR", "CHEK1", "WEE1", "PKMYT1"])
        mutation_type:    How to define context_gene-LOF (see get_mutant_wt_lines)
        cancer_type:      Optional cancer type filter

    Returns:
        Dict mapping checkpoint_gene → StratifiedResult

    Key threshold for combinatorial tier unlock:
        delta_dep < -0.15 in context_gene-LOF vs WT → combinatorial tier unlockable
        (Chronos: more negative = more essential)

    Note:
        This function DOWNLOADS DepMap data (CRISPR + mutations + sample_info).
        Use the cached versions — first call may be slow, subsequent calls are fast.
        If DepMap data is unavailable, returns empty dict with logged warning.
    """
    results: Dict[str, StratifiedResult] = {}

    try:
        crispr_df = load_crispr_gene_effect()
        mut_df = load_mutations()
        sample_df = load_sample_info()
    except Exception as e:
        logger.warning("get_checkpoint_dependency_in_context: data load failed: %s", e)
        return results

    try:
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=context_gene,
            mutation_df=mut_df,
            sample_info=sample_df,
            cancer_type=cancer_type,
            mutation_type=mutation_type,
        )
    except Exception as e:
        logger.warning("get_checkpoint_dependency_in_context: stratification failed: %s", e)
        return results

    if len(mutant_ids) < 3 or len(wt_ids) < 3:
        logger.warning(
            "get_checkpoint_dependency_in_context: insufficient lines "
            "(n_mut=%d, n_wt=%d) for %s",
            len(mutant_ids), len(wt_ids), context_gene,
        )
        return results

    # Filter CRISPR df to lines we have
    mut_lines_in_crispr = [m for m in mutant_ids if m in crispr_df.index]
    wt_lines_in_crispr  = [w for w in wt_ids     if w in crispr_df.index]

    if len(mut_lines_in_crispr) < 3 or len(wt_lines_in_crispr) < 3:
        logger.warning(
            "get_checkpoint_dependency_in_context: insufficient CRISPR-matched lines "
            "(n_mut=%d, n_wt=%d)",
            len(mut_lines_in_crispr), len(wt_lines_in_crispr),
        )
        return results

    for chk_gene in checkpoint_genes:
        chk_gene_upper = chk_gene.upper()
        if chk_gene_upper not in crispr_df.columns:
            logger.debug("Checkpoint gene %s not in CRISPR matrix", chk_gene_upper)
            continue

        mut_deps = crispr_df.loc[mut_lines_in_crispr, chk_gene_upper].dropna().values
        wt_deps  = crispr_df.loc[wt_lines_in_crispr,  chk_gene_upper].dropna().values

        if len(mut_deps) < 3 or len(wt_deps) < 3:
            continue

        mean_mut = float(mut_deps.mean())
        mean_wt  = float(wt_deps.mean())
        delta    = mean_mut - mean_wt  # negative = more essential in mutant

        # Two-sample t-test (Welch's)
        try:
            import scipy.stats as stats
            t_stat, p_val = stats.ttest_ind(mut_deps, wt_deps, equal_var=False)
        except ImportError:
            # Fallback: no scipy — p_value undefined
            p_val = float("nan")

        # Cohen's d
        pooled_std = float(
            (
                ((len(mut_deps) - 1) * mut_deps.std() ** 2 +
                 (len(wt_deps)  - 1) * wt_deps.std()  ** 2)
                / (len(mut_deps) + len(wt_deps) - 2)
            ) ** 0.5
        )
        cohend = (mean_mut - mean_wt) / pooled_std if pooled_std > 0 else 0.0

        results[chk_gene_upper] = StratifiedResult(
            checkpoint_gene=chk_gene_upper,
            context_gene=context_gene.upper(),
            mean_dep_mutant=round(mean_mut, 4),
            mean_dep_wt=round(mean_wt, 4),
            delta_dep=round(delta, 4),
            p_value=round(float(p_val), 6) if not _is_nan(p_val) else float("nan"),
            effect_size_cohend=round(cohend, 3),
            n_mutant=len(mut_deps),
            n_wt=len(wt_deps),
            is_more_essential_in_mutant=delta < -0.15,
            notes=(
                f"Δdep threshold for combinatorial unlock: < -0.15 "
                f"(current: {delta:+.3f})"
            ),
        )

    return results


def _is_nan(v) -> bool:
    """Safe NaN check that works for both float and non-float."""
    try:
        import math
        return math.isnan(v)
    except (TypeError, ValueError):
        return False
