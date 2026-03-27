"""
Orchestrator — wires all layers together for a single SL query.

Flow:
  1. Load (or reuse cached) DepMap matrices
  2. Stratify cell lines (mutant vs WT) for query gene
  3. Run SL engine → ranked SL partners
  4. Annotate pathways + framework cross-validation
  5. For each top partner, run drug mapper
  6. Assemble and return SLMapResult
"""
from __future__ import annotations

import logging
import time
from functools import lru_cache
from typing import Optional

import pandas as pd

from ..data.depmap_loader import (
    get_mutant_wt_lines,
    load_cna,
    load_crispr_gene_effect,
    load_expression,
    load_mutations,
    load_prism_meta,
    load_prism_viability,
    load_sample_info,
)
from .config import get_settings
from .drug_mapper import map_gene_to_drugs
from .models import SLMapResult, SLQueryInput
from .pathway_annotator import (
    annotate_frameworks,
    annotate_pathways,
    build_cross_validation,
)
from .sl_engine import SLEngine

logger = logging.getLogger(__name__)
cfg = get_settings()


# ── Singleton data store (loaded once, reused per process) ────────────────────

class DataStore:
    """Lazily loads and holds all DepMap matrices in memory."""

    _crispr: Optional[pd.DataFrame] = None
    _mutations: Optional[pd.DataFrame] = None
    _cna: Optional[pd.DataFrame] = None
    _expression: Optional[pd.DataFrame] = None
    _sample_info: Optional[pd.DataFrame] = None
    _prism: Optional[pd.DataFrame] = None
    _prism_meta: Optional[pd.DataFrame] = None

    @classmethod
    def ensure_loaded(cls, require_prism: bool = True) -> None:
        if cls._crispr is None:
            logger.info("Loading DepMap CRISPR matrix …")
            cls._crispr = load_crispr_gene_effect()
        if cls._mutations is None:
            logger.info("Loading mutation calls …")
            try:
                cls._mutations = load_mutations()
            except Exception as e:
                logger.warning("Mutation load failed: %s", e)
        if cls._cna is None:
            logger.info("Loading CNA matrix …")
            try:
                cls._cna = load_cna()
            except Exception as e:
                logger.warning("CNA load failed: %s", e)
        if cls._sample_info is None:
            logger.info("Loading sample info …")
            cls._sample_info = load_sample_info()
        if require_prism:
            if cls._prism is None:
                logger.info("Loading PRISM viability matrix …")
                try:
                    cls._prism = load_prism_viability()
                except Exception as e:
                    logger.warning("PRISM load failed (drug mapping may be limited): %s", e)
            if cls._prism_meta is None:
                try:
                    cls._prism_meta = load_prism_meta()
                except Exception as e:
                    logger.warning("PRISM meta load failed: %s", e)

    @classmethod
    def crispr(cls) -> pd.DataFrame:
        cls.ensure_loaded()
        return cls._crispr  # type: ignore

    @classmethod
    def mutations(cls) -> Optional[pd.DataFrame]:
        cls.ensure_loaded()
        return cls._mutations

    @classmethod
    def cna(cls) -> Optional[pd.DataFrame]:
        cls.ensure_loaded()
        return cls._cna

    @classmethod
    def sample_info(cls) -> pd.DataFrame:
        cls.ensure_loaded()
        return cls._sample_info  # type: ignore

    @classmethod
    def prism(cls) -> Optional[pd.DataFrame]:
        cls.ensure_loaded()
        return cls._prism

    @classmethod
    def prism_meta(cls) -> Optional[pd.DataFrame]:
        cls.ensure_loaded()
        return cls._prism_meta


# ── Main orchestration function ───────────────────────────────────────────────

def run_sl_analysis(query: SLQueryInput) -> SLMapResult:
    """
    Execute the full SL pipeline for the given query.
    This is the single entry point consumed by the API layer.
    """
    t0 = time.perf_counter()

    # 1. Ensure data is loaded
    DataStore.ensure_loaded()

    crispr = DataStore.crispr()
    sample_info = DataStore.sample_info()
    mutations = DataStore.mutations()
    cna = DataStore.cna()
    prism = DataStore.prism()
    prism_meta = DataStore.prism_meta()

    # 2. Stratify cell lines
    if mutations is None and query.mutation_type.value not in (
        "homozygous_deletion", "amplification"
    ):
        raise RuntimeError("Mutation table not available — cannot stratify cell lines.")

    mutant_ids, wt_ids = get_mutant_wt_lines(
        gene=query.gene,
        mutation_df=mutations if mutations is not None else pd.DataFrame(),
        sample_info=sample_info,
        cancer_type=query.cancer_type,
        mutation_type=query.mutation_type.value,
        cna_df=cna,
    )

    # If cancer-specific has too few lines, get pan-cancer lists as fallback
    MIN_LINES = cfg.min_cell_lines_cancer_specific
    if len(mutant_ids) + len(wt_ids) < MIN_LINES and query.cancer_type:
        logger.warning(
            "Only %d lines found for %s — fetching pan-cancer fallback lists",
            len(mutant_ids) + len(wt_ids),
            query.cancer_type,
        )
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=query.gene,
            mutation_df=mutations if mutations is not None else pd.DataFrame(),
            sample_info=sample_info,
            cancer_type=None,   # pan-cancer
            mutation_type=query.mutation_type.value,
            cna_df=cna,
        )

    # 3. SL engine
    engine = SLEngine(
        crispr_df=crispr,
        sample_info_df=sample_info,
        mutation_df=mutations,
        cna_df=cna,
        depmap_release=cfg.depmap_release,
    )
    input_ctx, partners, warnings_str = engine.compute_sl_partners(
        query=query,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )

    # 4. Pathway + framework annotation
    if query.include_pathway_context:
        partners = annotate_pathways(partners)
    partners = annotate_frameworks(query.gene, partners)

    cross_val = build_cross_validation(query.gene, partners)

    # 5. Drug mapping (top N partners only)
    top_partners_for_drug = partners[: query.top_n_partners]
    all_drug_pairs = []

    for partner in top_partners_for_drug:
        drug_pairs = map_gene_to_drugs(
            gene=partner.gene,
            sl_delta=partner.delta_dependency,
            sl_fdr=partner.fdr,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
            prism_df=prism,
            prism_meta=prism_meta,
            # GDSC note: gdsc_drugs_for_gene() exists in drug_mapper.py and
            # accepts a pre-loaded GDSC DataFrame (load via depmap_loader).
            # It is intentionally *not* auto-loaded here because the GDSC
            # Excel file is ~150 MB and most deployments prefer the lighter
            # PRISM-only path.  To enable GDSC enrichment:
            #   1. Call load_gdsc_viability() from depmap_loader.py
            #   2. Store it in DataStore (add _gdsc / gdsc() class methods)
            #   3. Pass it here instead of None
            gdsc_meta=None,
            top_n_drugs=query.top_n_drugs,
        )
        all_drug_pairs.extend(drug_pairs)

    # Re-sort all gene–drug pairs globally by rank_score
    all_drug_pairs.sort(key=lambda x: -x.rank_score)

    # 6. Assemble result
    warning_list = [w for w in warnings_str.split("\n") if w.strip()]

    result = SLMapResult(
        input_context=input_ctx,
        sl_partners=partners,
        gene_drug_pairs=all_drug_pairs,
        cross_validation=cross_val,
        notes=(
            f"Analysis completed in {(time.perf_counter() - t0)*1000:.0f} ms. "
            f"DepMap release: {cfg.depmap_release}."
        ),
        warnings=warning_list,
    )

    return result
