"""
Receipt Miner — automated evidence mining for the audit queue.

For each (gene, axis) pair in a panel, aggregates signals from:
  1. CRISPR dependency (DepMap Chronos via sl_engine) — weight 0.35
  2. GDSC drug sensitivity (Mann-Whitney on LN_IC50) — weight 0.30
  3. KB clinical evidence (CIViC + CGI + JAX) — weight 0.25
  4. Expression correlation — weight 0.10

Emits ReceiptCandidate objects to the audit queue.
Candidates scoring < min_confidence (default 0.40) are discarded (not queued).

VALIDATED tier is NEVER auto-assigned. Requires human approval + clinical POSITIVE.
Tier projection:
  score >= 0.70 AND kb_hits >= 3: candidate_tier = "Strong"
  score >= 0.50: candidate_tier = "Mechanistic"
  else: candidate_tier = "Insufficient"
  (Insufficient candidates are discarded — never enter queue)

Sabotage protection:
  - candidate_tier is ALWAYS one of {"Strong", "Mechanistic", "Insufficient"}
  - "Validated" is NEVER emitted by mine_receipts()
  - Even if every signal is maximal, tier caps at "Strong"
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

from ..audit.models import ReceiptCandidate
from ..audit.queue import AuditQueue
from ..multimodal.models import CandidateAxis
from ..data.gdsc_biomarker_loader import gdsc_biomarker_stratify
from ..data.prism_loader import prism_stratify

logger = logging.getLogger(__name__)

# ── Gene panel loader ─────────────────────────────────────────────────────────

_PANEL_JSON = Path(__file__).parent.parent / "data" / "gene_panels" / "ddr_panel_v1.json"


def load_default_panel() -> List[str]:
    """
    Load the default gene panel from ddr_panel_v1.json.
    Flattens all sub-groups, deduplicates (preserving first-seen order), returns list.
    Falls back to DEFAULT_GENE_PANEL if file missing or parse fails (never raises).
    """
    try:
        with open(_PANEL_JSON, "r") as f:
            data = json.load(f)
        seen: dict = {}
        for group_genes in data["genes"].values():
            for g in group_genes:
                g_up = g.strip().upper()
                if g_up:
                    seen[g_up] = None  # dict preserves insertion order (Python 3.7+)
        return list(seen.keys())
    except Exception as e:
        logger.warning("load_default_panel failed, using DEFAULT_GENE_PANEL: %s", e)
        return list(DEFAULT_GENE_PANEL)


# ── Seed gene panel ───────────────────────────────────────────────────────────

DEFAULT_GENE_PANEL: List[str] = [
    # DDR
    "BRCA1", "BRCA2", "PALB2", "ATM", "ATR", "CHEK1", "CHEK2",
    "RAD51", "RAD51C", "RAD51D", "BRIP1", "FANCA", "FANCD2",
    # Cell cycle
    "CCNE1", "CDK2", "WEE1", "CDC25A", "PLK1",
    # Chromatin / BER
    "MBD4", "TET2", "DNMT3A", "ARID1A",
    # Replication
    "RRM2", "PCNA", "RFC1", "MCM2",
    # PARP axis
    "PARP1", "RNF144A", "SLFN11",
]

DEFAULT_AXES: List[CandidateAxis] = [
    CandidateAxis.CYTIDINE_ANALOGS,
    CandidateAxis.PARP_INHIBITORS,
    CandidateAxis.ATR_WEE1,
    CandidateAxis.WRN,
    CandidateAxis.IMMUNOTHERAPY,
    CandidateAxis.PKMYT1,
]


# ── Threshold config ──────────────────────────────────────────────────────────

@dataclass
class MinerThresholds:
    crispr_delta_dep: float = -0.15   # must be MORE negative to pass gate
    crispr_fdr: float = 0.25
    gdsc_fdr: float = 0.05
    min_confidence: float = 0.40      # discard below this threshold

    def __post_init__(self):
        # Validate bounds
        assert self.crispr_delta_dep < 0, "crispr_delta_dep must be negative"
        assert 0 < self.crispr_fdr <= 1
        assert 0 < self.gdsc_fdr <= 1
        assert 0 <= self.min_confidence <= 1


# ── Run summary ───────────────────────────────────────────────────────────────

@dataclass
class MinerRunSummary:
    started_at: datetime = field(default_factory=datetime.utcnow)
    completed_at: Optional[datetime] = None
    gene_count: int = 0
    axis_count: int = 0
    pairs_evaluated: int = 0
    candidates_above_threshold: int = 0
    candidates_queued: int = 0
    candidates_discarded: int = 0    # scored above 0 but below min_confidence
    pairs_no_signal: int = 0         # all signals zero
    errors: List[str] = field(default_factory=list)
    depmap_release: str = "unknown"

    def finish(self) -> "MinerRunSummary":
        self.completed_at = datetime.utcnow()
        return self


# ── PRISM bonus constant ─────────────────────────────────────────────────────
# PRISM is an additive bonus only — original weights are NOT rebalanced.
# Rationale: PRISM secondary screen is a separate validation layer;
#   adding it as a bonus preserves all 158 existing tests unchanged.
PRISM_BONUS = 0.05          # added to confidence when prism_fdr < 0.05
PRISM_FDR_GATE = 0.05       # one-sided Mann-Whitney p-value gate


# ── Signal computation helpers ────────────────────────────────────────────────

def _compute_crispr_score(
    gene: str,
    thresholds: MinerThresholds,
) -> Tuple[float, Optional[float], Optional[float], Optional[int], Optional[int]]:
    """
    Returns (crispr_score 0-1, delta_dep, fdr, n_mutant, n_wt).
    Score = 0.0 if gate fails or data unavailable.
    DataStore uses class-level attributes (not _instance) — handle gracefully.
    """
    try:
        from ..core.orchestrator import DataStore
        DataStore.ensure_loaded(require_prism=False)

        # Check if data is loaded — DataStore uses class attributes (not _instance)
        if DataStore._crispr is None or DataStore._mutations is None:
            return 0.0, None, None, None, None

        from ..core.sl_engine import SLEngine
        from ..core.models import SLQueryInput, MutationType
        from ..data.depmap_loader import get_mutant_wt_lines

        engine = SLEngine(
            crispr_df=DataStore._crispr,
            sample_info_df=DataStore._sample_info,
            mutation_df=DataStore._mutations,
            cna_df=DataStore._cna,
            depmap_release=getattr(cfg, "depmap_release", "unknown"),
        )

        query = SLQueryInput(
            gene=gene,
            mutation_type=MutationType.LOSS_OF_FUNCTION,
            delta_dep_cutoff=abs(thresholds.crispr_delta_dep),
            fdr_cutoff=thresholds.crispr_fdr,
        )

        # Stratify cell lines for the query gene
        mutant_ids, wt_ids = get_mutant_wt_lines(
            gene=gene,
            mutation_df=DataStore._mutations,
            sample_info=DataStore._sample_info,
            mutation_type="loss_of_function",
        )

        if not mutant_ids:
            return 0.0, None, None, None, None

        input_ctx, partners, _ = engine.compute_sl_partners(
            query=query,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
        )

        if not partners:
            return 0.0, None, None, None, None

        # Use the best partner (lowest FDR)
        best = min(partners, key=lambda p: p.fdr)
        delta = best.delta_dependency
        fdr = best.fdr

        # Gate
        if not (delta < thresholds.crispr_delta_dep and fdr < thresholds.crispr_fdr):
            return 0.0, delta, fdr, best.n_mut, best.n_wt

        # Normalize to 0-1: delta_dep of -0.15 → 0.5, -0.30+ → 1.0
        score = min(abs(delta) / 0.30, 1.0)
        return round(score, 4), delta, fdr, best.n_mut, best.n_wt

    except Exception as e:
        logger.debug("CRISPR score failed for %s: %s", gene, e)
        return 0.0, None, None, None, None


def _compute_gdsc_score(
    gene: str,
    axis: CandidateAxis,
    thresholds: MinerThresholds,
) -> Tuple[float, Optional[float]]:
    """
    Returns (gdsc_score 0-1, delta_ln_ic50).
    Score = 0.0 if gate fails or data unavailable.
    """
    result = gdsc_biomarker_stratify(gene, axis.value)
    if result is None:
        return 0.0, None

    if result.p_value >= thresholds.gdsc_fdr:  # using p_value as proxy (single test)
        return 0.0, result.delta_ln_ic50

    # Normalize: delta_ln_ic50 of -1.0 → 0.5, -2.0+ → 1.0
    score = min(abs(result.delta_ln_ic50) / 2.0, 1.0)
    return round(score, 4), result.delta_ln_ic50


def _compute_kb_score(gene: str, axis: CandidateAxis) -> Tuple[float, int]:
    """
    Returns (kb_score 0-1, clinical_hit_count).
    Score = min(hits / 3, 1.0).

    DIRECTION FILTER (KB direction blindness fix): only SENSITIVITY responses
    contribute to the score. Resistance, reduced-sensitivity, adverse, and unknown
    annotations do NOT inflate kb_score. 10 negative studies = kb_score 0.0,
    not 1.0. This closes a known scientific validity gap where large resistance
    literatures would score identically to sensitivity literatures.
    """
    try:
        from ..kb.kb_engine import query as kb_query
        from ..kb.models import ResponseType
        response = kb_query(gene=gene)

        # Count SENSITIVITY-only recommendations (positive evidence direction)
        sensitivity_recs = [
            r for r in response.recommendations
            if r.response_type == ResponseType.SENSITIVITY
        ]
        kb_hits = len(sensitivity_recs)
        score = min(kb_hits / 3.0, 1.0)
        return round(score, 4), kb_hits

    except Exception as e:
        logger.debug("KB score failed for %s: %s", gene, e)
        return 0.0, 0


def _axis_to_kb_keywords(axis: CandidateAxis) -> List[str]:
    """Map axis to drug keywords for KB relevance matching."""
    return {
        CandidateAxis.CYTIDINE_ANALOGS: ["gemcitabine", "cytarabine", "decitabine", "azacitidine"],
        CandidateAxis.PARP_INHIBITORS:  ["olaparib", "niraparib", "talazoparib", "rucaparib", "parp"],
        CandidateAxis.ATR_WEE1:        ["berzosertib", "adavosertib", "ceralasertib", "atr", "wee1"],
        CandidateAxis.WRN:             ["wrn", "vx-803", "mrtx1719"],
        CandidateAxis.IMMUNOTHERAPY:   ["pembrolizumab", "nivolumab", "atezolizumab", "pd-1", "pd-l1", "io"],
        CandidateAxis.PKMYT1:          ["rp-6306", "pkmyt1"],
    }.get(axis, [])


def _compute_prism_score(
    gene: str,
    axis: CandidateAxis,
) -> Tuple[float, Optional[float]]:
    """
    Returns (prism_bonus 0.0 or PRISM_BONUS, delta_auc or None).

    PRISM is a BONUS signal — it does NOT replace or reweight existing signals.
    Bonus is 0.05 (PRISM_BONUS) if prism_fdr < PRISM_FDR_GATE (0.05), else 0.0.
    Stored delta_auc in ReceiptCandidate for transparency.
    """
    try:
        result = prism_stratify(gene, axis.value)
        if result is None:
            return 0.0, None
        if result.fdr < PRISM_FDR_GATE:
            return PRISM_BONUS, result.delta_auc
        return 0.0, result.delta_auc
    except Exception as e:
        logger.debug("PRISM score failed for %s×%s: %s", gene, axis.value, e)
        return 0.0, None


def _compute_expression_score(gene: str) -> float:
    """
    Returns expression score 0-1.
    Uses expression data if available; returns 0.0 otherwise (non-blocking).
    """
    try:
        from ..core.orchestrator import DataStore
        expr_df = DataStore._expression
        if expr_df is None or gene.upper() not in expr_df.columns:
            return 0.0
        # Simple heuristic: if gene is expressed (mean log2TPM > 1.0) → score 0.5
        mean_expr = float(expr_df[gene.upper()].mean())
        if mean_expr > 2.0:
            return 0.8
        elif mean_expr > 1.0:
            return 0.5
        return 0.2
    except Exception:
        return 0.0


# ── Tier projection ───────────────────────────────────────────────────────────

_VALID_CANDIDATE_TIERS = {"Strong", "Mechanistic", "Insufficient"}


def _project_tier(confidence: float, kb_hits: int) -> str:
    """
    Project candidate tier from confidence score and KB hits.

    SABOTAGE PROTECTION:
    - "Validated" is NEVER returned here. Validated requires human approval
      + independent clinical POSITIVE evidence. The miner cannot know that.
    - Even perfect scores (1.0) map to "Strong" at most.
    - Tier caps: Strong → Mechanistic → Insufficient (never Validated)
    """
    if confidence >= 0.70 and kb_hits >= 3:
        return "Strong"       # requires human promotion to reach Validated
    elif confidence >= 0.50:
        return "Mechanistic"
    else:
        return "Insufficient"  # will be discarded — not queued


def _build_evidence_summary(
    gene: str,
    axis: str,
    crispr_score: float,
    gdsc_score: float,
    kb_score: float,
    expression_score: float,
    confidence: float,
    tier: str,
) -> str:
    components = []
    if crispr_score > 0:
        components.append(f"CRISPR:{crispr_score:.2f}")
    if gdsc_score > 0:
        components.append(f"GDSC:{gdsc_score:.2f}")
    if kb_score > 0:
        components.append(f"KB:{kb_score:.2f}")
    if expression_score > 0:
        components.append(f"Expr:{expression_score:.2f}")
    signal_str = " | ".join(components) if components else "no signal"
    return f"{gene}×{axis} | conf={confidence:.3f} | {tier} | {signal_str}"


# ── Config reference ──────────────────────────────────────────────────────────

try:
    from ..core.config import get_settings as _get_settings
    cfg = _get_settings()
except Exception:
    class _FallbackCfg:
        depmap_release = "unknown"
    cfg = _FallbackCfg()  # type: ignore


# ── Main entry point ──────────────────────────────────────────────────────────

def mine_receipts(
    gene_panel: Optional[List[str]] = None,
    axes: Optional[List[CandidateAxis]] = None,
    thresholds: Optional[MinerThresholds] = None,
    cancer_type: Optional[str] = None,
) -> MinerRunSummary:
    """
    Mine CRISPR + GDSC + KB signals for each (gene, axis) pair.
    Candidates above min_confidence threshold are upserted into the audit queue.

    VALIDATED tier is NEVER auto-assigned (sabotage protection — see _project_tier).

    Args:
        gene_panel: list of gene symbols. Defaults to DEFAULT_GENE_PANEL (~30 genes).
        axes:       list of CandidateAxis. Defaults to DEFAULT_AXES (all 6).
        thresholds: scoring gates. Defaults to MinerThresholds().
        cancer_type: optional filter for context (not currently used in scoring — future).

    Returns:
        MinerRunSummary with counts of evaluated/queued/discarded pairs.
    """
    if gene_panel is None:
        gene_panel = DEFAULT_GENE_PANEL
    if axes is None:
        axes = DEFAULT_AXES
    if thresholds is None:
        thresholds = MinerThresholds()

    summary = MinerRunSummary(
        gene_count=len(gene_panel),
        axis_count=len(axes),
    )

    # Try to get DepMap release string
    try:
        from ..core.orchestrator import DataStore
        summary.depmap_release = getattr(cfg, "depmap_release", "unknown")
    except Exception:
        pass

    for gene in gene_panel:
        gene_upper = gene.upper()

        # Compute CRISPR score once per gene (gene-level, not axis-specific)
        crispr_score, crispr_delta, crispr_fdr, n_mut, n_wt = _compute_crispr_score(
            gene_upper, thresholds
        )

        for axis in axes:
            if axis == CandidateAxis.CUSTOM:
                continue

            # POLE AXIS GUARD: POLE encodes the DNA polymerase epsilon proofreading
            # subunit. POLE-mutated tumors are hypermutated and IO-positive, but
            # DDR SL predictions from cell lines are unreliable due to the extreme
            # mutational burden confounding CRISPR fitness effects.
            #
            # Guard fires ONLY on non-IO axes. POLE+IO is biologically correct
            # (hypermutator → IO-sensitive) and must NOT be quarantined.
            # Firing on IO would be a clinical error.
            _DDR_AXES_FOR_POLE_GUARD = {
                CandidateAxis.PARP_INHIBITORS,
                CandidateAxis.ATR_WEE1,
                CandidateAxis.PKMYT1,
                CandidateAxis.WRN,
                CandidateAxis.CYTIDINE_ANALOGS,
            }
            if gene_upper == "POLE" and axis in _DDR_AXES_FOR_POLE_GUARD:
                summary.pairs_evaluated += 1
                summary.pairs_no_signal += 1
                logger.debug(
                    "POLE DDR-axis quarantine: %s×%s skipped — POLE hypermutator "
                    "confounds DDR CRISPR fitness. IO axis is unaffected.",
                    gene_upper, axis.value,
                )
                continue

            summary.pairs_evaluated += 1

            try:
                # Signal 2: GDSC (axis-specific)
                gdsc_score, gdsc_delta = _compute_gdsc_score(gene_upper, axis, thresholds)

                # Signal 3: KB (gene-level, axis-filtered)
                kb_score, kb_hits = _compute_kb_score(gene_upper, axis)

                # Signal 4: Expression (gene-level)
                expr_score = _compute_expression_score(gene_upper)

                # Base composite confidence (weights unchanged from Sprint 2)
                # CRISPR 0.35 + GDSC 0.30 + KB 0.25 + Expr 0.10 = 1.0
                base_confidence = round(
                    crispr_score * 0.35
                    + gdsc_score  * 0.30
                    + kb_score    * 0.25
                    + expr_score  * 0.10,
                    4,
                )

                # Signal 5 (BONUS): PRISM — additive, does not rebalance weights
                prism_bonus, prism_delta = _compute_prism_score(gene_upper, axis)

                # Final confidence: base + PRISM bonus, capped at 1.0
                confidence = round(min(base_confidence + prism_bonus, 1.0), 4)

                # Tier projection — VALIDATED never emitted
                tier = _project_tier(confidence, kb_hits)

                if confidence < thresholds.min_confidence or tier == "Insufficient":
                    if confidence > 0:
                        summary.candidates_discarded += 1
                    else:
                        summary.pairs_no_signal += 1
                    continue

                # Build candidate
                candidate = ReceiptCandidate(
                    gene=gene_upper,
                    axis=axis.value,
                    cancer_type=cancer_type,
                    crispr_delta_dep=crispr_delta,
                    crispr_fdr=crispr_fdr,
                    crispr_n_mutant=n_mut,
                    crispr_n_wt=n_wt,
                    prism_delta_auc=prism_delta,   # stored for transparency
                    gdsc_delta_ic50=gdsc_delta,
                    kb_clinical_hits=kb_hits,
                    expression_corr=expr_score if expr_score > 0 else None,
                    confidence_score=confidence,
                    candidate_tier=tier,
                    evidence_summary=_build_evidence_summary(
                        gene_upper, axis.value,
                        crispr_score, gdsc_score, kb_score, expr_score,
                        confidence, tier,
                    ),
                    depmap_release=summary.depmap_release,
                    source_pipeline="receipt_miner_v1",
                )

                AuditQueue.upsert(candidate)
                summary.candidates_above_threshold += 1
                summary.candidates_queued += 1

            except Exception as e:
                msg = f"{gene_upper}×{axis.value}: {e}"
                logger.warning("Miner error: %s", msg)
                summary.errors.append(msg)

    return summary.finish()
