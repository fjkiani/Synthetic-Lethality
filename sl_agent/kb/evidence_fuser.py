"""
Evidence fusion and tiering engine.

Takes raw KBEvidence records from multiple sources and:
  1. Deduplicates by (gene, drug) key
  2. Counts how many KBs agree
  3. Upgrades tier when ≥2 KBs agree (consensus boost)
  4. Assigns confidence: high / medium / low
  5. Returns ranked DrugRecommendation list

Mars tiering contract:
  Tier A  — FDA-approved / guideline-standard + confirmed by ≥1 KB
  Tier B  — ≥2 KBs with clinical trial evidence (Phase II/III)
  Tier C  — Single KB OR Phase I only
  Tier D  — Preclinical / in-vitro / computational only
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .models import (
    DrugRecommendation,
    EvidenceTier,
    KBEvidence,
    KBQueryInput,
    ResponseType,
    SourceKB,
)

logger = logging.getLogger(__name__)

# Tier ranking (higher = stronger)
_TIER_RANK = {
    EvidenceTier.A: 4,
    EvidenceTier.B: 3,
    EvidenceTier.C: 2,
    EvidenceTier.D: 1,
    EvidenceTier.UNKNOWN: 0,
}

_RANK_TO_TIER = {v: k for k, v in _TIER_RANK.items()}


def _normalize_drug_key(drug: str) -> str:
    """Case-fold and strip salt/brand suffixes for dedup."""
    return (
        drug.lower()
        .replace(" hydrochloride", "")
        .replace(" mesylate", "")
        .replace(" tosylate", "")
        .strip()
    )


def _consensus_tier(tier: EvidenceTier, n_kbs: int) -> EvidenceTier:
    """
    Boost tier if multiple KBs agree.
    If ≥2 KBs say Tier C → promote to Tier B.
    If ≥2 KBs say Tier B → stays Tier B (already strong).
    """
    current_rank = _TIER_RANK.get(tier, 0)
    if n_kbs >= 2 and current_rank == _TIER_RANK[EvidenceTier.C]:
        return EvidenceTier.B
    return tier


def _confidence(tier: EvidenceTier, n_kbs: int, n_evidence: int) -> str:
    if tier == EvidenceTier.A:
        return "high"
    if tier == EvidenceTier.B and n_kbs >= 2:
        return "high"
    if tier == EvidenceTier.B or (tier == EvidenceTier.C and n_kbs >= 2):
        return "medium"
    if n_evidence >= 3:
        return "medium"
    return "low"


# ── Main fusion function ──────────────────────────────────────────────────────

def fuse_evidence(
    all_evidence: List[KBEvidence],
    query: KBQueryInput,
) -> List[DrugRecommendation]:
    """
    Fuse evidence from multiple KBs into ranked DrugRecommendations.

    Returns list sorted by tier (desc) then KB consensus count (desc).
    """
    if not all_evidence:
        return []

    # Group by (gene, normalized_drug_key)
    groups: Dict[Tuple[str, str], List[KBEvidence]] = defaultdict(list)
    for ev in all_evidence:
        key = (ev.gene.upper(), _normalize_drug_key(ev.drug))
        groups[key].append(ev)

    recommendations: List[DrugRecommendation] = []

    for (gene_key, drug_key), evidence_list in groups.items():
        # Filter by min_tier
        min_rank = _TIER_RANK.get(query.min_tier, 0)
        max_tier_in_group = max(
            (_TIER_RANK.get(e.tier, 0) for e in evidence_list), default=0
        )
        if max_tier_in_group < min_rank:
            continue

        # Response type: pick most common, prefer sensitivity
        response_counts: Dict[ResponseType, int] = defaultdict(int)
        for e in evidence_list:
            response_counts[e.response_type] += 1
        dominant_response = max(response_counts, key=lambda r: (
            response_counts[r],
            1 if r == ResponseType.SENSITIVITY else 0,
        ))

        # Skip resistance if not requested
        if dominant_response == ResponseType.RESISTANCE and not query.include_resistance:
            continue

        # Tier: take the strongest tier seen, then apply consensus boost
        best_tier_rank = max(_TIER_RANK.get(e.tier, 0) for e in evidence_list)
        best_tier = _RANK_TO_TIER.get(best_tier_rank, EvidenceTier.UNKNOWN)
        supporting_kbs = list({e.source_kb for e in evidence_list})
        final_tier = _consensus_tier(best_tier, len(supporting_kbs))

        # Drug name: pick the most common / longest
        drug_names = [e.drug for e in evidence_list]
        representative_drug = max(set(drug_names), key=lambda d: (drug_names.count(d), len(d)))

        # Drug class: prefer non-null
        drug_classes = [e.drug_class for e in evidence_list if e.drug_class]
        drug_class = drug_classes[0] if drug_classes else None

        # Collect PMIDs, trials, summaries, diseases
        pmids = list({p for e in evidence_list for p in (e.pmids or [])})
        trials = list({t for e in evidence_list for t in (e.clinical_trials or [])})
        summaries = [e.evidence_statement for e in evidence_list if e.evidence_statement][:5]
        diseases = list({e.disease for e in evidence_list if e.disease})

        # Variant (for context)
        variants = list({e.variant for e in evidence_list if e.variant})
        variant_str = variants[0] if len(variants) == 1 else (
            ", ".join(sorted(variants)) if variants else None
        )

        confidence = _confidence(final_tier, len(supporting_kbs), len(evidence_list))

        recommendations.append(
            DrugRecommendation(
                gene=gene_key,
                variant=variant_str,
                drug=representative_drug,
                drug_class=drug_class,
                response_type=dominant_response,
                tier=final_tier,
                confidence=confidence,
                supporting_kbs=supporting_kbs,
                num_evidence_items=len(evidence_list),
                pmids=pmids[:10],
                clinical_trials=trials[:5],
                evidence_summaries=summaries,
                cancer_types=diseases,
            )
        )

    # Sort: Tier (desc) → KB count (desc) → evidence count (desc)
    recommendations.sort(
        key=lambda r: (
            -_TIER_RANK.get(r.tier, 0),
            -len(r.supporting_kbs),
            -r.num_evidence_items,
        )
    )

    return recommendations[: query.top_n]
