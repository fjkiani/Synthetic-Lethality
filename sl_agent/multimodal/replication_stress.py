"""
Replication Stress (RS) Scoring Module — SL Agent v4.1

Computes a transparent, weighted RS score from patient genomic features.
Used exclusively to modulate tier promotion for ATR/WEE1 and PKMYT1 axes.

Design constraints (DO NOT RELAX):
  - MBD4 alone (0.05) → RS = "low". Cannot reach "moderate" or "high".
  - RS can NEVER promote an axis directly to "Validated".
  - RS can only promote Mechanistic → Strong, AND only when CRISPR or pharma
    evidence already exists (RS is a co-factor, not a standalone promoter).
  - All weights and thresholds are explicitly documented here.

Documented weights (sum to 1.10 — intentionally > 1 to allow partial feature
scores while keeping normalized output in 0.0–1.0 range):
  CCNE1_amp   0.30  strongest single RS driver (replication origin over-firing)
  MYC_amp     0.20  oncogene-induced replication stress
  ARID1A_LOF  0.20  chromatin remodeling loss → fork instability
  MSI_H       0.15  dMMR → elevated RS, fork stalling at microsatellites
  TP53_LOF    0.10  G1/S checkpoint loss → replication without arrest
  high_ploidy 0.10  genome doubling → centrosome stress, mitotic RS
  MBD4_LOF    0.05  BER deficiency → modest fork lesion accumulation

Level thresholds:
  score ≥ 0.40 → "high"
  score ≥ 0.20 → "moderate"
  score >  0.00 → "low"
  score == 0.00 → "none"

Validation cases (enforced by test_rs_benchmark.py):
  MBD4 alone                           → 0.05  → "low"
  MBD4 + TP53                          → 0.15  → "low"
  MBD4 + TP53 + ARID1A                 → 0.35  → "moderate"
  CCNE1 + TP53                         → 0.40  → "high"
  CCNE1 + MYC + ARID1A                 → 0.70  → "high"  (capped at 1.0)
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Optional

if TYPE_CHECKING:
    from .models import EvidenceMatrix, CombinationVulnerability

from pydantic import BaseModel, Field


# ── Feature set ───────────────────────────────────────────────────────────────

class RSFeatureSet(BaseModel):
    """
    Binary + optional continuous genomic features used to compute RS score.
    All boolean fields default False (absent).
    """
    # Binary features (each contributes its documented weight if True)
    ccne1_amplified: bool = Field(False, description="CCNE1 amp → premature S-phase / origin over-firing")
    myc_amplified:   bool = Field(False, description="MYC amp → transcription-replication collision, RS")
    arid1a_lof:      bool = Field(False, description="Chromatin remodeling loss → fork instability")
    tp53_lof:        bool = Field(False, description="G1/S checkpoint loss → cells enter S-phase with damage")
    mbd4_lof:        bool = Field(False, description="BER deficiency → modest fork lesion accumulation")
    msi_h:           bool = Field(False, description="dMMR / MSI-H → fork stalling at microsatellites")
    high_ploidy:     bool = Field(False, description="Aneuploidy / genome doubling → centrosome + mitotic RS")

    # Optional continuous features (reserved for future quantitative integration)
    tmb:          Optional[float] = Field(None, description="Tumor mutational burden (mut/Mb)")
    ploidy_score: Optional[float] = Field(None, description="Numeric ploidy value if available")


# ── RS score output ───────────────────────────────────────────────────────────

class RSScore(BaseModel):
    """
    Computed Replication Stress score with full audit trail.
    Transparent: every weight, contributing feature, and threshold is documented.
    """
    score: float = Field(..., description="Weighted sum, capped at 1.0, floor 0.0")
    level: str   = Field(..., description="'high' / 'moderate' / 'low' / 'none'")
    contributing_features: List[str]       = Field(default_factory=list)
    weights_applied:       Dict[str, float] = Field(default_factory=dict)
    rationale: str = Field(..., description="One-line human rationale for the score and level")


# ── Documented weights ────────────────────────────────────────────────────────

_RS_WEIGHTS: Dict[str, float] = {
    "ccne1_amplified": 0.30,
    "myc_amplified":   0.20,
    "arid1a_lof":      0.20,
    "msi_h":           0.15,
    "tp53_lof":        0.10,
    "high_ploidy":     0.10,
    "mbd4_lof":        0.05,   # deliberately lowest — BER modest RS contributor
}

# Level thresholds
_HIGH_THRESHOLD     = 0.40
_MODERATE_THRESHOLD = 0.20


# ── Core function ─────────────────────────────────────────────────────────────

def compute_rs_score(features: RSFeatureSet) -> RSScore:
    """
    Compute a transparent weighted RS score from binary genomic features.

    Algorithm:
      1. For each True boolean feature, accumulate its documented weight.
      2. Cap result at 1.0, floor at 0.0.
      3. Assign level: high / moderate / low / none.
      4. Record all contributing features and weights for audit.

    CRITICAL CONSTRAINTS:
      - MBD4 alone: score = 0.05 → "low"   (verified by test)
      - MBD4 + TP53: score = 0.15 → "low"   (verified by test)
      - MBD4 + TP53 + ARID1A: score = 0.35 → "moderate" (verified by test)
      - CCNE1 + TP53: score = 0.40 → "high"  (verified by test)
    """
    raw_score = 0.0
    contributing: List[str] = []
    weights_applied: Dict[str, float] = {}

    feature_dict = {
        "ccne1_amplified": features.ccne1_amplified,
        "myc_amplified":   features.myc_amplified,
        "arid1a_lof":      features.arid1a_lof,
        "msi_h":           features.msi_h,
        "tp53_lof":        features.tp53_lof,
        "high_ploidy":     features.high_ploidy,
        "mbd4_lof":        features.mbd4_lof,
    }

    for feature_name, is_present in feature_dict.items():
        if is_present:
            w = _RS_WEIGHTS[feature_name]
            raw_score += w
            contributing.append(feature_name)
            weights_applied[feature_name] = w

    # Cap at 1.0
    score = min(1.0, max(0.0, raw_score))

    # Level assignment
    if score >= _HIGH_THRESHOLD:
        level = "high"
    elif score >= _MODERATE_THRESHOLD:
        level = "moderate"
    elif score > 0.0:
        level = "low"
    else:
        level = "none"

    # Build rationale
    if not contributing:
        rationale = "No RS features present. RS score = 0 (none)."
    else:
        feature_str = " + ".join(contributing)
        rationale = (
            f"RS features: {feature_str}. "
            f"Weighted sum = {raw_score:.2f} → capped score = {score:.2f} → level = '{level}'. "
            f"{'MBD4 contributes modest RS (0.05 weight); co-mutations required to reach moderate/high.' if features.mbd4_lof and not features.ccne1_amplified else ''}"
        ).strip()

    return RSScore(
        score=round(score, 4),
        level=level,
        contributing_features=contributing,
        weights_applied=weights_applied,
        rationale=rationale,
    )


# ── Convenience predicates ────────────────────────────────────────────────────

def rs_is_high(rs: Optional[RSScore]) -> bool:
    return rs is not None and rs.level == "high"


def rs_is_moderate_or_higher(rs: Optional[RSScore]) -> bool:
    return rs is not None and rs.level in ("moderate", "high")


def rs_is_actionable(rs: Optional[RSScore]) -> bool:
    """True if RS score is high enough to be considered in tier promotion logic."""
    return rs_is_high(rs)


# ── Combinatorial vulnerability assessment ────────────────────────────────────

# Genes known to drive BER/DDR deficiency relevant to RS
_BER_DDR_GENES = {"MBD4", "MLH1", "MSH2", "MSH6", "BRCA1", "BRCA2"}

# Checkpoint axes that can form combinatorial vulnerability
_CHECKPOINT_AXES_FOR_COMBINATION = {
    "atr_wee1": ["ATR", "CHK1", "WEE1"],
    "pkmyt1":   ["PKMYT1"],
}

def assess_combinatorial_vulnerability(
    gene: str,
    rs_score: "RSScore",
    evidence_matrix: "EvidenceMatrix",
) -> Optional["CombinationVulnerability"]:
    """
    Check if a gene's BER/DDR deficiency creates a combinatorial vulnerability
    with checkpoint kinases, where neither pathway alone is SL.

    Logic:
      1. Gene must be in known BER/DDR deficiency set (or rs_score > 0)
      2. ATR/WEE1 axis must be "Mechanistic candidate only" individually (not already Strong/Validated)
      3. There must be a mechanistic chain: gene_LOF → fork_stress → checkpoint_dependency
      4. Identify exact data gaps required to unlock promotion

    Returns CombinationVulnerability or None if not applicable.

    CONSTRAINT: This function ONLY identifies the vulnerability — it does NOT modify tiers.
    Tier promotion still requires independent CRISPR or pharma evidence via fuse_matrix().
    """
    # Late import to avoid circular dependency
    from .models import CombinationVulnerability

    gene_upper = gene.upper()

    # Must have some RS contribution
    if rs_score.score == 0.0:
        return None

    # Find ATR/WEE1 row in matrix
    atr_row = None
    for row in evidence_matrix.rows:
        if row.axis.value == "atr_wee1":
            atr_row = row
            break

    if atr_row is None:
        return None

    # Both pathways must be below Strong for a "combinatorial" story to make sense
    atr_tier = atr_row.recommendation_tier or "Unknown"
    if atr_tier in ("Validated SL therapeutic lever", "Strong candidate dependency axis"):
        # Already realized — still worth reporting, but as REALIZED combinatorial
        combined_tier = atr_tier
        data_gaps: List[str] = []
        convergence = "replication_fork_collapse"
    else:
        # Not yet realized — identify gaps
        combined_tier = None
        data_gaps = _identify_data_gaps(gene_upper, atr_row)
        convergence = "mitotic_catastrophe"

    # Build mechanistic rationale
    if gene_upper == "MBD4":
        rationale = (
            f"MBD4 LOF → BER collapse → 5mC>T lesions → replication fork stalling → "
            f"ATR/CHK1 checkpoint dependency → if checkpoint inhibited → mitotic catastrophe. "
            f"RS score = {rs_score.score:.2f} ({rs_score.level}). "
            f"{'Data gaps remain before combinatorial tier can be promoted.' if data_gaps else 'Combinatorial tier realized by available evidence.'}"
        )
    elif gene_upper == "ARID1A":
        rationale = (
            f"ARID1A LOF → chromatin access failure → fork instability → "
            f"ATR dependency (Williamson 2016). RS score = {rs_score.score:.2f} ({rs_score.level}). "
            f"{'Combinatorial tier unlocked by available CRISPR + functional data.' if not data_gaps else 'Additional evidence needed.'}"
        )
    else:
        rationale = (
            f"{gene_upper} LOF contributes to replication stress (RS={rs_score.score:.2f}, {rs_score.level}), "
            f"creating potential checkpoint dependency. "
            f"{'Data gaps identified for combinatorial tier promotion.' if data_gaps else 'Evidence sufficient for combinatorial assessment.'}"
        )

    return CombinationVulnerability(
        pathway_a=_gene_to_pathway(gene_upper),
        pathway_a_gene=gene_upper,
        pathway_b="checkpoint_inhibition",
        pathway_b_targets=["ATR", "CHK1", "WEE1"],
        convergence_mechanism=convergence,
        pathway_a_alone_tier=_get_gene_alone_tier(gene_upper),
        pathway_b_alone_tier=atr_tier,
        combined_tier=combined_tier,
        data_gaps=data_gaps,
        rs_score=rs_score.score,
        rs_level=rs_score.level,
        rationale=rationale,
    )


def _gene_to_pathway(gene: str) -> str:
    """Map gene to its deficiency pathway label."""
    pathway_map = {
        "MBD4":  "BER_deficiency",
        "ARID1A": "chromatin_remodeling_loss",
        "BRCA1": "HR_deficiency",
        "BRCA2": "HR_deficiency",
        "MLH1":  "MMR_deficiency",
        "MSH2":  "MMR_deficiency",
        "MSH6":  "MMR_deficiency",
        "CCNE1": "replication_stress_amplification",
    }
    return pathway_map.get(gene, f"{gene}_deficiency")


def _get_gene_alone_tier(gene: str) -> str:
    """Conservative standalone tier for pathway A gene (without combination)."""
    # These are the tier the gene's OWN axis gets individually
    # For BER/chromatin genes, standalone checkpoint dependency is typically Mechanistic
    return "Mechanistic candidate only"


def _identify_data_gaps(gene: str, atr_row) -> List[str]:
    """Identify what data is missing to unlock the combinatorial tier."""
    from .models import ModalityStatus
    gaps = []
    cells = atr_row.cells()

    if cells["crispr"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"ATR/WEE1 CRISPR dependency in {gene}-LOF lines "
            f"(DepMap stratification: Chronos Δdep in {gene}-LOF vs WT)"
        )
    if cells["prism"].status != ModalityStatus.POSITIVE and cells["gdsc"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"ATR/WEE1 drug sensitivity stratified by {gene} status "
            f"(GDSC/PRISM: ceralasertib / adavosertib / berzosertib AUC in {gene}-LOF vs WT)"
        )
    if cells["in_vitro"].status != ModalityStatus.POSITIVE:
        gaps.append(
            f"Isogenic {gene}-KO + ATRi synergy experiment "
            f"(Bliss/Loewe combination index in {gene}-KO vs WT)"
        )
    if not gaps:
        # All present — tier should have been promoted already
        gaps = []  # REALIZED

    return gaps
