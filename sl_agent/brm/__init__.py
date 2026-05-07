"""
sl_agent.brm — Brain Metastasis Targetability Engine v2

Discovery-first: expression-stratified SL from DepMap 24Q4,
FDR-weighted novelty scoring, 3-section partitioned output.
"""

from sl_agent.brm.models import (
    BrMRowClass,
    ScoreBasis,
    ExploitMode,
    BrMEvidenceClass,
    CascadeStep,
    TargetNodeRole,
    ContextSpecificity,
    SLPartner,
    BrMGeneRow,
    BrMTargetabilityMatrix,
    BrMQueryInput,
    RUO_DISCLAIMER,
)
from sl_agent.brm.matrix_builder import build_brm_matrix

__all__ = [
    "BrMRowClass",
    "ScoreBasis",
    "ExploitMode",
    "BrMEvidenceClass",
    "CascadeStep",
    "TargetNodeRole",
    "ContextSpecificity",
    "SLPartner",
    "BrMGeneRow",
    "BrMTargetabilityMatrix",
    "BrMQueryInput",
    "RUO_DISCLAIMER",
    "build_brm_matrix",
]
