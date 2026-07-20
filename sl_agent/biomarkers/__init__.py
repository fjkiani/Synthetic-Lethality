"""
Biomarker subsystem — prognostic & diagnostic capabilities for the v3 spine.

Doctrine (carried from the PI3K/mTOR audit):
  - A biomarker is "validated" ONLY if it reproduces in an independent,
    non-overlapping cohort. Same-cohort reprocessing is a sensitivity check.
  - Everything else is "discovery_only" and labelled as such.
  - Every number in a shipped interpretation must trace to a written receipt.
  - RESEARCH USE ONLY. Not validated for clinical decision-making.

This subsystem is *advisory*: biomarker results attach to the multimodal
evidence surface as context and NEVER override synthetic-lethality drug
routing.
"""
from .models import (  # noqa: F401
    BiomarkerType,
    BiomarkerStatus,
    ValidationReceipt,
    BiomarkerModel,
    BiomarkerResult,
    BiomarkerContext,
)
from .registry import BiomarkerRegistry  # noqa: F401

__all__ = [
    "BiomarkerType",
    "BiomarkerStatus",
    "ValidationReceipt",
    "BiomarkerModel",
    "BiomarkerResult",
    "BiomarkerContext",
    "BiomarkerRegistry",
]
