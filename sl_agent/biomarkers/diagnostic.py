"""
Diagnostic biomarker engine (subtype classifier wrapper).

First instance: consensusOV (Chen et al. 2018) HGSOC molecular subtypes
(differentiated / immunoreactive / mesenchymal / proliferative).

We reuse the published classifier verbatim (R). This module builds a
BiomarkerModel from the R-produced concordance receipts and carries the
documented composition/TME confound as an explicit caveat — the subtype call
is a *composition-aware* classifier, not a claim of intrinsic epithelial subtype.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

from .models import BiomarkerModel, BiomarkerStatus, BiomarkerType, ValidationReceipt

COMPOSITION_CAVEAT = (
    "HGSOC transcriptomic subtypes are partly driven by tumour cellular composition / "
    "microenvironment and subclonal mixing (Geistlinger 2020; Tanis 2026); treat as a "
    "composition-aware classifier, not an intrinsic epithelial subtype."
)


def model_from_receipts(
    cancer: str,
    method_id: str,
    receipts_json: str | Path,
    source_pmid: Optional[str] = None,
    source_doi: Optional[str] = None,
    source_note: str = "",
) -> BiomarkerModel:
    """
    Build a diagnostic BiomarkerModel from an R-produced concordance receipts JSON.

    Expected JSON shape (written by w3 script):
      {
        "gene_panel_size": int,             # classifier feature count (informational)
        "cohorts": [{"cohort_id":..., "n":..., "subtype_counts": {...}}, ...],
        "concordance": {"cohort_a":..., "cohort_b":..., "metric_name":"adjusted_rand"|"pct_agree",
                        "metric_value":..., "ci95_low":..., "ci95_high":...,
                        "is_independent": true},
        "interpretation": "...",
        "caveats": [...]
      }
    Cross-cohort concordance on INDEPENDENT cohorts is the external-replication
    analogue for a classifier.
    """
    d = json.loads(Path(receipts_json).read_text())
    receipts: List[ValidationReceipt] = []

    for c in d.get("cohorts", []):
        receipts.append(ValidationReceipt(
            kind="train" if c.get("is_reference") else "internal_cv",
            cohort_id=c["cohort_id"],
            n=int(c["n"]),
            is_independent_of_training=False,
            method_detail="subtype distribution",
            extra={"subtype_counts": c.get("subtype_counts", {})},
        ))

    conc = d.get("concordance")
    ext = None
    if conc:
        ext = ValidationReceipt(
            kind="external_replication",
            cohort_id=f"{conc.get('cohort_a')}__vs__{conc.get('cohort_b')}",
            n=int(conc.get("n", 0)),
            is_independent_of_training=bool(conc.get("is_independent", False)),
            metric_name=conc.get("metric_name"),
            metric_value=conc.get("metric_value"),
            ci95_low=conc.get("ci95_low"),
            ci95_high=conc.get("ci95_high"),
            method_detail="cross-cohort subtype concordance",
        )
        receipts.append(ext)

    status = (
        BiomarkerStatus.VALIDATED
        if ext is not None and ext.is_independent_of_training and (ext.metric_value or 0) > 0
        else BiomarkerStatus.DISCOVERY_ONLY
    )

    caveats = list(d.get("caveats", []))
    if COMPOSITION_CAVEAT not in caveats:
        caveats.append(COMPOSITION_CAVEAT)

    return BiomarkerModel(
        cancer=cancer,
        biomarker_type=BiomarkerType.DIAGNOSTIC,
        method_id=method_id,
        source_pmid=source_pmid,
        source_doi=source_doi,
        source_note=source_note,
        gene_panel=[],  # classifier is a fixed published model, not a gene list we own
        status=status,
        receipts=receipts,
        interpretation=d.get("interpretation", ""),
        caveats=caveats,
    )
