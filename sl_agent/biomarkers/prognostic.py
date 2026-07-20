"""
Prognostic biomarker engine (cancer-agnostic).

Reproduces a published survival signature (first instance: Riester et al. 2014,
JNCI, OV) as a penalised-Cox risk score. Feature selection happens INSIDE CV
folds; the primary reported metric is out-of-fold / external, never in-sample.

Heavy lifting (curated multi-cohort data, model fit) runs in R via the worker
scripts and writes JSON receipts; this module builds BiomarkerModel objects
from those receipts and applies a frozen score to new expression.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .models import BiomarkerModel, BiomarkerStatus, BiomarkerType, ValidationReceipt


def apply_risk_score(
    expr_z: "np.ndarray",
    gene_index: Dict[str, int],
    coefficients: Dict[str, float],
) -> "np.ndarray":
    """
    Apply a frozen linear (Cox) risk score to z-scored expression.

    expr_z       : samples x genes matrix (z-scored within cohort)
    gene_index   : HUGO -> column index in expr_z
    coefficients : HUGO -> Cox coefficient (log-HR)
    Returns a per-sample risk score (higher = worse prognosis).
    """
    score = np.zeros(expr_z.shape[0], dtype=float)
    used = 0
    for gene, coef in coefficients.items():
        j = gene_index.get(gene)
        if j is None:
            continue
        score += coef * expr_z[:, j]
        used += 1
    if used == 0:
        raise ValueError("No signature genes found in the provided expression matrix.")
    return score


def model_from_receipts(
    cancer: str,
    method_id: str,
    receipts_json: str | Path,
    source_pmid: Optional[str] = None,
    source_doi: Optional[str] = None,
    source_note: str = "",
) -> BiomarkerModel:
    """
    Build a prognostic BiomarkerModel from an R-produced receipts JSON.

    Expected JSON shape (written by w1/w2 scripts):
      {
        "gene_panel": [...],
        "train": {"cohort_id":..., "n":..., "n_events":..., "metric_name":"c_index",
                  "metric_value":..., "seed":..., "epv":...},
        "internal_cv": {... optional ...},
        "external": {"cohort_id":..., "n":..., "n_events":..., "metric_name":"c_index",
                     "metric_value":..., "ci95_low":..., "ci95_high":..., "p_value":...,
                     "ph_test_p":..., "is_independent": true/false},
        "interpretation": "...",
        "caveats": [...]
      }
    """
    d = json.loads(Path(receipts_json).read_text())
    receipts: List[ValidationReceipt] = []

    def _mk(kind: str, blk: dict) -> ValidationReceipt:
        return ValidationReceipt(
            kind=kind,
            cohort_id=blk["cohort_id"],
            n=int(blk["n"]),
            n_events=blk.get("n_events"),
            is_independent_of_training=bool(blk.get("is_independent", kind == "train" and False)),
            metric_name=blk.get("metric_name"),
            metric_value=blk.get("metric_value"),
            ci95_low=blk.get("ci95_low"),
            ci95_high=blk.get("ci95_high"),
            p_value=blk.get("p_value"),
            ph_test_p=blk.get("ph_test_p"),
            epv=blk.get("epv"),
            seed=blk.get("seed"),
            method_detail=blk.get("method_detail"),
            extra=blk.get("extra", {}),
        )

    if "train" in d:
        receipts.append(_mk("train", d["train"]))
    if "internal_cv" in d:
        receipts.append(_mk("internal_cv", d["internal_cv"]))
    if "external" in d:
        ext = _mk("external_replication", d["external"])
        receipts.append(ext)
    if "sensitivity_check" in d:
        receipts.append(_mk("sensitivity_check", d["sensitivity_check"]))

    ext = next((r for r in receipts if r.kind == "external_replication"), None)
    status = (
        BiomarkerStatus.VALIDATED
        if ext is not None and ext.is_independent_of_training
        else BiomarkerStatus.DISCOVERY_ONLY
    )

    return BiomarkerModel(
        cancer=cancer,
        biomarker_type=BiomarkerType.PROGNOSTIC,
        method_id=method_id,
        source_pmid=source_pmid,
        source_doi=source_doi,
        source_note=source_note,
        gene_panel=d.get("gene_panel", []),
        status=status,
        receipts=receipts,
        interpretation=d.get("interpretation", ""),
        caveats=d.get("caveats", []),
    )
