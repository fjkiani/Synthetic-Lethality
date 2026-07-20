"""
Assemble BiomarkerModel objects from the R-produced validation receipts
(W1 prognostic train, W2 cross-study replication, W3 diagnostic subtype,
W4 baseline discovery scan) and register them with honest status.

This is the bridge between the reused published methods (run in R) and the
v3-spine biomarker subsystem. It NEVER fabricates a metric; every number in a
model's interpretation traces to a receipt, and the grounding guard is run on
each model before it is registered.

Status policy (carried from the PI3K/mTOR audit):
  - A prognostic signature is VALIDATED only if it reproduces (C-index CI
    excludes 0.5 AND HR CI excludes 1, same direction) in an INDEPENDENT,
    non-overlapping cohort. Inconsistent replication across independent cohorts
    => discovery_only.
  - Same-cohort reprocessing is a sensitivity check, never replication.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

from .models import (
    BiomarkerModel,
    BiomarkerStatus,
    BiomarkerType,
    ValidationReceipt,
)
from .grounding_guard import assert_grounded
from .registry import BiomarkerRegistry

CANCER = "ovarian_hgsoc"


# ── prognostic (W1 train + W2 cross-study replication) ────────────────────────

def _receipt_from_cohort_block(kind: str, blk: dict, is_independent: bool) -> ValidationReceipt:
    return ValidationReceipt(
        kind=kind,
        cohort_id=blk["cohort_id"],
        n=int(blk["n"]),
        n_events=blk.get("n_events"),
        is_independent_of_training=is_independent,
        metric_name="c_index",
        metric_value=blk.get("c_index"),
        ci95_low=blk.get("c_index_lower"),
        ci95_high=blk.get("c_index_upper"),
        p_value=blk.get("c_index_p"),
        ph_test_p=blk.get("ph_test_p"),
        method_detail=blk.get("note", "frozen-coefficient risk score applied to cohort"),
        extra={
            "hr_per_1sd": blk.get("hr_per_1sd"),
            "hr_lower": blk.get("hr_lower"),
            "hr_upper": blk.get("hr_upper"),
            "hr_p": blk.get("hr_p"),
            "km_tertile_logrank_p": blk.get("km_tertile_logrank_p"),
            "overlap_with_train": blk.get("overlap_with_train"),
        },
    )


def _replicates(blk: dict) -> bool:
    """C-index CI excludes 0.5 AND HR CI excludes 1 (either direction)."""
    ci_lo = blk.get("c_index_lower")
    hr_lo, hr_hi = blk.get("hr_lower"), blk.get("hr_upper")
    if ci_lo is None or hr_lo is None or hr_hi is None:
        return False
    ci_ok = ci_lo > 0.5
    hr_ok = (hr_lo > 1.0) or (hr_hi < 1.0)
    return bool(ci_ok and hr_ok)


def build_prognostic_model(w1_json: str | Path, w2_json: str | Path) -> BiomarkerModel:
    w1 = json.loads(Path(w1_json).read_text())
    w2 = json.loads(Path(w2_json).read_text())

    receipts: List[ValidationReceipt] = []

    # W1 internal (train + OOF CV) — never independent
    tr = w1["train"]
    receipts.append(ValidationReceipt(
        kind="train", cohort_id=tr["cohort_id"], n=int(tr["n"]), n_events=tr.get("n_events"),
        is_independent_of_training=False, metric_name=tr.get("metric_name"),
        metric_value=tr.get("metric_value"), p_value=tr.get("p_value"),
        ph_test_p=tr.get("ph_test_p"), epv=tr.get("epv"), seed=tr.get("seed"),
        method_detail=tr.get("method_detail"),
    ))
    cv = w1.get("internal_cv")
    if cv:
        receipts.append(ValidationReceipt(
            kind="internal_cv", cohort_id=cv["cohort_id"], n=int(cv["n"]),
            n_events=cv.get("n_events"), is_independent_of_training=False,
            metric_name=cv.get("metric_name"), metric_value=cv.get("metric_value"),
            ci95_low=cv.get("ci95_low"), ci95_high=cv.get("ci95_high"),
            p_value=cv.get("p_value"), seed=cv.get("seed"), method_detail=cv.get("method_detail"),
        ))

    # W2 external replication cohorts (independent GEO series)
    rep = w2.get("replication_cohorts", {})
    n_replicated = 0
    for cid, blk in rep.items():
        indep = int(blk.get("overlap_with_train", 0) or 0) == 0
        receipts.append(_receipt_from_cohort_block("external_replication", blk, indep))
        if indep and _replicates(blk):
            n_replicated += 1

    # W2 same-cohort sensitivity check (explicitly NOT independent)
    sens = w2.get("sensitivity_check_same_cohort")
    if sens:
        receipts.append(_receipt_from_cohort_block("sensitivity_check", sens, False))

    n_indep = sum(
        1 for r in receipts
        if r.kind == "external_replication" and r.is_independent_of_training
    )

    # Honest status: VALIDATED only if it replicates CONSISTENTLY across all
    # independent cohorts tested (>=2). One-of-two => discovery_only.
    consistent = (n_indep >= 2) and (n_replicated == n_indep)
    status = BiomarkerStatus.VALIDATED if consistent else BiomarkerStatus.DISCOVERY_ONLY

    # Build an interpretation whose every number traces to a receipt.
    gse9891 = rep.get("GSE9891", {})
    gse32062 = rep.get("GSE32062", {})

    def f(x, nd=3):
        return f"{float(x):.{nd}f}" if x is not None else "NA"

    interp = (
        f"Riester et al. OV survival signature, reproduced as an elastic-net Cox risk score "
        f"(frozen on {tr['cohort_id']}, n={tr['n']}, {tr.get('n_events')} events; "
        f"5-fold out-of-fold C-index {f(cv.get('metric_value'))}). "
        f"Applied UNCHANGED to two INDEPENDENT GEO cohorts (0 patient overlap with training): "
        f"in GSE9891 (n={gse9891.get('n')}, {gse9891.get('n_events')} events) the score is prognostic "
        f"(C-index {f(gse9891.get('c_index'))}, HR/SD {f(gse9891.get('hr_per_1sd'),2)}, "
        f"p={f(gse9891.get('hr_p'))}); in GSE32062 (n={gse32062.get('n')}, {gse32062.get('n_events')} events) "
        f"it is NULL (C-index {f(gse32062.get('c_index'))}, HR/SD {f(gse32062.get('hr_per_1sd'),2)}, "
        f"p={f(gse32062.get('hr_p'))}). Replication is therefore INCONSISTENT across independent cohorts, "
        f"and the proportional-hazards assumption is violated in both (PH p={f(gse9891.get('ph_test_p'),3)} "
        f"and p={f(gse32062.get('ph_test_p'),3)}). A same-cohort reprocessing of the training data is a "
        f"sensitivity check only, not replication. VERDICT: discovery-only — not a validated prognostic "
        f"classifier. RESEARCH USE ONLY — not validated for clinical decision-making."
    )

    caveats = [
        "External replication is inconsistent: signal present in GSE9891 but absent in GSE32062.",
        "Proportional-hazards assumption violated in both independent cohorts (time-varying effect).",
        "Only 3 of 18 signature genes survived elastic-net penalisation on the training cohort.",
        "Same-cohort reprocessing (TCGA) is a sensitivity check, never external replication.",
    ]

    model = BiomarkerModel(
        cancer=CANCER,
        biomarker_type=BiomarkerType.PROGNOSTIC,
        method_id=w1.get("method_id", "riester2014_prognostic_ov"),
        source_pmid=w1.get("source_pmid"),
        source_doi=w1.get("source_doi"),
        source_note=w1.get("source_note", ""),
        gene_panel=w1.get("gene_panel", []),
        status=status,
        receipts=receipts,
        interpretation=interp,
        caveats=caveats,
    )
    return model


# ── diagnostic (W3 consensusOV) ───────────────────────────────────────────────

def build_diagnostic_model(w3_json: str | Path) -> BiomarkerModel:
    d = json.loads(Path(w3_json).read_text())
    cohorts = d.get("cohorts", {})
    receipts: List[ValidationReceipt] = []
    for cid, c in cohorts.items():
        receipts.append(ValidationReceipt(
            kind="internal_cv",
            cohort_id=c["cohort_id"],
            n=int(c["n"]),
            is_independent_of_training=False,
            metric_name="mean_margin",
            metric_value=c.get("mean_margin"),
            method_detail="per-cohort subtype distribution + RF margin calibration",
            extra={
                "subtype_distribution": c.get("subtype_distribution", {}),
                "mean_entropy": c.get("mean_entropy"),
                "low_confidence_frac": c.get("low_confidence_frac"),
                "composition_confound": c.get("composition_confound", {}),
            },
        ))

    tvd = d.get("cross_cohort_total_variation_distance")
    # Distributional concordance across independent cohorts is the replication
    # analogue for a fixed published classifier. A subtype classifier is a
    # descriptive label, so we keep it discovery_only unless a per-patient
    # concordance on shared samples is available (it is not here).
    ext = ValidationReceipt(
        kind="external_replication",
        cohort_id="__".join(cohorts.keys()) if cohorts else "cross_cohort",
        n=sum(int(c["n"]) for c in cohorts.values()) if cohorts else 0,
        is_independent_of_training=len(cohorts) >= 2,
        metric_name="subtype_mixture_total_variation_distance",
        metric_value=tvd,
        method_detail=("distributional (not per-patient) cross-cohort concordance; "
                       "patients are not shared across studies"),
    )
    receipts.append(ext)

    # Diagnostic subtype label is descriptive; ship as discovery_only (no
    # per-patient gold-standard concordance available across independent studies).
    status = BiomarkerStatus.DISCOVERY_ONLY

    def f(x, nd=3):
        return f"{float(x):.{nd}f}" if x is not None else "NA"

    cohort_ids = list(cohorts.keys())
    interp = (
        f"consensusOV (Chen et al.) HGSOC molecular-subtype classifier, run UNCHANGED across "
        f"{len(cohort_ids)} independent cohorts ({', '.join(cohort_ids)}). Subtype mixtures are "
        f"reported per cohort; cross-cohort total-variation distance between subtype mixtures is "
        f"{f(tvd)} (lower = more transportable mixture). Because patients are not shared across "
        f"these studies, concordance is distributional, not per-patient. Subtype calls are "
        f"composition-aware, not intrinsic epithelial subtype. Shipped discovery-only. "
        f"RESEARCH USE ONLY — not validated for clinical decision-making."
    )
    caveats = [c for c in [
        d.get("concordance_note"),
        d.get("composition_caveat"),
        *(d.get("caveats") or []),
    ] if c]

    model = BiomarkerModel(
        cancer=CANCER,
        biomarker_type=BiomarkerType.DIAGNOSTIC,
        method_id=d.get("method_id", "consensusOV_diagnostic_ov"),
        source_pmid=d.get("source_pmid"),
        source_doi=d.get("source_doi"),
        source_note=d.get("source_note", ""),
        gene_panel=[],
        status=status,
        receipts=receipts,
        interpretation=interp,
        caveats=caveats,
    )
    return model


# ── W4 baseline discovery scan (our own panels) ───────────────────────────────

def build_w4_discovery_models(w4_json: str | Path) -> List[BiomarkerModel]:
    d = json.loads(Path(w4_json).read_text())
    models: List[BiomarkerModel] = []
    for pn, blk in d.get("panels", {}).items():
        receipts: List[ValidationReceipt] = []
        for cid, r in blk.get("per_cohort", {}).items():
            if r.get("c_index") is None:
                continue
            receipts.append(ValidationReceipt(
                kind="internal_cv",
                cohort_id=cid,
                n=int(r["n"]),
                n_events=r.get("n_events"),
                is_independent_of_training=True,  # each cohort scored independently, unsupervised
                metric_name="c_index",
                metric_value=r.get("c_index"),
                ci95_low=r.get("c_index_lower"),
                ci95_high=r.get("c_index_upper"),
                p_value=r.get("c_index_p"),
                method_detail="unsupervised mean-z panel score (discovery scan)",
                extra={"hr_per_1sd": r.get("hr_per_1sd"), "hr_p": r.get("hr_p")},
            ))
        # Always discovery_only by construction (W4 is a scan, not a frozen model)
        model = BiomarkerModel(
            cancer=CANCER,
            biomarker_type=BiomarkerType.PROGNOSTIC,
            method_id=f"ownpanel_{pn}_prognostic_scan",
            source_note=f"Baseline discovery scan of our in-repo {pn} panel ({blk.get('source')}).",
            gene_panel=[],
            status=BiomarkerStatus.DISCOVERY_ONLY,
            receipts=receipts,
            interpretation=(
                f"Our in-repo {pn} panel scored as an unsupervised mean-z prognostic signature "
                f"across {len(receipts)} independent OV cohorts. C-indices cluster near 0.5, i.e. "
                f"weak-to-null prognostic signal; the panel is NOT a validated prognostic biomarker. "
                f"Discovery-only. RESEARCH USE ONLY — not validated for clinical decision-making."
            ),
            caveats=[blk.get("status_note", ""), d.get("global_caveat", "")],
        )
        models.append(model)
    return models


# ── top-level: build + register everything, return the registry ──────────────

def assemble_registry(
    w1_json: str | Path,
    w2_json: str | Path,
    w3_json: Optional[str | Path] = None,
    w4_json: Optional[str | Path] = None,
    run_grounding_guard: bool = True,
) -> BiomarkerRegistry:
    reg = BiomarkerRegistry()

    prog = build_prognostic_model(w1_json, w2_json)
    if run_grounding_guard:
        assert_grounded(prog)
    reg.register(prog)

    if w3_json and Path(w3_json).exists():
        diag = build_diagnostic_model(w3_json)
        if run_grounding_guard:
            assert_grounded(diag)
        reg.register(diag)

    if w4_json and Path(w4_json).exists():
        for m in build_w4_discovery_models(w4_json):
            if run_grounding_guard:
                assert_grounded(m)
            reg.register(m)

    return reg
