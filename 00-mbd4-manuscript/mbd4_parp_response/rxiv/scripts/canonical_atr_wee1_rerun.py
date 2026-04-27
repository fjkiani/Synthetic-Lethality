#!/usr/bin/env python3
"""
Single controlled rerun for MBD4 manuscript ATR/WEE1 (ceralasertib / adavosertib) stats.
Outputs one canonical CSV + JSON. Does not modify manuscript prose.

WT rules:
- Methods-aligned: WT = no somatic MBD4 mutation in DepMap OmicsSomaticMutations (any MBD4 call).
- Alternate: WT = all cell lines not in MBD4 True-LOF set (includes MBD4 non-LOF mutants in WT).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[5]
CACHE = REPO / "oncology-coPilot/src/components/orchestrator/Analysis/SyntheticLethality/.cache/depmap"
OUT_DIR = REPO / "publications/00-mbd4-manuscript/mbd4_parp_response/artifacts/canonical_atr_wee1_rerun_20260405"

DEPMAP_RELEASE = "24Q2 (Methods); local parquet bundle under SyntheticLethality/.cache/depmap (no embedded release tag in files)"
GDSC_RELEASE = "GDSC2 (parquet DATASET=GDSC2; Sanger/CancerRxGene GDSC2 screen bundled with DepMap cache)"


def pooled_cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    v1, v2 = np.var(a, ddof=1), np.var(b, ddof=1)
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    if sp == 0 or not np.isfinite(sp):
        return float("nan")
    return float((np.mean(a) - np.mean(b)) / sp)


def mw_less(lof: np.ndarray, wt: np.ndarray) -> Tuple[float, float]:
    """One-sided Mann-Whitney U: LOF < WT (greater sensitivity = lower LN_IC50 etc.)."""
    lof, wt = np.asarray(lof, dtype=float), np.asarray(wt, dtype=float)
    lof, wt = lof[np.isfinite(lof)], wt[np.isfinite(wt)]
    if len(lof) == 0 or len(wt) == 0:
        return float("nan"), float("nan")
    _, p = stats.mannwhitneyu(lof, wt, alternative="less")
    return float(np.mean(lof) - np.mean(wt)), float(p)


def row(
    drug_name: str,
    dataset: str,
    metric: str,
    n_lof: int,
    n_wt: int,
    delta: float,
    p: float,
    d: float,
    cohort_definition: str,
) -> Dict[str, Any]:
    return {
        "drug_name": drug_name,
        "dataset": dataset,
        "metric": metric,
        "n_lof": int(n_lof),
        "n_wt": int(n_wt),
        "delta": delta,
        "p_value_full_precision": p,
        "cohens_d_full_precision": d,
        "cohort_definition": cohort_definition,
        "depmap_release": DEPMAP_RELEASE,
        "gdsc_release": GDSC_RELEASE,
    }


def load_base() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model = pd.read_parquet(CACHE / "Model.parquet").reset_index()
    if "ModelID" not in model.columns and "index" in model.columns:
        model = model.rename(columns={"index": "ModelID"})
    omics = pd.read_parquet(CACHE / "OmicsSomaticMutations.parquet")
    gdsc = pd.read_parquet(CACHE / "GDSC2.parquet")
    gdsc = gdsc.dropna(subset=["COSMIC_ID"])
    gdsc["COSMIC_INT"] = pd.to_numeric(gdsc["COSMIC_ID"], errors="coerce").fillna(0).astype(int)
    model["COSMIC_INT"] = pd.to_numeric(model["COSMICID"], errors="coerce").fillna(0).astype(int)
    model = model[model["COSMIC_INT"] > 0]
    return model, omics, gdsc


def mbd4_sets(omics: pd.DataFrame) -> Tuple[set, set, set]:
    gene_col = "HugoSymbol" if "HugoSymbol" in omics.columns else "Gene"
    mbd4 = omics[omics[gene_col] == "MBD4"]
    if "LikelyLoF" not in mbd4.columns:
        raise RuntimeError("OmicsSomaticMutations missing LikelyLoF")
    lof_models = set(mbd4[mbd4["LikelyLoF"] == True]["ModelID"])
    any_mbd4_models = set(mbd4["ModelID"])
    return lof_models, any_mbd4_models, set(mbd4[mbd4["LikelyLoF"] == True]["ModelID"])


def tp53_mut_models(omics: pd.DataFrame) -> set:
    gene_col = "HugoSymbol" if "HugoSymbol" in omics.columns else "Gene"
    t = omics[omics[gene_col] == "TP53"]
    return set(t["ModelID"])


def cosmic_for_models(model: pd.DataFrame, model_ids: set) -> set:
    return set(model[model["ModelID"].isin(model_ids)]["COSMIC_INT"])


def msi_cosmic(model: pd.DataFrame) -> set:
    col = "ModelSubtypeFeatures"
    if col not in model.columns:
        return set()
    m = model[col].astype(str).str.contains("MSI", case=False, na=False)
    return set(model.loc[m, "COSMIC_INT"])


def bowel_cosmic(model: pd.DataFrame) -> set:
    d = model["OncotreePrimaryDisease"].astype(str)
    m = d.str.contains("Bowel|Colorectal", case=False, na=False)
    return set(model.loc[m, "COSMIC_INT"])


def filter_drug(gdsc: pd.DataFrame, pattern: str) -> pd.DataFrame:
    return gdsc[gdsc["DRUG_NAME"].str.contains(pattern, case=False, na=False)].copy()


def analyze_subset(
    drug_df: pd.DataFrame,
    metric: str,
    lof_cosmic: set,
    wt_cosmic: set,
    row_filter: Optional[Callable[[pd.DataFrame], pd.Series]] = None,
) -> Tuple[int, int, float, float, float, np.ndarray, np.ndarray]:
    d = drug_df
    if row_filter is not None:
        d = d[row_filter(d)]
    is_lof = d["COSMIC_INT"].isin(lof_cosmic)
    is_wt = d["COSMIC_INT"].isin(wt_cosmic)
    lof_v = d.loc[is_lof, metric].dropna().astype(float).values
    wt_v = d.loc[is_wt, metric].dropna().astype(float).values
    delta, p = mw_less(lof_v, wt_v)
    dval = pooled_cohens_d(lof_v, wt_v)
    return len(lof_v), len(wt_v), delta, p, dval, lof_v, wt_v


def main() -> int:
    if not CACHE.is_dir():
        print("RERUN BLOCKED: DepMap cache missing at", CACHE, file=sys.stderr)
        return 2

    model, omics, gdsc = load_base()
    lof_models, any_mbd4_models, _ = mbd4_sets(omics)
    lof_cosmic = cosmic_for_models(model, lof_models)
    all_mut_cosmic = cosmic_for_models(model, any_mbd4_models)
    # Alternate WT: all DepMap models with COSMIC minus LOF cosmic (not manuscript-primary)
    all_cosmic = set(model["COSMIC_INT"])
    wt_alt_cosmic = all_cosmic - lof_cosmic
    # Methods WT: no somatic MBD4
    wt_methods_cosmic = all_cosmic - all_mut_cosmic

    tp53_models = tp53_mut_models(omics)
    tp53_cosmic = cosmic_for_models(model, tp53_models)
    msi_c = msi_cosmic(model)
    bowel_c = bowel_cosmic(model)

    rows: List[Dict[str, Any]] = []

    def add_analysis(
        drug_label: str,
        drug_pattern: str,
        metrics: List[str],
        wt_label: str,
        wt_cosmic: set,
        cohort_suffix: str,
        row_filter: Optional[Callable[[pd.DataFrame], pd.Series]] = None,
    ) -> None:
        sub = filter_drug(gdsc, drug_pattern)
        for metric in metrics:
            if metric not in sub.columns:
                continue
            n_lof, n_wt, delta, p, d, _, _ = analyze_subset(
                sub, metric, lof_cosmic, wt_cosmic, row_filter=row_filter
            )
            rows.append(
                row(
                    drug_label,
                    "GDSC2",
                    metric,
                    n_lof,
                    n_wt,
                    delta,
                    p,
                    d,
                    cohort_suffix + f"; WT_rule={wt_label}",
                )
            )

    # Primary + alternate WT (ceralasertib, adavosertib)
    cera_metrics = ["LN_IC50", "AUC", "Z_SCORE"]
    add_analysis(
        "ceralasertib",
        r"Ceralasertib|AZD6738",
        cera_metrics,
        "Methods_no_somatic_MBD4",
        wt_methods_cosmic,
        "MBD4 True-LOF (LikelyLoF=True) vs WT; primary drug screen rows",
    )
    add_analysis(
        "ceralasertib",
        r"Ceralasertib|AZD6738",
        cera_metrics,
        "Alternate_non_LOF_only_WT_includes_MBD4_nonLoF_mutants",
        wt_alt_cosmic,
        "MBD4 True-LOF vs comparator pool excluding only True-LOF (NOT Methods primary)",
    )

    add_analysis(
        "adavosertib",
        r"MK-1775",
        ["LN_IC50"],
        "Methods_no_somatic_MBD4",
        wt_methods_cosmic,
        "MBD4 True-LOF vs WT; adavosertib = MK-1775 in GDSC2",
    )
    add_analysis(
        "adavosertib",
        r"MK-1775",
        ["LN_IC50"],
        "Alternate_non_LOF_only_WT_includes_MBD4_nonLoF_mutants",
        wt_alt_cosmic,
        "MBD4 True-LOF vs alternate WT pool",
    )

    # MSI purge (ceralasertib LN_IC50 primary WT)
    sub_cera = filter_drug(gdsc, r"Ceralasertib|AZD6738")

    def not_msi(df: pd.DataFrame) -> pd.Series:
        return ~df["COSMIC_INT"].isin(msi_c)

    n_lof, n_wt, delta, p, d, _, _ = analyze_subset(
        sub_cera, "LN_IC50", lof_cosmic, wt_methods_cosmic, row_filter=not_msi
    )
    rows.append(
        row(
            "ceralasertib",
            "GDSC2",
            "LN_IC50",
            n_lof,
            n_wt,
            delta,
            p,
            d,
            "Stress MSI purge: exclude ModelSubtypeFeatures containing 'MSI' from both arms; Methods WT",
        )
    )

    # TP53 stratified (ceralasertib): LOF∩TP53 vs WT∩TP53 (Methods WT)
    lof_tp53 = lof_cosmic & tp53_cosmic
    wt_tp53 = wt_methods_cosmic & tp53_cosmic
    for metric in cera_metrics:
        n_lof, n_wt, delta, p, d, _, _ = analyze_subset(
            sub_cera, metric, lof_tp53, wt_tp53, row_filter=None
        )
        rows.append(
            row(
                "ceralasertib",
                "GDSC2",
                metric,
                n_lof,
                n_wt,
                delta,
                p,
                d,
                "Stress TP53: MBD4-LOF∩TP53-mut vs MBD4-WT(Methods)∩TP53-mut; TP53-mut = any somatic TP53 in OmicsSomaticMutations",
            )
        )

    # Lineage: Bowel / non-Bowel (ceralasertib LN_IC50, Methods WT)

    def is_bowel_df(df: pd.DataFrame) -> pd.Series:
        return df["COSMIC_INT"].isin(bowel_c)

    def is_non_bowel_df(df: pd.DataFrame) -> pd.Series:
        return ~df["COSMIC_INT"].isin(bowel_c)

    n_lof, n_wt, delta, p, d, _, _ = analyze_subset(
        sub_cera, "LN_IC50", lof_cosmic, wt_methods_cosmic, row_filter=is_bowel_df
    )
    rows.append(
        row(
            "ceralasertib",
            "GDSC2",
            "LN_IC50",
            n_lof,
            n_wt,
            delta,
            p,
            d,
            "Stress Lineage Bowel: OncotreePrimaryDisease matches Bowel|Colorectal; Methods WT",
        )
    )

    n_lof, n_wt, delta, p, d, _, _ = analyze_subset(
        sub_cera, "LN_IC50", lof_cosmic, wt_methods_cosmic, row_filter=is_non_bowel_df
    )
    rows.append(
        row(
            "ceralasertib",
            "GDSC2",
            "LN_IC50",
            n_lof,
            n_wt,
            delta,
            p,
            d,
            "Stress Lineage non-Bowel: complement of Bowel|Colorectal filter; Methods WT",
        )
    )

    # Leave-one-out on LOF lines with ceralasertib LN_IC50 (Methods WT)
    is_lof_base = sub_cera["COSMIC_INT"].isin(lof_cosmic)
    is_wt_base = sub_cera["COSMIC_INT"].isin(wt_methods_cosmic)
    wt_vals_full = sub_cera.loc[is_wt_base, "LN_IC50"].dropna().astype(float).values
    lof_lines = sorted(
        sub_cera.loc[is_lof_base & sub_cera["LN_IC50"].notna(), "COSMIC_INT"].unique().tolist()
    )
    ps: List[float] = []
    for drop in lof_lines:
        lof_v = sub_cera.loc[
            is_lof_base & (sub_cera["COSMIC_INT"] != drop), "LN_IC50"
        ].dropna().astype(float).values
        _, p_i = mw_less(lof_v, wt_vals_full)
        if np.isfinite(p_i):
            ps.append(p_i)
    if ps:
        rows.append(
            row(
                "ceralasertib",
                "GDSC2",
                "LN_IC50_LOO_max_p",
                len(lof_lines),
                len(wt_vals_full),
                float(np.nan),
                float(max(ps)),
                float(np.nan),
                f"Leave-one-out over {len(lof_lines)} LOF cosmids with LN_IC50; reported p=max iteration p (one-sided MWU LOF<WT); all_lt_0.05={all(x < 0.05 for x in ps)}",
            )
        )
        rows.append(
            row(
                "ceralasertib",
                "GDSC2",
                "LN_IC50_LOO_min_p",
                len(lof_lines),
                len(wt_vals_full),
                float(np.nan),
                float(min(ps)),
                float(np.nan),
                "Leave-one-out summary: min iteration p",
            )
        )

    def _json_safe(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {k: _json_safe(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_json_safe(v) for v in obj]
        if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return None
        return obj

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "canonical_atr_wee1_rerun.csv"
    json_path = OUT_DIR / "canonical_atr_wee1_rerun.json"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(_json_safe(rows), f, indent=2)

    print(str(csv_path))
    print(str(json_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
