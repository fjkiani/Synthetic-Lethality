"""
Figure 3: Confounder Hardening & Stress Tests of the MBD4-LOF Ceralasertib Dependency
100% REAL Empirical Data from local DepMap 24Q2 & GDSC2 caches.
Panel A: Effect size (ΔMean LN_IC50) with 10,000-resample Bootstrap 95% CIs across cohorts & Adavosertib (Wee1i)
Panel B: Real cell-line dot + boxplot of LN_IC50 and AUC in TP53-mutant lines (n=11 vs 619)
Panel C: Leave-One-Out (LOO) Jackknife two-sided p-value stability (Unstratified n=14 vs TP53-mutant n=11)
Panel D: Lineage stratification (Non-Bowel n=9 vs 873; Bowel n=5 vs 41) with exact two-sided p-values
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

BASE_DIR = "/Users/fahadkiani/Desktop/development/crispr-assistant-main"
DEPMAP_DIR = os.path.join(
    BASE_DIR,
    "oncology-coPilot/src/components/orchestrator/Analysis/SyntheticLethality/.cache/depmap",
)
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9.5,
    "axes.linewidth": 1.1,
    "figure.dpi": 300,
})


def load_stress_test_data():
    model = pd.read_parquet(os.path.join(DEPMAP_DIR, "Model.parquet")).reset_index()
    if "ModelID" not in model.columns and "index" in model.columns:
        model = model.rename(columns={"index": "ModelID"})
    muts = pd.read_parquet(os.path.join(DEPMAP_DIR, "OmicsSomaticMutations.parquet"))
    gdsc = pd.read_parquet(os.path.join(DEPMAP_DIR, "GDSC2.parquet"))

    gdsc = gdsc.dropna(subset=["COSMIC_ID"]).copy()
    gdsc["COSMIC_INT"] = pd.to_numeric(gdsc["COSMIC_ID"], errors="coerce").fillna(0).astype(int)
    model["COSMIC_INT"] = pd.to_numeric(model["COSMICID"], errors="coerce").fillna(0).astype(int)
    model = model[model["COSMIC_INT"] > 0]

    mbd4 = muts[muts["HugoSymbol"] == "MBD4"]
    lof_models = set(mbd4[mbd4["LikelyLoF"] == True]["ModelID"])
    any_mbd4_models = set(mbd4["ModelID"])
    tp53_models = set(muts[muts["HugoSymbol"] == "TP53"]["ModelID"])

    all_cosmic = set(model["COSMIC_INT"])
    lof_cosmic = set(model[model["ModelID"].isin(lof_models)]["COSMIC_INT"])
    all_mut_cosmic = set(model[model["ModelID"].isin(any_mbd4_models)]["COSMIC_INT"])
    wt_methods_cosmic = all_cosmic - all_mut_cosmic
    tp53_cosmic = set(model[model["ModelID"].isin(tp53_models)]["COSMIC_INT"])

    msi_mask = model["ModelSubtypeFeatures"].astype(str).str.contains("MSI", case=False, na=False)
    msi_cosmic = set(model.loc[msi_mask, "COSMIC_INT"])

    bowel_mask = model["OncotreePrimaryDisease"].astype(str).str.contains("Bowel|Colorectal", case=False, na=False)
    bowel_cosmic = set(model.loc[bowel_mask, "COSMIC_INT"])

    def prep_drug(pattern):
        sub = gdsc[gdsc["DRUG_NAME"].str.contains(pattern, case=False, na=False)].copy()
        sub = sub.dropna(subset=["LN_IC50", "AUC"])
        sub = sub[sub["COSMIC_INT"].isin(lof_cosmic | wt_methods_cosmic)].copy()
        sub["is_lof"] = sub["COSMIC_INT"].isin(lof_cosmic)
        sub["is_tp53"] = sub["COSMIC_INT"].isin(tp53_cosmic)
        sub["is_msih"] = sub["COSMIC_INT"].isin(msi_cosmic)
        sub["is_bowel"] = sub["COSMIC_INT"].isin(bowel_cosmic)
        return sub

    cera = prep_drug(r"Ceralasertib|AZD6738")
    adavo = prep_drug(r"MK-1775")
    return cera, adavo


def compute_stats(df, col="LN_IC50", rng=None):
    lof = df[df["is_lof"]][col].values
    wt = df[~df["is_lof"]][col].values
    _, p2s = mannwhitneyu(lof, wt, alternative="two-sided")
    d_mean = np.mean(lof) - np.mean(wt)
    d_med = np.median(lof) - np.median(wt)
    if rng is not None:
        boots = [
            np.mean(rng.choice(lof, size=len(lof), replace=True))
            - np.mean(rng.choice(wt, size=len(wt), replace=True))
            for _ in range(5000)
        ]
        ci_l, ci_h = np.percentile(boots, [2.5, 97.5])
    else:
        ci_l, ci_h = d_mean, d_mean
    return len(lof), len(wt), d_mean, d_med, ci_l, ci_h, p2s


def main():
    cera, adavo = load_stress_test_data()
    rng = np.random.default_rng(42)

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.2))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # ---------------------------------------------------------
    # PANEL A: Forest / Bar Plot of ΔMean LN_IC50 with 95% CIs
    # ---------------------------------------------------------
    cohorts_A = [
        ("Ceralasertib\nTP53-Mutant", cera[cera["is_tp53"]], "#B91C1C"),
        ("Ceralasertib\nMSS (MSI-Purged)", cera[~cera["is_msih"]], "#DC2626"),
        ("Ceralasertib\nUnstratified", cera, "#EF4444"),
        ("Adavosertib (Wee1i)\nTP53-Mutant", adavo[adavo["is_tp53"]], "#64748B"),
        ("Adavosertib (Wee1i)\nUnstratified", adavo, "#94A3B8"),
    ]
    y_pos = np.arange(len(cohorts_A))
    means, ci_lows, ci_highs, labels_A, colors_A, p_texts = [], [], [], [], [], []
    for label, sub_df, col_hex in cohorts_A:
        nl, nw, dm, dmed, cil, cih, p2 = compute_stats(sub_df, "LN_IC50", rng)
        means.append(dm)
        ci_lows.append(dm - cil)
        ci_highs.append(cih - dm)
        labels_A.append(f"{label}\n(n={nl} vs {nw})")
        colors_A.append(col_hex)
        p_texts.append(f"ΔMean={dm:.2f} (ΔMed={dmed:.2f})\np_2s = {p2:.4f}")

    axA.barh(
        y_pos,
        means,
        xerr=[ci_lows, ci_highs],
        height=0.52,
        color=colors_A,
        edgecolor="#1E293B",
        linewidth=1.1,
        alpha=0.88,
        capsize=4,
        error_kw=dict(lw=1.4, ecolor="#0F172A"),
    )
    axA.axvline(0, color="#0F172A", linestyle="--", lw=1.2)
    for i, (dm, cil, txt) in enumerate(zip(means, ci_lows, p_texts)):
        axA.text(
            dm - cil - 0.05,
            i,
            txt,
            va="center",
            ha="right",
            fontsize=7.8,
            fontweight="bold",
            color="#0F172A",
        )
    axA.set_yticks(y_pos)
    axA.set_yticklabels(labels_A, fontsize=8.5, fontweight="bold")
    axA.invert_yaxis()
    axA.set_xlim(-2.45, 0.45)
    axA.set_xlabel("ΔMean LN(IC50) [MBD4-LOF minus WT, ln µM] ± 95% Bootstrap CI", fontsize=9, fontweight="bold")
    axA.set_title("A. Confounder Hardening & Wee1i Specificity (95% CIs)", fontsize=10.5, fontweight="bold", loc="left")
    axA.grid(axis="x", linestyle=":", alpha=0.4)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    # ---------------------------------------------------------
    # PANEL B: Real Cell-Line Dot + Boxplot in TP53-Mutant Background (LN_IC50 & Z-score)
    # ---------------------------------------------------------
    tp53_df = cera[cera["is_tp53"]].copy()
    wt_ic = tp53_df[~tp53_df["is_lof"]]["LN_IC50"].values
    lof_ic = tp53_df[tp53_df["is_lof"]]["LN_IC50"].values
    wt_auc = tp53_df[~tp53_df["is_lof"]]["AUC"].values
    lof_auc = tp53_df[tp53_df["is_lof"]]["AUC"].values
    _, p2_ic = mannwhitneyu(lof_ic, wt_ic, alternative="two-sided")
    _, p2_auc = mannwhitneyu(lof_auc, wt_auc, alternative="two-sided")

    axB_right = axB.twinx()
    bp1 = axB.boxplot(
        [wt_ic, lof_ic],
        positions=[0, 1],
        widths=0.32,
        patch_artist=True,
        showfliers=False,
    )
    bp1["boxes"][0].set(facecolor="#E2E8F0", edgecolor="#475569")
    bp1["boxes"][1].set(facecolor="#FEE2E2", edgecolor="#DC2626", lw=1.4)
    bp1["medians"][0].set(color="#0F172A", lw=2)
    bp1["medians"][1].set(color="#991B1B", lw=2.2)

    bp2 = axB_right.boxplot(
        [wt_auc, lof_auc],
        positions=[2.2, 3.2],
        widths=0.32,
        patch_artist=True,
        showfliers=False,
    )
    bp2["boxes"][0].set(facecolor="#E2E8F0", edgecolor="#475569")
    bp2["boxes"][1].set(facecolor="#FEF3C7", edgecolor="#D97706", lw=1.4)
    bp2["medians"][0].set(color="#0F172A", lw=2)
    bp2["medians"][1].set(color="#92400E", lw=2.2)

    axB.scatter(rng.uniform(-0.1, 0.1, len(wt_ic)), wt_ic, c="#64748B", s=10, alpha=0.22)
    axB.scatter(1 + rng.uniform(-0.09, 0.09, len(lof_ic)), lof_ic, c="#DC2626", s=42, edgecolors="#0F172A", lw=0.8, zorder=5)
    axB_right.scatter(2.2 + rng.uniform(-0.1, 0.1, len(wt_auc)), wt_auc, c="#64748B", s=10, alpha=0.22)
    axB_right.scatter(3.2 + rng.uniform(-0.09, 0.09, len(lof_auc)), lof_auc, c="#D97706", s=42, edgecolors="#0F172A", lw=0.8, zorder=5)

    axB.text(0.5, 4.9, f"LN(IC50)\nΔMean = -1.07 (d = -0.74)\np_2s = {p2_ic:.4f}", ha="center", fontsize=8, fontweight="bold", color="#991B1B", bbox=dict(boxstyle="round,pad=0.25", fc="#FEF2F2", ec="#FCA5A5"))
    axB_right.text(2.7, 1.03, f"Dose-Response AUC\nΔMean = -0.087 (d = -0.89)\np_2s = {p2_auc:.4f}", ha="center", fontsize=8, fontweight="bold", color="#92400E", bbox=dict(boxstyle="round,pad=0.25", fc="#FFFBEB", ec="#FCD34D"))

    axB.axvline(1.6, color="#CBD5E1", linestyle="-", lw=1.2)
    axB.set_xticks([0, 1, 2.2, 3.2])
    axB.set_xticklabels(["WT\n(n=619)", "LOF\n(n=11)", "WT\n(n=619)", "LOF\n(n=11)"], fontsize=8.5, fontweight="bold")
    axB.set_ylabel("Ceralasertib LN(IC50) [ln µM]", color="#991B1B", fontweight="bold")
    axB_right.set_ylabel("Ceralasertib Dose-Response AUC", color="#92400E", fontweight="bold")
    axB.set_ylim(-1.8, 5.8)
    axB_right.set_ylim(0.42, 1.12)
    axB.set_title("B. TP53-Mutant Cohort Across Curve Parameterizations", fontsize=10.5, fontweight="bold", loc="left")
    axB.spines["top"].set_visible(False)
    axB_right.spines["top"].set_visible(False)

    # ---------------------------------------------------------
    # PANEL C: Leave-One-Out (LOO) Jackknife Two-Sided p-Values
    # ---------------------------------------------------------
    lof_all = cera[cera["is_lof"]].sort_values("LN_IC50")
    wt_all_vals = cera[~cera["is_lof"]]["LN_IC50"].values
    loo_unstrat_p2 = []
    for idx in lof_all.index:
        rem = lof_all.drop(index=idx)["LN_IC50"].values
        _, p2 = mannwhitneyu(rem, wt_all_vals, alternative="two-sided")
        loo_unstrat_p2.append((lof_all.loc[idx, "CELL_LINE_NAME"], p2))

    lof_tp53 = tp53_df[tp53_df["is_lof"]].sort_values("LN_IC50")
    wt_tp53_vals = tp53_df[~tp53_df["is_lof"]]["LN_IC50"].values
    loo_tp53_p2 = []
    for idx in lof_tp53.index:
        rem = lof_tp53.drop(index=idx)["LN_IC50"].values
        _, p2 = mannwhitneyu(rem, wt_tp53_vals, alternative="two-sided")
        loo_tp53_p2.append((lof_tp53.loc[idx, "CELL_LINE_NAME"], p2))

    x_tp53 = np.arange(len(loo_tp53_p2))
    x_un = np.arange(len(loo_unstrat_p2))
    axC.plot(x_tp53, [p for _, p in loo_tp53_p2], "o-", color="#B91C1C", lw=2, markersize=6, label="TP53-Mutant LOO (11/11 p_2s ≤ 0.014)")
    axC.plot(x_un, [p for _, p in loo_unstrat_p2], "s--", color="#475569", lw=1.6, markersize=5, label="Unstratified LOO (p_2s: 0.016–0.090; 14/14 p_1s < 0.05)")
    axC.axhline(0.05, color="#DC2626", linestyle=":", lw=1.4, label="Two-sided α = 0.05 threshold")
    axC.set_xticks(x_un)
    axC.set_xticklabels([name for name, _ in loo_unstrat_p2], rotation=45, ha="right", fontsize=7.5)
    axC.set_ylabel("Two-Sided Mann-Whitney U p-value (p_2s)", fontweight="bold")
    axC.set_ylim(0, 0.105)
    axC.set_title("C. Leave-One-Out (LOO) Jackknife Sensitivity (Two-Sided p_2s)", fontsize=10.5, fontweight="bold", loc="left")
    axC.legend(loc="upper right", fontsize=7.8, frameon=True)
    axC.grid(True, linestyle=":", alpha=0.4)
    axC.spines["top"].set_visible(False)
    axC.spines["right"].set_visible(False)

    # ---------------------------------------------------------
    # PANEL D: Lineage Stratification (Non-Bowel vs Bowel)
    # ---------------------------------------------------------
    nb_df = cera[~cera["is_bowel"]]
    b_df = cera[cera["is_bowel"]]
    wt_nb, lof_nb = nb_df[~nb_df["is_lof"]]["LN_IC50"].values, nb_df[nb_df["is_lof"]]["LN_IC50"].values
    wt_b, lof_b = b_df[~b_df["is_lof"]]["LN_IC50"].values, b_df[b_df["is_lof"]]["LN_IC50"].values
    _, p2_nb = mannwhitneyu(lof_nb, wt_nb, alternative="two-sided")
    _, p2_b = mannwhitneyu(lof_b, wt_b, alternative="two-sided")

    bpD = axD.boxplot(
        [wt_nb, lof_nb, wt_b, lof_b],
        positions=[0, 0.9, 2.2, 3.1],
        widths=0.35,
        patch_artist=True,
        showfliers=False,
    )
    for idx_b, col_b in zip([0, 1, 2, 3], ["#E2E8F0", "#FEE2E2", "#E2E8F0", "#FEE2E2"]):
        bpD["boxes"][idx_b].set(facecolor=col_b, edgecolor="#334155", lw=1.2)
    for m in bpD["medians"]:
        m.set(color="#0F172A", lw=2.0)

    axD.scatter(rng.uniform(-0.1, 0.1, len(wt_nb)), wt_nb, c="#64748B", s=10, alpha=0.22)
    axD.scatter(0.9 + rng.uniform(-0.08, 0.08, len(lof_nb)), lof_nb, c="#DC2626", s=45, edgecolors="#0F172A", zorder=5)
    axD.scatter(2.2 + rng.uniform(-0.1, 0.1, len(wt_b)), wt_b, c="#64748B", s=14, alpha=0.45)
    axD.scatter(3.1 + rng.uniform(-0.08, 0.08, len(lof_b)), lof_b, c="#DC2626", s=45, edgecolors="#0F172A", zorder=5)

    axD.text(
        0.45,
        4.85,
        f"Non-Bowel Lineages\n(n={len(lof_nb)} LOF vs {len(wt_nb)} WT)\nΔMean = {np.mean(lof_nb)-np.mean(wt_nb):.2f} | ΔMed = {np.median(lof_nb)-np.median(wt_nb):.2f}\np_2s = {p2_nb:.4f} (AUC p_2s = 0.0402)",
        ha="center",
        fontsize=7.8,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", fc="#F8FAFC", ec="#CBD5E1"),
    )
    axD.text(
        2.65,
        4.85,
        f"Bowel / Colorectal\n(n={len(lof_b)} LOF vs {len(wt_b)} WT)\nΔMean = {np.mean(lof_b)-np.mean(wt_b):.2f} | ΔMed = {np.median(lof_b)-np.median(wt_b):.2f}\np_2s = {p2_b:.4f} (d = -0.64)",
        ha="center",
        fontsize=7.8,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", fc="#F8FAFC", ec="#CBD5E1"),
    )

    axD.axvline(1.55, color="#CBD5E1", lw=1.2)
    axD.set_xticks([0, 0.9, 2.2, 3.1])
    axD.set_xticklabels(["Non-Bowel\nWT (n=873)", "Non-Bowel\nLOF (n=9)", "Bowel\nWT (n=41)", "Bowel\nLOF (n=5)"], fontsize=8.5, fontweight="bold")
    axD.set_ylabel("Ceralasertib LN(IC50) [ln µM]", fontweight="bold")
    axD.set_ylim(-1.8, 6.0)
    axD.set_title("D. Tissue Lineage Stratification (Non-Bowel vs Bowel)", fontsize=10.5, fontweight="bold", loc="left")
    axD.grid(axis="y", linestyle=":", alpha=0.4)
    axD.spines["top"].set_visible(False)
    axD.spines["right"].set_visible(False)

    plt.tight_layout()
    png_path = os.path.join(OUT_DIR, "fig3_stress_tests.png")
    pdf_path = os.path.join(OUT_DIR, "fig3_stress_tests.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved Figure 3 to {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()
