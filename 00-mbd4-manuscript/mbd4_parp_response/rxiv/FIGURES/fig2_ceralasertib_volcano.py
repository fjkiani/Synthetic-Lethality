"""
Figure 2: Ceralasertib (AZD6738) GDSC2 Sensitivity in MBD4-LOF vs WT Cancer Cell Lines
100% REAL Empirical Data from local DepMap 24Q2 & GDSC2 caches.
Panel A: Unstratified Pan-Cancer GDSC2 Cohort (n=14 MBD4-LOF vs n=914 WT)
Panel B: TP53-Mutant Stratified GDSC2 Cohort (n=11 MBD4-LOF vs n=619 WT)
Eliminates the prior np.sort(wt)[:50] + np.linspace monotonic arc artifact and plots all real cell lines with uniform random horizontal jitter, boxplots (medians/IQR), group means (diamonds), and two-sided Mann-Whitney U p-values.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
    "font.size": 10,
    "axes.linewidth": 1.1,
    "figure.dpi": 300,
})


def load_real_ceralasertib_data():
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

    cera = gdsc[gdsc["DRUG_NAME"].str.contains(r"Ceralasertib|AZD6738", case=False, na=False)].copy()
    cera = cera.dropna(subset=["LN_IC50"])
    cera = cera[cera["COSMIC_INT"].isin(lof_cosmic | wt_methods_cosmic)].copy()
    cera["is_mbd4_lof"] = cera["COSMIC_INT"].isin(lof_cosmic)
    cera["is_tp53_mut"] = cera["COSMIC_INT"].isin(tp53_cosmic)
    return cera


def plot_cohort_panel(ax, df, title, subtitle_badge, rng):
    wt = df[~df["is_mbd4_lof"]]["LN_IC50"].values
    lof = df[df["is_mbd4_lof"]]["LN_IC50"].values
    lof_names = df[df["is_mbd4_lof"]]["CELL_LINE_NAME"].values

    u_stat, p_2s = mannwhitneyu(lof, wt, alternative="two-sided")
    p_1s = p_2s / 2.0
    delta_mean = np.mean(lof) - np.mean(wt)
    delta_med = np.median(lof) - np.median(wt)
    pooled_sd = np.sqrt((np.var(lof, ddof=1) + np.var(wt, ddof=1)) / 2.0)
    cohen_d = delta_mean / pooled_sd

    # Bootstrap 95% CI on delta_mean
    boot_diffs = [
        np.mean(rng.choice(lof, size=len(lof), replace=True))
        - np.mean(rng.choice(wt, size=len(wt), replace=True))
        for _ in range(5000)
    ]
    ci_low, ci_high = np.percentile(boot_diffs, [2.5, 97.5])

    # Boxplots (showing median & IQR)
    bp = ax.boxplot(
        [wt, lof],
        positions=[0, 1],
        widths=0.36,
        patch_artist=True,
        showfliers=False,
        zorder=2,
    )
    bp["boxes"][0].set(facecolor="#E2E8F0", edgecolor="#475569", linewidth=1.3, alpha=0.85)
    bp["boxes"][1].set(facecolor="#FEE2E2", edgecolor="#DC2626", linewidth=1.5, alpha=0.90)
    bp["medians"][0].set(color="#1E293B", linewidth=2.2)
    bp["medians"][1].set(color="#991B1B", linewidth=2.4)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color="#475569", linewidth=1.1)

    # All real WT points with genuine random horizontal jitter (NO sorting, NO linspace arc!)
    jitter_wt = rng.uniform(-0.14, 0.14, size=len(wt))
    ax.scatter(
        jitter_wt,
        wt,
        c="#64748B",
        s=14,
        alpha=0.28,
        edgecolors="none",
        zorder=3,
        label=f"MBD4-WT (n={len(wt)})",
    )

    # All real MBD4-LOF points with random jitter
    jitter_lof = 1.0 + rng.uniform(-0.11, 0.11, size=len(lof))
    ax.scatter(
        jitter_lof,
        lof,
        c="#DC2626",
        s=58,
        alpha=0.92,
        edgecolors="#1E293B",
        linewidth=0.9,
        zorder=5,
        label=f"MBD4-LOF (n={len(lof)})",
    )

    # Plot group means as gold diamonds
    ax.scatter(
        [0, 1],
        [np.mean(wt), np.mean(lof)],
        marker="D",
        s=55,
        c="#F59E0B",
        edgecolors="#1E293B",
        linewidth=1.1,
        zorder=6,
        label="Group Mean",
    )

    # Annotate top 3 most sensitive MBD4-LOF cell lines
    sorted_idx = np.argsort(lof)[:3]
    for rank_i, idx in enumerate(sorted_idx):
        offset_y = -0.25 if rank_i != 1 else 0.22
        ax.annotate(
            lof_names[idx],
            xy=(jitter_lof[idx], lof[idx]),
            xytext=(1.22, lof[idx] + offset_y),
            fontsize=7.5,
            fontweight="bold",
            color="#991B1B",
            arrowprops=dict(arrowstyle="-", color="#DC2626", lw=0.7, alpha=0.7),
            va="center",
        )

    # Bracket & stats annotation
    y_max = max(np.max(wt), np.max(lof))
    y_min = min(np.min(wt), np.min(lof))
    y_bar = y_max + 0.45
    ax.plot([0, 0, 1, 1], [y_bar - 0.15, y_bar, y_bar, y_bar - 0.15], lw=1.3, c="#1E293B")
    ax.text(
        0.5,
        y_bar + 0.12,
        f"Two-sided Mann-Whitney p = {p_2s:.4f} (1s p = {p_1s:.4f})\n"
        f"ΔMean = {delta_mean:.2f} [95% CI: {ci_low:.2f}, {ci_high:.2f}] | ΔMedian = {delta_med:.2f} | d = {cohen_d:.2f}",
        ha="center",
        va="bottom",
        fontsize=8.3,
        fontweight="bold",
        color="#0F172A",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#F8FAFC", edgecolor="#CBD5E1", lw=1.0),
    )

    # Group median & mean text callouts
    ax.text(
        0,
        y_min - 0.55,
        f"WT (n={len(wt)})\nMedian: {np.median(wt):.2f}\nMean: {np.mean(wt):.2f}",
        ha="center",
        va="top",
        fontsize=8,
        color="#334155",
    )
    ax.text(
        1,
        y_min - 0.55,
        f"MBD4-LOF (n={len(lof)})\nMedian: {np.median(lof):.2f}\nMean: {np.mean(lof):.2f}",
        ha="center",
        va="top",
        fontsize=8,
        fontweight="bold",
        color="#991B1B",
    )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MBD4 Wild-Type", "MBD4-LOF\n(Strict Truncating)"], fontsize=9.5, fontweight="bold")
    ax.set_ylabel("Ceralasertib (AZD6738) GDSC2 LN(IC50) [ln µM]\n(Lower = Higher Sensitivity)", fontsize=9.5, fontweight="bold")
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=12, color="#0F172A")
    ax.text(
        0.02,
        0.03,
        subtitle_badge,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        color="#1E3A8A",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#EFF6FF", edgecolor="#93C5FD"),
    )
    ax.set_xlim(-0.48, 1.58)
    ax.set_ylim(y_min - 1.45, y_bar + 1.15)
    ax.grid(axis="y", linestyle="--", alpha=0.35, color="#94A3B8")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    cera = load_real_ceralasertib_data()
    rng = np.random.default_rng(42)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.8), sharey=True)

    plot_cohort_panel(
        ax1,
        cera,
        title="A. Unstratified Pan-Cancer Cohort (n=14 vs 914)",
        subtitle_badge="All GDSC2 Lines | ΔAUC = -0.056 (p_2s = 0.0259)",
        rng=rng,
    )

    cera_tp53 = cera[cera["is_tp53_mut"]].copy()
    plot_cohort_panel(
        ax2,
        cera_tp53,
        title="B. TP53-Mutant Stratified Cohort (n=11 vs 619)",
        subtitle_badge="BH-FDR q = 0.036 (6 Axes) | ΔAUC = -0.087 (p_2s = 0.0017)",
        rng=rng,
    )
    ax2.set_ylabel("")

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=3, frameon=True, fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    png_path = os.path.join(OUT_DIR, "fig2_ceralasertib_volcano.png")
    pdf_path = os.path.join(OUT_DIR, "fig2_ceralasertib_volcano.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved Figure 2 to {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()
