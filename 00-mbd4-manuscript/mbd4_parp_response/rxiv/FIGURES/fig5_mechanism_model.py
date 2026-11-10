"""
Figure 5: Mechanistic Model of Upstream MBD4 Glycosylase Loss vs Downstream BER Trapping & ATR Dependency
Addresses Reviewer 1 & Reviewer 2 Item 5 critiques:
- Replaces "U:G Mismatches" with "T:G Mismatches at Deaminated 5mCpG Sites"
- Contrasts Upstream MBD4 Glycosylase Loss (no nicked SSB / no PARP1 trapping) vs Downstream XRCC1/PolB Block
- Removes the skull glyph and legacy product watermark from the figure export
- Increases font size and legibility across all nodes
- Reports exact two-sided Mann-Whitney U p-values (p_2s = 0.0060 TP53-stratified; p_2s = 0.0430 unstratified; PARPi p_2s = 0.841; TOST p = 0.0199)
"""
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "figure.dpi": 300,
})


def draw_box(ax, x, y, w, h, title, body, fc, ec, title_color="#0F172A", body_color="#1E293B", lw=1.8):
    box = mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        x,
        y + h * 0.20,
        title,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color=title_color,
        zorder=4,
    )
    ax.text(
        x,
        y - h * 0.16,
        body,
        ha="center",
        va="center",
        fontsize=9.5,
        color=body_color,
        linespacing=1.32,
        zorder=4,
    )


def draw_arrow(ax, x1, y1, x2, y2, color="#334155", lw=2.0, style="-|>"):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=16),
        zorder=2,
    )


def main():
    fig, ax = plt.subplots(figsize=(13.5, 7.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title Banner
    ax.text(
        0.5,
        0.96,
        "Biochemical & Pharmacogenomic Divergence in MBD4-LOF Base Excision Repair Deficiency",
        ha="center",
        va="center",
        fontsize=13.5,
        fontweight="bold",
        color="#0F172A",
    )

    # Top Root Node: Upstream MBD4 Glycosylase Inactivation
    draw_box(
        ax,
        x=0.50,
        y=0.83,
        w=0.68,
        h=0.14,
        title="UPSTREAM GLYCOSYLASE LOSS: Biallelic MBD4 Inactivation (Truncating LOF, NMD p_2s = 0.0367)",
        body="Failure to excise thymine from T:G mismatches at deaminated 5mCpG dinucleotides\n"
             "and failure to excise misincorporated 5-substituted pyrimidine analogs (Gemcitabine / 5-FU)",
        fc="#FEF2F2",
        ec="#DC2626",
        title_color="#991B1B",
        lw=2.2,
    )

    # Branch Arrows from Root to Three Columns
    draw_arrow(ax, 0.30, 0.76, 0.17, 0.64, color="#1D4ED8", lw=2.2)
    draw_arrow(ax, 0.50, 0.76, 0.50, 0.64, color="#B91C1C", lw=2.2)
    draw_arrow(ax, 0.70, 0.76, 0.83, 0.64, color="#475569", lw=2.2)

    # ---------------------------------------------------------
    # COLUMN 1 (LEFT): AXIS 1 — CYTIDINE ANALOG BOTTLENECK
    # ---------------------------------------------------------
    draw_box(
        ax,
        x=0.17,
        y=0.53,
        w=0.29,
        h=0.19,
        title="AXIS 1: Pyrimidine Clearance Defect",
        body="• Intact MBD4 excises dFdC / 5-FU\n"
             "• MBD4-LOF leaves analogs & T:G\n"
             "  mismatches embedded in template\n"
             "• Polymerase δ/ε stalling at fork",
        fc="#EFF6FF",
        ec="#2563EB",
        title_color="#1E40AF",
    )
    draw_arrow(ax, 0.17, 0.43, 0.17, 0.31, color="#2563EB", lw=2.2)
    draw_box(
        ax,
        x=0.17,
        y=0.19,
        w=0.29,
        h=0.20,
        title="ACTIONABLE: Gemcitabine (dFdC)",
        body="• SL-XGBRanker Zero-Shot Rank #1\n"
             "• Loads replication-fork lesions\n"
             "• Sequential staggered dosing with\n"
             "  ATRi (ATRiUM / TRESR precedent)",
        fc="#DBEAFE",
        ec="#1D4ED8",
        title_color="#1E3A8A",
        lw=2.2,
    )

    # Cross-link arrow from Axis 1 (Fork Stalling) to Axis 2 (RPA-ssDNA / ATR)
    draw_arrow(ax, 0.315, 0.53, 0.355, 0.53, color="#991B1B", lw=2.0)

    # ---------------------------------------------------------
    # COLUMN 2 (CENTER): AXIS 2 — REPLICATION STRESS & ATR DEPENDENCY
    # ---------------------------------------------------------
    draw_box(
        ax,
        x=0.50,
        y=0.53,
        w=0.29,
        h=0.19,
        title="AXIS 2: Replication-Fork Stress",
        body="• Persistent T:G & dFdC lesions uncouple\n"
             "  helicase from replicative polymerase\n"
             "• Accumulation of RPA-coated ssDNA\n"
             "• TP53 loss removes G1/S arrest",
        fc="#FEF2F2",
        ec="#DC2626",
        title_color="#991B1B",
    )
    draw_arrow(ax, 0.50, 0.43, 0.50, 0.31, color="#DC2626", lw=2.2)
    draw_box(
        ax,
        x=0.50,
        y=0.19,
        w=0.29,
        h=0.20,
        title="SYNTHETIC LETHAL: Ceralasertib (ATRi)",
        body="• TP53-Mutant: ΔLN_IC50 = -1.07 (d = -0.74)\n"
             "  p_2s = 0.0060 (BH-FDR q = 0.036)\n"
             "• Unstratified: ΔLN_IC50 = -0.73, p_2s = 0.0430\n"
             "• MSS-Only: ΔLN_IC50 = -0.91, p_2s = 0.0307",
        fc="#FEE2E2",
        ec="#B91C1C",
        title_color="#7F1D1D",
        lw=2.2,
    )

    # ---------------------------------------------------------
    # COLUMN 3 (RIGHT): AXIS 3 — ABSENCE OF SSB & PARPi DECOUPLING
    # ---------------------------------------------------------
    draw_box(
        ax,
        x=0.83,
        y=0.53,
        w=0.29,
        h=0.19,
        title="AXIS 3: No Nicked SSB Substrate",
        body="• Unlike XRCC1/Polβ loss (downstream),\n"
             "  upstream MBD4 loss prevents APE1\n"
             "  backbone incision (no 5'-dRP nick)\n"
             "• Zero PARP1 upregulation (TOST p = 0.020)",
        fc="#F8FAFC",
        ec="#64748B",
        title_color="#334155",
    )
    draw_arrow(ax, 0.83, 0.43, 0.83, 0.31, color="#64748B", lw=2.2)
    draw_box(
        ax,
        x=0.83,
        y=0.19,
        w=0.29,
        h=0.20,
        title="CONFIRMED NULL: Single-Agent PARPi",
        body="• No PARP1-DNA trapping substrate\n"
             "• Direct 5-Drug GDSC2 Evaluation (n=8 vs 478):\n"
             "  Talazoparib (p=0.97), Olaparib (p=0.52),\n"
             "  Composite PARPi Z-score p_2s = 0.841",
        fc="#F1F5F9",
        ec="#475569",
        title_color="#1E293B",
        lw=2.2,
    )

    plt.tight_layout()
    png_path = os.path.join(OUT_DIR, "fig5_mechanism_model.png")
    pdf_path = os.path.join(OUT_DIR, "fig5_mechanism_model.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved Figure 5 to {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()
