#!/usr/bin/env python3
"""
Figure 5: Dual-Axis Therapeutic Vulnerability Model
Pathway diagram: MBD4 LOF → BER stress → fork vulnerability → Cytidine + ATRi axes.
PARP1 upregulation shown as adaptive biomarker (side panel), NOT as a therapeutic axis.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(12, 9))
ax.set_xlim(0, 12)
ax.set_ylim(-1.0, 9)
ax.axis('off')

def add_box(ax, x, y, w, h, text, color, fontsize=9, fontweight='bold', textcolor='white'):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                          facecolor=color, edgecolor='white', linewidth=2, zorder=5)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight, color=textcolor, zorder=6)

def add_arrow(ax, start, end, color='#2C3E50', style='->', lw=2.5):
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                               connectionstyle='arc3,rad=0'))

# === TOP: MBD4 LOF ===
add_box(ax, 4.5, 7.8, 3, 0.7, 'MBD4 Loss-of-Function', '#8E44AD', fontsize=12)

# === DOWN: BER Deficiency ===
add_arrow(ax, (6, 7.8), (6, 7.3))
add_box(ax, 4.2, 6.5, 3.6, 0.7, 'BER Pathway Deficiency\n(Glycosylase loss)', '#9B59B6', fontsize=10)

# === DOWN: Unresolved Base Damage ===
add_arrow(ax, (6, 6.5), (6, 6.0))
add_box(ax, 3.5, 5.2, 5.0, 0.7, 'Unprocessed U:G Mismatches\nat CpG Sites', '#7D3C98', fontsize=10)

# === FORK: Two therapeutic axes diverge ===

# --- LEFT BRANCH: AXIS 1 — Cytidine Analog SL ---
add_arrow(ax, (4.5, 5.2), (2.5, 4.7), color='#27AE60', lw=2.5)
add_box(ax, 0.5, 3.8, 4.0, 0.7, 'Accumulated U:G\nMismatch Substrate', '#27AE60', fontsize=10)

add_arrow(ax, (2.5, 3.8), (2.5, 3.2), color='#27AE60')
add_box(ax, 0.3, 2.3, 4.4, 0.8, 'Gemcitabine / Cytarabine\nExploit BER Substrate', '#1E8449', fontsize=10)

add_arrow(ax, (2.5, 2.3), (2.5, 1.7), color='#27AE60')
add_box(ax, 0.1, 0.8, 4.8, 0.8, 'AXIS 1: CYTIDINE ANALOGS\n(Validated — isogenic + rescue)', '#145A32', fontsize=11, textcolor='#D5F5E3')

# --- RIGHT BRANCH: AXIS 2 — ATR Inhibition ---
add_arrow(ax, (7.5, 5.2), (9.5, 4.7), color='#2980B9', lw=2.5)
add_box(ax, 7.5, 3.8, 4.0, 0.7, 'Fork-Blocking Lesions\n→ Replication Stress', '#2980B9', fontsize=10)

add_arrow(ax, (9.5, 3.8), (9.5, 3.2), color='#2980B9')
add_box(ax, 7.3, 2.3, 4.4, 0.8, 'ATR Checkpoint\nDependency', '#2471A3', fontsize=10)

add_arrow(ax, (9.5, 2.3), (9.5, 1.7), color='#2980B9')
add_box(ax, 7.1, 0.8, 4.8, 0.8, 'AXIS 2: ATR INHIBITION\n(p=0.003, d=−0.74; 4/4 stress)', '#1A5276', fontsize=11, textcolor='#D6EAF8')

# === CENTER: Convergence arrow ===
ax.annotate('', xy=(6, 0.2), xytext=(2.5, 0.8),
            arrowprops=dict(arrowstyle='->', color='#27AE60', lw=2, linestyle='--'))
ax.annotate('', xy=(6, 0.2), xytext=(9.5, 0.8),
            arrowprops=dict(arrowstyle='->', color='#2980B9', lw=2, linestyle='--'))

# Convergence point
convergence = FancyBboxPatch((4.2, -0.5), 3.6, 0.6, boxstyle="round,pad=0.1",
                              facecolor='#1C2833', edgecolor='#F39C12', linewidth=3, zorder=7)
ax.add_patch(convergence)
ax.text(6, -0.2, '☠ REPLICATION FORK COLLAPSE', ha='center', va='center',
        fontsize=10, fontweight='bold', color='#F39C12', zorder=8)

# === SIDE PANEL: PARP1 Biomarker (NOT an axis) ===
# Small dashed box to the right, labeled as biomarker
parp_box = FancyBboxPatch((0.2, -0.5), 3.5, 0.5, boxstyle="round,pad=0.1",
                           facecolor='#F5F5F5', edgecolor='#95A5A6', linewidth=1.5,
                           linestyle='--', zorder=5)
ax.add_patch(parp_box)
ax.text(1.95, -0.25, 'No PARP1 Upregulation\n(No Selective Transcriptional Induction)', ha='center', va='center',
        fontsize=8, fontweight='normal', color='#7F8C8D', fontstyle='italic', zorder=6)

# Dashed arrow from BER stress to PARP1 biomarker box
ax.annotate('', xy=(1.95, -0.0), xytext=(3.5, 5.2),
            arrowprops=dict(arrowstyle='->', color='#95A5A6', lw=1.5, linestyle=':', 
                           connectionstyle='arc3,rad=-0.3'))

# === Branch labels ===
ax.text(2.5, 4.9, 'BER substrate\naccumulation', ha='center', fontsize=8,
        color='#27AE60', fontstyle='italic')
ax.text(9.5, 4.9, 'Constitutive\nfork stalling', ha='center', fontsize=8,
        color='#2980B9', fontstyle='italic')

# === Title ===
ax.text(6, 9.0, 'MBD4-LOF Dual Therapeutic Framework: Cytidine + ATRi',
        ha='center', fontsize=15, fontweight='bold', color='#2C3E50')

# === Evidence badges ===
ax.text(7.8, -0.7, 'RUO — Research Use Only', fontsize=8, color='#7F8C8D', fontstyle='italic')

plt.tight_layout()
plt.savefig('fig5_mechanism_model.png', dpi=300, bbox_inches='tight')
plt.savefig('fig5_mechanism_model.pdf', bbox_inches='tight')
print("Figure 5 saved — Cytidine/ATRi dual strategy, PARP1 demoted to biomarker.")
