#!/usr/bin/env python3
"""
Figure 3: Four Confound Stress Tests
4-panel figure confirming ceralasertib signal is MBD4-specific.

Panel C p-values are loaded from the frozen LOO receipt (matches canonical max p, not rounded literals).
"""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

_BUNDLE = Path(__file__).resolve().parent.parent.parent
_LOO_RECEIPT = (
    _BUNDLE
    / "artifacts/canonical_atr_wee1_rerun_20260405"
    / "ceralasertib_ln_ic50_leave_one_out_MANUSCRIPT_RECEIPT.json"
)
with open(_LOO_RECEIPT, encoding="utf-8") as _f:
    _loo = json.load(_f)
_loo_iters = sorted(
    _loo["iterations_sorted_by_cosmic_id_then_index"],
    key=lambda x: x["iteration_index_1_based"],
)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# --- Panel A: MSI purge (DepMap ModelSubtypeFeatures) ---
ax = axes[0, 0]
groups = ['Full cohort\n(n=14 LOF)', 'MSS only\n(n=10 LOF)']
deltas = [-0.74, -0.91]
pvals = [0.021, 0.015]
colors_a = ['#E74C3C', '#C0392B']

bars = ax.bar(groups, deltas, color=colors_a, width=0.5, edgecolor='white', linewidth=1.5)
for i, (bar, p) in enumerate(zip(bars, pvals)):
    y = bar.get_height()
    symbol = '**' if p < 0.01 else '*'
    ax.text(bar.get_x() + bar.get_width()/2, y - 0.05,
            f'p={p}\n{symbol}', ha='center', va='top', fontsize=9, fontweight='bold', color='white')

ax.set_ylabel('Δ LN_IC50 (LOF − WT)', fontsize=10, fontweight='bold')
ax.set_title('A. MSI Purge (DepMap ModelSubtypeFeatures)', fontsize=12, fontweight='bold')
ax.axhline(0, color='gray', linewidth=0.8, linestyle='--')
ax.set_ylim(-1.2, 0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(0.95, 0.05, 'Signal\nSTRONGER', transform=ax.transAxes, ha='right', va='bottom',
        fontsize=10, fontweight='bold', color='#27AE60',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1', edgecolor='#27AE60'))

# --- Panel B: TP53 Hijack Check ---
ax = axes[0, 1]
groups = ['MBD4-WT\nTP53-mut\n(n=619)', 'MBD4-LOF\nTP53-mut\n(n=11)']
means = [2.127, 1.058]
colors_b = ['#D5DBDB', '#E74C3C']

bars = ax.bar(groups, means, color=colors_b, width=0.5, edgecolor='white', linewidth=1.5)
ax.annotate('', xy=(1, means[1]), xytext=(0, means[0]),
            arrowprops=dict(arrowstyle='->', color='#C0392B', lw=2.5))
ax.text(0.5, (means[0] + means[1])/2 + 0.1,
        'Δ = −1.07\np = 0.003 **\nd = −0.74',
        ha='center', fontsize=9, fontweight='bold', color='#C0392B',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDEDEC'))

ax.set_ylabel('GDSC2 LN_IC50', fontsize=10, fontweight='bold')
ax.set_title('B. TP53 Hijack Check', fontsize=12, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(0.95, 0.05, '>1 log-unit\nbeyond TP53', transform=ax.transAxes, ha='right', va='bottom',
        fontsize=10, fontweight='bold', color='#27AE60',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1', edgecolor='#27AE60'))

# --- Panel C: Leave-One-Out ---
ax = axes[1, 0]
loo_pvals = [
    float(x["p_value_one_sided_mwu_lof_less_wt"]) for x in _loo_iters
]
loo_labels = [f'#{x["iteration_index_1_based"]}' for x in _loo_iters]

colors_c = ['#E74C3C' if p < 0.05 else '#F39C12' for p in loo_pvals]
bars = ax.barh(loo_labels, loo_pvals, color=colors_c, height=0.6, edgecolor='white')

ax.axvline(0.05, color='#E74C3C', linewidth=2, linestyle='--', label='α = 0.05')
ax.axvline(0.10, color='#F39C12', linewidth=1.5, linestyle=':', label='α = 0.10')

ax.set_xlabel('p-value (one-sided Mann-Whitney)', fontsize=10, fontweight='bold')
ax.set_ylabel('LOO iteration', fontsize=10, fontweight='bold')
ax.set_title('C. Leave-One-Out Robustness', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, 0.12)
ax.text(0.05, 0.95, '14/14\nrobust', transform=ax.transAxes, ha='left', va='top',
        fontsize=11, fontweight='bold', color='#27AE60',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1', edgecolor='#27AE60'))

# --- Panel D: Lineage Trap ---
ax = axes[1, 1]
lineages = ['Bowel\n(n=5 LOF)', 'Non-Bowel\n(n=9 LOF)', 'All\n(n=14 LOF)']
deltas_d = [-0.72, -0.87, -0.74]
pvals_d = [0.13, 0.025, 0.021]
colors_d = ['#F39C12', '#E67E22', '#E74C3C']

bars = ax.bar(lineages, deltas_d, color=colors_d, width=0.5, edgecolor='white', linewidth=1.5)
for bar, p in zip(bars, pvals_d):
    y = bar.get_height()
    symbol = '*' if p < 0.05 else 'ns' if p > 0.10 else '†'
    ax.text(bar.get_x() + bar.get_width()/2, y - 0.04,
            f'p={p}\n{symbol}', ha='center', va='top', fontsize=9, fontweight='bold', color='white')

ax.set_ylabel('Δ LN_IC50 (LOF − WT)', fontsize=10, fontweight='bold')
ax.set_title('D. Lineage Trap', fontsize=12, fontweight='bold')
ax.axhline(0, color='gray', linewidth=0.8, linestyle='--')
ax.set_ylim(-1.2, 0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(0.95, 0.05, '8 lineages\nrepresented', transform=ax.transAxes, ha='right', va='bottom',
        fontsize=10, fontweight='bold', color='#27AE60',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1', edgecolor='#27AE60'))

plt.suptitle('Confound Stress Tests — Ceralasertib Signal Validation',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('fig3_stress_tests.png', dpi=300, bbox_inches='tight')
plt.savefig('fig3_stress_tests.pdf', bbox_inches='tight')
print("Figure 3 saved.")
