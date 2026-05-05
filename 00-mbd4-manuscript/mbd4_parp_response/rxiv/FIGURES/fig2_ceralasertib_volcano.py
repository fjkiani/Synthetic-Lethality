#!/usr/bin/env python3
"""
Figure 2: Ceralasertib Sensitivity — MBD4-LOF vs WT
Box/strip plot of GDSC2 LN_IC50 values.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# --- Data from GDSC2 stratification results ---
# Load actual extracted GDSC2 Ceralasertib dataset distributions
import os
wt_cera_path = os.path.join(os.path.dirname(__file__), 'wt_cera.npy')
lof_cera_path = os.path.join(os.path.dirname(__file__), 'lof_cera.npy')

if os.path.exists(wt_cera_path) and os.path.exists(lof_cera_path):
    wt_data = np.load(wt_cera_path)
    lof_data = np.load(lof_cera_path)
else:
    raise FileNotFoundError("CRITICAL ERROR: Data cache is missing. Extract the true data from GDSC2/Omics first.")

fig, ax = plt.subplots(figsize=(6, 7))

# Box plots
bp = ax.boxplot([wt_data, lof_data], positions=[1, 2], widths=0.5,
                patch_artist=True, showfliers=False,
                medianprops=dict(color='black', linewidth=2),
                whiskerprops=dict(linewidth=1.5),
                capprops=dict(linewidth=1.5))

bp['boxes'][0].set_facecolor('#D5DBDB')
bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_facecolor('#E74C3C')
bp['boxes'][1].set_alpha(0.5)

# Strip plot for LOF (individual points)
jitter = np.linspace(-0.08, 0.08, len(lof_data))
ax.scatter(np.full(len(lof_data), 2) + jitter, lof_data,
           color='#C0392B', s=60, alpha=0.8, edgecolors='white', linewidth=0.8, zorder=5)

# WT subsample scatter (deterministic top 50 to avoid random)
wt_sample = np.sort(wt_data)[:50]
jitter_wt = np.linspace(-0.15, 0.15, len(wt_sample))
ax.scatter(np.full(len(wt_sample), 1) + jitter_wt, wt_sample,
           color='#7F8C8D', s=15, alpha=0.3, zorder=4)

# Significance bracket
y_max = max(wt_data.max(), lof_data.max()) + 0.5
ax.plot([1, 1, 2, 2], [y_max, y_max + 0.15, y_max + 0.15, y_max],
        color='black', linewidth=1.5)
ax.text(1.5, y_max + 0.25, 'p = 0.021 *', ha='center', fontsize=11, fontweight='bold')

# Delta annotation
ax.annotate('', xy=(2.6, np.median(lof_data)), xytext=(2.6, np.median(wt_data)),
            arrowprops=dict(arrowstyle='<->', color='#E74C3C', lw=2))
ax.text(2.75, (np.median(lof_data) + np.median(wt_data)) / 2,
        'Δ = −0.74\nd = −0.50',
        ha='left', va='center', fontsize=10, color='#C0392B', fontweight='bold')

ax.set_xticks([1, 2])
ax.set_xticklabels(['WT\n(n=914)', 'MBD4 LOF\n(n=14)'], fontsize=12, fontweight='bold')
ax.set_ylabel('Ceralasertib LN_IC50 (GDSC2)', fontsize=13, fontweight='bold')
ax.set_title('MBD4-LOF Sensitivity to ATR Inhibition', fontsize=14, fontweight='bold', pad=15)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=11)

plt.tight_layout()
plt.savefig('fig2_ceralasertib_volcano.png', dpi=300, bbox_inches='tight')
plt.savefig('fig2_ceralasertib_volcano.pdf', bbox_inches='tight')
print("Figure 2 saved.")
