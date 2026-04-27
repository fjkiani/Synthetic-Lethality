#!/usr/bin/env python3
"""
Figure 1: Multimodal Evidence Matrix Heatmap
MBD4-LOF therapeutic vulnerability assessment across 6 axes × 7 modalities.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

import json
import os

# Try local dir, then check upward if needed
script_dir = os.path.dirname(os.path.abspath(__file__))
json_path = os.path.join(script_dir, 'final_evidence_matrix_output.json')
if not os.path.exists(json_path):
    print(f"Warning: {json_path} not found locally, using default structure.")

with open(json_path, 'r') as f:
    data = json.load(f)

axes = data['axes']
modalities = data['modalities']
matrix = data['matrix']
tiers = data['tiers']

# Color map matching string states directly
colors = {
    'MISSING': '#E8E8E8',    # gray
    'POSITIVE': '#E74C3C',   # red
    'NEGATIVE': '#3498DB',   # blue
    'MIXED': '#F39C12',      # amber
    'CONFOUNDED': '#9B59B6', # purple
}

fig, ax = plt.subplots(figsize=(10, 6))

for i in range(len(axes)):
    for j in range(len(modalities)):
        val = matrix[i][j]
        rect = plt.Rectangle((j, len(axes) - 1 - i), 1, 1,
                              facecolor=colors[val], edgecolor='white', linewidth=2)
        ax.add_patch(rect)
        ax.text(j + 0.5, len(axes) - 1 - i + 0.5, val,
                ha='center', va='center', fontsize=7, fontweight='bold',
                color='white' if val in ['POSITIVE', 'NEGATIVE', 'CONFOUNDED'] else '#333')

ax.set_xlim(0, len(modalities))
ax.set_ylim(0, len(axes))
ax.set_xticks([j + 0.5 for j in range(len(modalities))])
ax.set_xticklabels(modalities, fontsize=8, ha='center')
ax.set_yticks([i + 0.5 for i in range(len(axes))])
ax.set_yticklabels(reversed(axes), fontsize=9)
ax.tick_params(axis='both', length=0)

# Tier annotations on right side
for i, tier in enumerate(tiers):
    color = {'Validated SL': '#27AE60', 'Strong': '#E67E22', 'Mechanistic': '#F1C40F', 'Sentinel': '#1ABC9C',
             'Negative': '#E74C3C', 'Confounded': '#9B59B6', 'Insufficient': '#BDC3C7'}.get(tier, '#BDC3C7')
    ax.text(len(modalities) + 0.2, len(axes) - 1 - i + 0.5, tier,
            ha='left', va='center', fontsize=8, fontweight='bold', color=color)

ax.set_title('MBD4-LOF Multimodal Evidence Matrix', fontsize=14, fontweight='bold', pad=15)
ax.text(len(modalities) + 0.2, len(axes) + 0.3, 'Tier', fontsize=9, fontweight='bold', color='#666')

# Legend
legend_patches = [mpatches.Patch(facecolor=colors[k], edgecolor='gray', label=k) for k in ['POSITIVE', 'NEGATIVE', 'MIXED', 'MISSING', 'CONFOUNDED']]
ax.legend(handles=legend_patches, loc='lower center', bbox_to_anchor=(0.5, -0.18),
          ncol=5, fontsize=8, frameon=False)

plt.tight_layout()
plt.savefig('fig1_evidence_matrix.png', dpi=300, bbox_inches='tight')
plt.savefig('fig1_evidence_matrix.pdf', bbox_inches='tight')
print("Figure 1 saved.")
