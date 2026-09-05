#!/usr/bin/env python3
"""Render the MBD4-LOF Triaxis therapeutic model and CrisPRO routing schema."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

OUT = Path(__file__).resolve().parent
fig, ax = plt.subplots(figsize=(15, 10))
ax.set_xlim(0, 15)
ax.set_ylim(0, 10)
ax.axis("off")

def box(x, y, w, h, text, color, size=10, edge="white"):
    patch = FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.12",facecolor=color,edgecolor=edge,linewidth=2)
    ax.add_patch(patch)
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=size,fontweight="bold",color="white",wrap=True)

def arrow(start,end,color="#34495E",style="->",lw=2.2):
    ax.annotate("",xy=end,xytext=start,arrowprops=dict(arrowstyle=style,color=color,lw=lw))

ax.text(7.5,9.65,"MBD4-LOF Triaxis Therapeutic Model",ha="center",fontsize=18,fontweight="bold",color="#17202A")
box(5.3,8.65,4.4,0.65,"MBD4 loss-of-function", "#6C3483", 13)
arrow((7.5,8.65),(7.5,8.05))
box(4.5,7.3,6.0,0.7,"Failed methyl-CpG mismatch excision", "#884EA0", 12)
arrow((6.0,7.3),(2.6,6.5),"#239B56")
arrow((7.5,7.3),(7.5,6.5),"#2471A3")
arrow((9.0,7.3),(12.4,6.5),"#B9770E")
box(0.5,5.55,4.2,0.85,"AXIS 1\nUnresolved BER substrate", "#239B56", 11)
box(5.4,5.55,4.2,0.85,"AXIS 2\nFork stalling + checkpoint dependence", "#2471A3", 11)
box(10.3,5.55,4.2,0.85,"AXIS 3\nCpG>TpG hypermutator state", "#B9770E", 11)
arrow((2.6,5.55),(2.6,4.8),"#239B56")
arrow((7.5,5.55),(7.5,4.8),"#2471A3")
arrow((12.4,5.55),(12.4,4.8),"#B9770E")
box(0.5,3.85,4.2,0.85,"Gemcitabine / Cytarabine\nsubstrate loading", "#196F3D", 10)
box(5.4,3.85,4.2,0.85,"Ceralasertib / Adavosertib\ncheckpoint removal", "#1A5276", 10)
box(10.3,3.85,4.2,0.85,"Altered peptides → candidate neoantigens\nimmune recognition", "#9C640C", 10)
arrow((2.6,3.85),(4.8,2.75),"#239B56")
arrow((7.5,3.85),(6.3,2.75),"#2471A3")
box(3.5,1.75,4.2,0.85,"SYNERGISTIC_COMBINATION_CANDIDATE\ncoordinate fork failure", "#17202A", 10, "#F4D03F")
arrow((12.4,3.85),(12.4,2.75),"#B9770E")
box(10.3,1.75,4.2,0.85,"IMMUNOTHERAPY_CHECKPOINT_BLOCKADE\nhypermutator-confirmed route", "#784212", 10, "#F5B041")
box(0.5,0.35,5.2,0.75,"PARP1 exploratory gate: 7.41 log1p TPM Q75\nindependent BRCA/HRD routes preserved", "#626567", 9, "#D5D8DC")
box(9.3,0.35,5.2,0.75,"Evidence tiers remain distinct\nAxis 1 causal | Axis 2 pharmacogenomic | Axis 3 human retrospective", "#626567", 9, "#D5D8DC")
fig.tight_layout()
for name in ["crispro_routing_flowchart.png","Figure_4.png","fig5_mechanism_model.png"]:
    fig.savefig(OUT/name,dpi=300,bbox_inches="tight")
for name in ["crispro_routing_flowchart.pdf","fig5_mechanism_model.pdf"]:
    fig.savefig(OUT/name,bbox_inches="tight")
print("Rendered MBD4-LOF Triaxis therapeutic model")
