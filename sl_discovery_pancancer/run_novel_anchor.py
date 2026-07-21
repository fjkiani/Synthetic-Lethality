"""
run_novel_anchor.py — NOVEL SL discovery via anchor-and-extend.

For each literature anchor driver, run genome-wide CRISPR differential-dependency
(SLEngine) to discover NEW partners, drop any partner already in the known panel,
then run the full 5-layer battery (sl_battery) on the top novel partners.

Writes checkpointed CSVs so a timeout never loses completed anchors.
Usage: python run_novel_anchor.py        (runs all anchors)
       python run_novel_anchor.py PTEN ARID1A   (subset)
"""
import sys, json, time
from pathlib import Path
import pandas as pd, numpy as np

sys.path.insert(0, "/mnt/shared-workspace/sl_discovery")
sys.path.insert(0, "/mnt/shared-workspace/repo/Synthetic-Lethality")

from sl_agent.core.sl_engine import SLEngine
from sl_agent.core.models import SLQueryInput, MutationType
from sl_battery import Battery, bh_fdr, LOF_EFFECTS
from known_sl_panel import NOVEL_ANCHORS, KNOWN_SL_PANEL

CACHE = "/mnt/shared-workspace/depmap_cache"
OUT = Path("/mnt/shared-workspace/sl_discovery/results"); OUT.mkdir(parents=True, exist_ok=True)

# partners already known (exclude from "novel")
KNOWN_PARTNERS = {(e["driver"], e["partner"]) for e in KNOWN_SL_PANEL}
# common pan-essential paralog partners we consider "known mechanism" per anchor
KNOWN_BY_DRIVER = {}
for e in KNOWN_SL_PANEL:
    KNOWN_BY_DRIVER.setdefault(e["driver"], set()).add(e["partner"])


def sets_for(muts, universe, driver, mode):
    g = muts[muts["gene_symbol"] == driver]
    if mode == "LOF":
        mut = set(g[g["effect"].isin(LOF_EFFECTS)]["ModelID"])
    else:
        drv = set(g[(g["effect"] == "missense") & (g.get("cancer_driver", False) == True)]["ModelID"])
        mut = drv if len(drv & universe) >= 5 else set(g[g["effect"] == "missense"]["ModelID"])
    mut &= universe
    wt = universe - set(g["ModelID"])
    return mut, wt


def main():
    t0 = time.time()
    want = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    print("loading data...", flush=True)
    crispr = pd.read_parquet(f"{CACHE}/crispr_gene_effect.parquet")
    crispr.columns = [c.split(" (")[0] if " (" in c else c for c in crispr.columns]
    model = pd.read_parquet(f"{CACHE}/depmap_model.parquet")
    muts = pd.read_parquet("/mnt/shared-workspace/sl_discovery/mutations_in_crispr.parquet").rename(columns={"ACH": "ModelID"})
    eng = SLEngine(crispr_df=crispr, sample_info_df=model, mutation_df=muts, depmap_release="24Q4")
    bat = Battery(crispr=crispr, muts=muts, model=model)
    universe = set(crispr.index)
    print(f"engine+battery ready ({time.time()-t0:.1f}s)", flush=True)

    anchors = [a for a in NOVEL_ANCHORS if (want is None or a["driver"] in want)]
    disc_rows, battery_rows = [], []
    for a in anchors:
        drv, mode = a["driver"], a["mode"]
        mut, wt = sets_for(muts, universe, drv, mode)
        if len(mut) < 5:
            print(f"[skip] {drv} n_mut={len(mut)} <5", flush=True); continue
        q = SLQueryInput(gene=drv, mutation_type=MutationType.ANY, cancer_type="pan_cancer",
                         top_n_partners=60, fdr_cutoff=0.25, delta_dep_cutoff=0.15, include_codependency=True)
        ts = time.time()
        _, partners, _ = eng.compute_sl_partners(q, mutant_ids=list(mut), wt_ids=list(wt))
        # rank, drop known partners for THIS driver, drop broadly-essential noise
        novel = []
        for p in partners:
            if p.gene in KNOWN_BY_DRIVER.get(drv, set()):
                continue
            be = bat.broadly_essential(p.gene)
            disc_rows.append(dict(driver=drv, mode=mode, partner=p.gene,
                                  delta=round(p.delta_dependency, 4), fdr=round(p.fdr, 5),
                                  cohens_d=round(p.effect_size_cohend, 3),
                                  codep_r=round(p.codependency_r, 3) if p.codependency_r is not None else None,
                                  pathway=p.pathway, broadly_ess=round(be, 3), known_partner=False))
            if p.fdr <= 0.10 and be < 0.90:   # promote credible novel leads to battery
                novel.append(p.gene)
        print(f"[{drv}] {len(partners)} partners, {len(novel)} novel leads -> battery ({time.time()-ts:.1f}s)", flush=True)
        # full battery on top novel leads (cap 8 per anchor for compute)
        for gene in novel[:8]:
            r = bat.run_axis(drv, gene, mode); r["known"] = False
            battery_rows.append(r)
        pd.DataFrame(disc_rows).to_csv(OUT/"novel_discovery.csv", index=False)
        if battery_rows:
            pd.DataFrame(battery_rows).to_csv(OUT/"novel_battery.csv", index=False)

    dd = pd.DataFrame(disc_rows)
    dd.to_csv(OUT/"novel_discovery.csv", index=False)
    bb = pd.DataFrame(battery_rows)
    if len(bb):
        bb["L1_fdr"] = bh_fdr(bb["L1_p"].values)
        bb.to_csv(OUT/"novel_battery.csv", index=False)
        bb.to_json(OUT/"novel_battery.json", orient="records", indent=1)
    print(f"DONE anchors={len(anchors)} discovery_rows={len(dd)} battery_rows={len(bb)} in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
