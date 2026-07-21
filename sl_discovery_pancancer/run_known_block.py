"""
run_known_block.py <BLOCK>  — validate the known-SL panel pairs for one block
through the full 5-layer battery. Writes per-block CSV + JSON checkpoints so a
timeout never loses completed work.

Usage: python run_known_block.py A     (blocks: A, B, C)
"""
import sys, json, time
from pathlib import Path
import pandas as pd

sys.path.insert(0, "/mnt/shared-workspace/sl_discovery")
sys.path.insert(0, "/mnt/shared-workspace/repo/Synthetic-Lethality")

from sl_battery import Battery, bh_fdr
from known_sl_panel import KNOWN_SL_PANEL

BLOCK = sys.argv[1] if len(sys.argv) > 1 else "A"
CACHE = "/mnt/shared-workspace/depmap_cache"
OUT = Path("/mnt/shared-workspace/sl_discovery/results")
OUT.mkdir(parents=True, exist_ok=True)
ckpt = OUT / f"known_block_{BLOCK}.csv"

t0 = time.time()
print(f"[{BLOCK}] loading data...", flush=True)
crispr = pd.read_parquet(f"{CACHE}/crispr_gene_effect.parquet")
crispr.columns = [c.split(" (")[0] if " (" in c else c for c in crispr.columns]
model = pd.read_parquet(f"{CACHE}/depmap_model.parquet")
muts = pd.read_parquet("/mnt/shared-workspace/sl_discovery/mutations_in_crispr.parquet")
muts = muts.rename(columns={"ACH": "ModelID"})
bat = Battery(crispr=crispr, muts=muts, model=model)
print(f"[{BLOCK}] battery ready ({time.time()-t0:.1f}s). MSI={len(bat.msi_ids)} TP53={len(bat.tp53_ids)}", flush=True)

pairs = [e for e in KNOWN_SL_PANEL if e["block"] == BLOCK]
print(f"[{BLOCK}] {len(pairs)} pairs to validate", flush=True)

rows = []
for i, e in enumerate(pairs, 1):
    ts = time.time()
    r = bat.run_axis(e["driver"], e["partner"], e["mode"])
    r["drug_axis"] = e["drug_axis"]
    r["cn_event"] = e["cn_event"]
    r["citation"] = e["citation"]
    r["known"] = True
    rows.append(r)
    # per-pair checkpoint
    pd.DataFrame(rows).to_csv(ckpt, index=False)
    print(f"[{BLOCK}] {i}/{len(pairs)} {e['driver']}->{e['partner']} ({e['mode']}) "
          f"n_mut={r.get('n_mut_total')} L1_p={r.get('L1_p')} "
          f"survives={r.get('survives_battery')} reclass={r.get('reclassification')} "
          f"({time.time()-ts:.1f}s)", flush=True)

df = pd.DataFrame(rows)
# panel-level BH-FDR on primary p across this block
if "L1_p" in df.columns:
    df["L1_fdr"] = bh_fdr(df["L1_p"].values)
df.to_csv(ckpt, index=False)
df.to_json(OUT / f"known_block_{BLOCK}.json", orient="records", indent=1)
print(f"[{BLOCK}] DONE {len(df)} pairs in {time.time()-t0:.1f}s -> {ckpt}", flush=True)
