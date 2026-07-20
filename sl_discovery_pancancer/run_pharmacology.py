"""
run_pharmacology.py — GDSC2 biomarker-stratified pharmacology confirmation for
targetable known-SL axes, keyed on SANGER_MODEL_ID (the fix from the prior phase).

For each (driver, drug-axis) it computes:
  - LOF vs WT one-sided Wilcoxon on LN_IC50 (LOF more sensitive = lower), Cohen's d, BH-FDR
  - MSI-purge re-test (exclude MSI/hypermutator lines)
  - falsification: mechanistically-sparing agents must NOT show the effect

Pan-cancer (not lineage-restricted). Writes checkpointed CSV.
"""
import sys, json, time
from pathlib import Path
import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, "/mnt/shared-workspace/sl_discovery")
CACHE = Path("/mnt/shared-workspace/depmap_cache")
OUT = Path("/mnt/shared-workspace/sl_discovery/results"); OUT.mkdir(parents=True, exist_ok=True)
LOF = {"frameshift", "nonsense", "ess_splice", "start_lost", "stop_lost"}

# driver -> {axis drugs (regex), mode, falsification drugs}
AXES = {
    "PTEN":  dict(mode="LOF", drugs=["AZD8186", "Dactolisib", "AZD6482", "Buparlisib", "Pictilisib", "GSK2636771"],
                  falsify=["Taselisib", "Alpelisib"], label="PI3K-beta"),
    "BRCA1": dict(mode="LOF", drugs=["Olaparib", "Niraparib", "Talazoparib", "Rucaparib", "Veliparib"],
                  falsify=[], label="PARP"),
    "BRCA2": dict(mode="LOF", drugs=["Olaparib", "Niraparib", "Talazoparib", "Rucaparib", "Veliparib"],
                  falsify=[], label="PARP"),
    "ATM":   dict(mode="LOF", drugs=["Ceralasertib", "AZD6738", "Berzosertib", "VE-822", "Elimusertib"],
                  falsify=[], label="ATR"),
    "ARID1A":dict(mode="LOF", drugs=["Tazemetostat", "GSK126", "EPZ-6438", "PF-06821497"],
                  falsify=[], label="EZH2"),
    "STK11": dict(mode="LOF", drugs=["AZD8055", "Rapamycin", "Temsirolimus", "Dactolisib", "OSI-027"],
                  falsify=[], label="mTOR"),
    "KEAP1": dict(mode="LOF", drugs=["CB-839", "Telaglenastat", "BPTES"],
                  falsify=[], label="glutaminase"),
    "VHL":   dict(mode="LOF", drugs=["PT2385", "Belzutifan", "PT2977"],
                  falsify=[], label="HIF2A"),
    "RB1":   dict(mode="LOF", drugs=["Palbociclib", "Ribociclib", "Abemaciclib"],
                  falsify=[], label="CDK4/6", expect_resistant=True),  # RB-loss -> RESISTANT to CDK4/6i (direction check)
}


def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return 0.0
    sp = np.sqrt(((na-1)*a.std(ddof=1)**2 + (nb-1)*b.std(ddof=1)**2)/(na+nb-2))
    return float((a.mean()-b.mean())/sp) if sp > 0 else 0.0


def wilcox(mut, wt, alt="less"):
    if len(mut) < 5 or len(wt) < 5: return 1.0
    try: return float(stats.mannwhitneyu(mut, wt, alternative=alt)[1])
    except ValueError: return 1.0


def main():
    t0 = time.time()
    model = pd.read_parquet(CACHE/"depmap_model.parquet")
    muts = pd.read_parquet("/mnt/shared-workspace/sl_discovery/mutations_in_crispr.parquet").rename(columns={"ACH":"ModelID"})
    ach2sanger = model["SangerModelID"].dropna().to_dict()
    msf = model["ModelSubtypeFeatures"].astype(str)
    msi_sanger = {ach2sanger[a] for a in model.index[msf.str.contains("MSI", case=False, na=False)] if a in ach2sanger}
    print(f"loaded model+muts ({time.time()-t0:.1f}s); MSI Sanger n={len(msi_sanger)}", flush=True)

    g2 = pd.read_excel(CACHE/"gdsc2.xlsx")
    print(f"GDSC2 {g2.shape}; drugs={g2['DRUG_NAME'].nunique()}", flush=True)

    rows = []
    for drv, cfg in AXES.items():
        g = muts[muts["gene_symbol"] == drv]
        lof_ach = set(g[g["effect"].isin(LOF)]["ModelID"])
        anyalt_ach = set(g["ModelID"])
        all_ach = set(model.index)
        lof_s = {ach2sanger[a] for a in lof_ach if a in ach2sanger}
        wt_s = {ach2sanger[a] for a in (all_ach - anyalt_ach) if a in ach2sanger}
        drugs = cfg["drugs"] + cfg["falsify"]
        for drug in drugs:
            sub = g2[g2["DRUG_NAME"].str.contains(drug, case=False, na=False)]
            if sub.empty:
                continue
            m = sub[sub["SANGER_MODEL_ID"].isin(lof_s)]["LN_IC50"].dropna().astype(float).values
            w = sub[sub["SANGER_MODEL_ID"].isin(wt_s)]["LN_IC50"].dropna().astype(float).values
            if len(m) < 5 or len(w) < 5:
                continue
            is_fals = drug in cfg["falsify"]
            alt = "greater" if cfg.get("expect_resistant") else "less"
            p = wilcox(m, w, alt=alt)
            # MSI purge
            mp = sub[sub["SANGER_MODEL_ID"].isin(lof_s - msi_sanger)]["LN_IC50"].dropna().astype(float).values
            wp = sub[sub["SANGER_MODEL_ID"].isin(wt_s - msi_sanger)]["LN_IC50"].dropna().astype(float).values
            p_msi = wilcox(mp, wp, alt=alt) if len(mp) >= 5 and len(wp) >= 5 else float("nan")
            rows.append(dict(driver=drv, axis_label=cfg["label"], drug=drug,
                             role=("falsification" if is_fals else "axis"),
                             direction=("resistant" if cfg.get("expect_resistant") else "sensitive"),
                             n_mut=len(m), n_wt=len(w),
                             delta_ln_ic50=round(float(m.mean()-w.mean()), 4),
                             cohens_d=round(cohens_d(m, w), 3), p_one_sided=p,
                             p_msipurge=p_msi))
            print(f"  {drv:6s} {drug:14s} [{'FALS' if is_fals else 'axis'}] "
                  f"n={len(m)}/{len(w)} d={cohens_d(m,w):.2f} p={p:.4f}", flush=True)

    df = pd.DataFrame(rows)
    if len(df):
        df["fdr_bh"] = multipletests(df["p_one_sided"].clip(0,1).fillna(1), method="fdr_bh")[1]
    df.to_csv(OUT/"pharmacology_v2.csv", index=False)
    df.to_json(OUT/"pharmacology_v2.json", orient="records", indent=1)
    print(f"DONE {len(df)} drug-tests in {time.time()-t0:.1f}s -> {OUT}/pharmacology_v2.csv", flush=True)


if __name__ == "__main__":
    main()
