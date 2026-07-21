"""
sl_battery.py — Shared 5-layer synthetic-lethality validation battery.

Replicates the canonical MBD4 confound stress-test doctrine
(canonical_atr_wee1_rerun.py + fig3_stress_tests.py) as a reusable harness,
operating on the DepMap 24Q4 CRISPR Chronos matrix (rows=ModelID, cols=genes).

The 5 layers (a candidate must be reported with ALL five populated, or an
explicit reason + capability-build note):
  1. CRISPR primary        — one-sided Mann-Whitney (mutant more dependent), Cohen's d, BH-FDR
  2. MSI purge             — re-test excluding MSI/hypermutator lines (ModelSubtypeFeatures~"MSI")
  3. TP53-hijack strat     — re-test within TP53-mutant lines only
  4. Leave-one-out (LOO)   — drop each mutant line; report WORST-case p (max_p)
  5. Lineage trap          — per-lineage split (OncotreePrimaryDisease); flag single-lineage signals
  + positional/CN flag     — cytoband cluster co-dependency (handled in annotation, not here)
  + broadly-essential flag — frac lines Chronos < -0.5 >= 0.50

Doctrine (from modality_fuser): a failing stress test RECLASSIFIES a candidate
(records the failing test); it does NOT delete it.

All tests: alternative="less" (mutant Chronos more negative = more dependent).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats

LOF_EFFECTS = {"frameshift", "nonsense", "ess_splice", "start_lost", "stop_lost"}
MIN_N = 5


def pooled_cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0.0
    sp = np.sqrt(((na-1)*a.std(ddof=1)**2 + (nb-1)*b.std(ddof=1)**2) / (na+nb-2))
    return float((a.mean()-b.mean())/sp) if sp > 0 else 0.0


def mw_less(lof: np.ndarray, wt: np.ndarray):
    """One-sided Mann-Whitney (lof < wt). Returns (delta_mean, p)."""
    if len(lof) < MIN_N or len(wt) < MIN_N:
        return float("nan"), 1.0
    delta = float(np.mean(lof) - np.mean(wt))
    try:
        _, p = stats.mannwhitneyu(lof, wt, alternative="less")
    except ValueError:
        p = 1.0
    return delta, float(p)


class Battery:
    """Holds the CRISPR matrix, mutation frame, and model metadata; runs the battery."""

    def __init__(self, crispr: pd.DataFrame, muts: pd.DataFrame, model: pd.DataFrame):
        # crispr: index=ModelID, columns=gene (ENTREZ suffix already stripped)
        self.crispr = crispr
        self.universe = set(crispr.index)
        self.muts = muts  # columns incl gene_symbol, effect, coding, ModelID (ACH renamed)
        self.model = model  # index=ModelID
        # Precompute MSI + TP53 + lineage maps on ModelID
        msf = model["ModelSubtypeFeatures"].astype(str) if "ModelSubtypeFeatures" in model.columns else pd.Series("", index=model.index)
        self.msi_ids = set(model.index[msf.str.contains("MSI", case=False, na=False)])
        tp = muts[muts["gene_symbol"] == "TP53"]
        self.tp53_ids = set(tp["ModelID"]) & self.universe
        self.lineage = model["OncotreePrimaryDisease"].astype(str) if "OncotreePrimaryDisease" in model.columns else pd.Series("NA", index=model.index)

    # ---- mutant / WT set construction ----
    def sets_for(self, driver: str, mode: str):
        g = self.muts[self.muts["gene_symbol"] == driver]
        if mode == "LOF":
            mut = set(g[g["effect"].isin(LOF_EFFECTS)]["ModelID"])
        else:  # ACT: driver-annotated missense if >=5 else any missense
            drv = set(g[(g["effect"] == "missense") & (g.get("cancer_driver", False) == True)]["ModelID"])
            mut = drv if len(drv & self.universe) >= MIN_N else set(g[g["effect"] == "missense"]["ModelID"])
        mut &= self.universe
        wt = self.universe - set(g["ModelID"])  # no alteration of gene at all
        return mut, wt

    # ---- single differential test for one partner gene ----
    def _test_partner(self, partner: str, mut_ids: set, wt_ids: set):
        if partner not in self.crispr.columns:
            return None
        col = self.crispr[partner]
        lof_v = col.reindex(list(mut_ids)).dropna().astype(float).values
        wt_v = col.reindex(list(wt_ids)).dropna().astype(float).values
        delta, p = mw_less(lof_v, wt_v)
        d = pooled_cohens_d(lof_v, wt_v)
        return dict(n_mut=len(lof_v), n_wt=len(wt_v), delta=delta, p=p, cohens_d=d)

    def broadly_essential(self, partner: str) -> float:
        if partner not in self.crispr.columns:
            return float("nan")
        col = self.crispr[partner].dropna().astype(float)
        return float((col < -0.5).mean())

    # ---- the full 5-layer battery for one (driver, partner) axis ----
    def run_axis(self, driver: str, partner: str, mode: str) -> dict:
        mut, wt = self.sets_for(driver, mode)
        out = dict(driver=driver, partner=partner, mode=mode,
                   n_mut_total=len(mut), n_wt_total=len(wt))

        # Layer 1: CRISPR primary
        prim = self._test_partner(partner, mut, wt)
        if prim is None:
            out["status"] = "partner_not_in_crispr"
            return out
        out.update({f"L1_{k}": v for k, v in prim.items()})

        # Layer 2: MSI purge (exclude MSI lines from both arms)
        mut_np, wt_np = mut - self.msi_ids, wt - self.msi_ids
        r = self._test_partner(partner, mut_np, wt_np)
        out.update({"L2_msipurge_p": r["p"] if r else float("nan"),
                    "L2_msipurge_delta": r["delta"] if r else float("nan"),
                    "L2_msipurge_n_mut": r["n_mut"] if r else 0})

        # Layer 3: TP53-stratified (within TP53-mutant lines)
        mut_t, wt_t = mut & self.tp53_ids, wt & self.tp53_ids
        r = self._test_partner(partner, mut_t, wt_t)
        out.update({"L3_tp53_p": r["p"] if r else float("nan"),
                    "L3_tp53_delta": r["delta"] if r else float("nan"),
                    "L3_tp53_n_mut": r["n_mut"] if r else 0})

        # Layer 4: leave-one-out (drop each mutant, worst-case p)
        loo_ps = []
        mut_list = sorted(mut)
        if len(mut_list) >= MIN_N + 1:
            for drop in mut_list:
                r = self._test_partner(partner, set(mut_list) - {drop}, wt)
                if r:
                    loo_ps.append(r["p"])
        out["L4_loo_max_p"] = float(max(loo_ps)) if loo_ps else float("nan")
        out["L4_loo_n_iter"] = len(loo_ps)

        # Layer 5: lineage trap (largest-lineage vs rest; flag single-lineage)
        lin_of_mut = self.lineage.reindex(list(mut)).dropna()
        lineage_flag = None
        if len(lin_of_mut) >= MIN_N:
            top_lin = lin_of_mut.value_counts().index[0]
            n_top = int((lin_of_mut == top_lin).sum())
            frac_top = n_top / len(lin_of_mut)
            # test within top lineage only
            mut_L = {m for m in mut if self.lineage.get(m, "NA") == top_lin}
            wt_L = {w for w in wt if self.lineage.get(w, "NA") == top_lin}
            r = self._test_partner(partner, mut_L, wt_L)
            out.update({"L5_top_lineage": top_lin, "L5_frac_in_top": round(frac_top, 3),
                        "L5_within_lineage_p": r["p"] if r else float("nan"),
                        "L5_within_lineage_n_mut": r["n_mut"] if r else 0})
            lineage_flag = frac_top >= 0.80  # >=80% of mutants from one lineage = lineage-confined
        out["flag_lineage_confined"] = bool(lineage_flag) if lineage_flag is not None else None

        # Extra flags
        out["flag_broadly_essential"] = self.broadly_essential(partner) >= 0.50
        out["broadly_essential_frac"] = round(self.broadly_essential(partner), 3)

        # ---- Verdict / tier (reclassify, never delete) ----
        out.update(self._verdict(out))
        return out

    @staticmethod
    def _verdict(o: dict) -> dict:
        """Assign survives-battery flags + a reclassification reason if it fails."""
        p1 = o.get("L1_p", 1.0)
        reasons = []
        survives = True
        # primary significance (nominal gate; FDR applied at panel level)
        if not (p1 < 0.05):
            survives = False; reasons.append("primary_not_sig(p>=0.05)")
        # MSI purge must survive (signal not MMR-driven)
        if not (o.get("L2_msipurge_p", 1.0) < 0.10) and o.get("L2_msipurge_n_mut", 0) >= MIN_N:
            reasons.append("fails_msi_purge")
            survives = False
        # TP53 strat: informative only if enough TP53-mut mutants
        if o.get("L3_tp53_n_mut", 0) >= MIN_N and not (o.get("L3_tp53_p", 1.0) < 0.10):
            reasons.append("attenuates_in_tp53_strat")
        # LOO robustness
        if not np.isnan(o.get("L4_loo_max_p", np.nan)) and not (o["L4_loo_max_p"] < 0.10):
            reasons.append("fails_leave_one_out(max_p>=0.10)")
            survives = False
        # lineage confinement
        if o.get("flag_lineage_confined") is True:
            reasons.append("lineage_confined")
        if o.get("flag_broadly_essential") is True:
            reasons.append("broadly_essential")
            survives = False
        return dict(survives_battery=bool(survives),
                    reclassification="; ".join(reasons) if reasons else "none")


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR for a list/array of p-values."""
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if ok.sum() > 0:
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out
