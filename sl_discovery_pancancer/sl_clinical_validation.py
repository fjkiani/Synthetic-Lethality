"""
sl_clinical_validation.py
=========================

External-validation / trust-scoring layer for pan-cancer synthetic-lethality
candidates. This is the "next layer" on top of the 5-layer confound battery
(sl_battery.py): it does NOT claim clinical proof (that requires prospective
patient/trial data we do not have). Instead it grades each battery-surviving
axis on a three-rung external-validation ladder and assigns a trust tier.

Rungs
-----
Rung 1  Benchmark A -- known-SL recovery: precision/recall of the pipeline
        against a curated *clinically anchored* SL truth set (drug/trial-backed
        pairs). PARP/BRCA is scored separately as modality-limited because
        CRISPR KO != PARP trapping.
Rung 2  Benchmark B -- cross-screen reproducibility: re-run the L1 differential
        dependency test on the independent DepMap per-screen effect matrix
        (ScreenGeneEffect). NOTE: this shares most cell lines with the
        integrated Chronos discovery matrix, so it is an *integration-
        robustness* check, not an independent-cohort replication. Labeled as
        such everywhere.
Rung 3  Clinical actionability + precedent: map each targetable node to real
        drugs + clinical phase (Broad Drug Repurposing Hub), clinical gene-drug
        annotations + FDA/EMA labels (clinpgx / PharmGKB-derived), and (for
        headline axes) literature/trial precedent.

Trust tiers
-----------
T1 Clinically anchored : survives battery + replicates in independent screen
                         + targetable node has an approved/late-phase drug
                         + clinical annotation/precedent.
T2 Reproduced+druggable: survives battery + replicates + druggable partner
                         (any phase), no strong clinical precedent.
T3 Reproduced,not-drugg: survives battery + replicates, no clinical-stage drug.
T4 Discovery-only      : survives battery but fails replication or is single-
                         modality. Explicitly discovery-only.

All statistics reuse sl_battery (same one-sided MWU, Cohen's d, BH-FDR).
Nothing is fabricated; every clinical claim traces to a dataset row or PMID.
"""

from __future__ import annotations
import numpy as np
import pandas as pd

from sl_battery import mw_less, pooled_cohens_d, bh_fdr, LOF_EFFECTS, MIN_N


# ---------------------------------------------------------------------------
# Rung 2 -- cross-screen reproducibility (integration-robustness check)
# ---------------------------------------------------------------------------
def build_independent_screen(screen_csv: str, screenmap_csv: str) -> pd.DataFrame:
    """Load DepMap ScreenGeneEffect (per-screen) and collapse ScreenID->ModelID.

    Multiple screens per model are averaged. Gene columns are de-suffixed of
    the ' (ENTREZ)' tag to match the discovery matrix and mutation frame.
    Returns a ModelID-indexed effect matrix.
    """
    screen = pd.read_csv(screen_csv, index_col=0)
    screen.columns = [c.split(" (")[0] for c in screen.columns]
    smap = pd.read_csv(screenmap_csv)
    s2m = dict(zip(smap["ScreenID"], smap["ModelID"]))
    screen.index = screen.index.map(s2m)
    screen = screen[~screen.index.isna()]
    return screen.groupby(level=0).mean()


def sets_for(muts: pd.DataFrame, universe: set, driver: str, mode: str):
    """Reproduce sl_battery.Battery.sets_for mutant/WT construction.

    LOF  -> LOF effect set.
    ACT  -> driver-annotated missense if >=MIN_N in universe else any missense.
    CN_* -> handled by caller (expression/CN proxy), returns (None, None).
    """
    g = muts[muts["gene_symbol"] == driver]
    if mode == "LOF":
        mut = set(g[g["effect"].isin(LOF_EFFECTS)]["ModelID"])
    elif mode == "ACT":
        drv = set(g[(g["effect"] == "missense") & (g.get("cancer_driver", False) == True)]["ModelID"])
        mut = drv if len(drv & universe) >= MIN_N else set(g[g["effect"] == "missense"]["ModelID"])
    else:
        return None, None
    wt = universe - mut
    return mut, wt


def replicate_axis(indep: pd.DataFrame, muts: pd.DataFrame, driver: str,
                   partner: str, mode: str) -> dict:
    """Re-run L1 (one-sided MWU mutant-more-dependent) on the independent screen.

    Returns dict with n_mut, n_wt, delta, p, cohens_d, and status
    ('insufficient_n' if not enough mutant/WT lines or partner missing).
    CN_* modes are marked 'na_cn' (proxy CN not re-derivable in this matrix).
    """
    out = dict(driver=driver, partner=partner, mode=mode,
               screen_n_mut=0, screen_n_wt=0, screen_delta=np.nan,
               screen_p=np.nan, screen_d=np.nan, screen_status="insufficient_n")
    if mode.startswith("CN"):
        out["screen_status"] = "na_cn"
        return out
    universe = set(indep.index)
    if partner not in indep.columns:
        out["screen_status"] = "partner_absent"
        return out
    mut, wt = sets_for(muts, universe, driver, mode)
    if mut is None:
        out["screen_status"] = "na_mode"
        return out
    mut_u = [x for x in (mut & universe)]
    wt_u = [x for x in (wt & universe)]
    mv = indep.loc[mut_u, partner].dropna()
    wv = indep.loc[wt_u, partner].dropna()
    out["screen_n_mut"], out["screen_n_wt"] = int(len(mv)), int(len(wv))
    if len(mv) < MIN_N or len(wv) < MIN_N:
        return out
    delta, p = mw_less(mv.values, wv.values)
    out["screen_delta"] = float(delta)
    out["screen_p"] = float(p)
    out["screen_d"] = float(pooled_cohens_d(mv.values, wv.values))
    out["screen_status"] = "tested"
    return out


def score_reproducibility(rep_df: pd.DataFrame, disc_df: pd.DataFrame,
                          p_gate: float = 0.05) -> pd.DataFrame:
    """Add replication verdict + cross-screen concordance to per-axis table.

    - screen_replicated: True if tested, screen_p<p_gate, and effect direction
      (sign of screen_delta) matches discovery L1_delta (both negative =
      mutant-more-dependent). BH-FDR is computed across tested axes.
    - sign_concordant: sign(screen_delta)==sign(discovery delta).
    """
    df = rep_df.merge(
        disc_df[["driver", "partner", "mode", "L1_delta", "L1_p", "L1_cohens_d"]],
        on=["driver", "partner", "mode"], how="left")
    tested = df["screen_status"] == "tested"
    df["screen_fdr"] = np.nan
    if tested.any():
        df.loc[tested, "screen_fdr"] = bh_fdr(df.loc[tested, "screen_p"].values)
    df["sign_concordant"] = (
        np.sign(df["screen_delta"]) == np.sign(df["L1_delta"])) & tested
    df["screen_replicated"] = (
        tested & (df["screen_p"] < p_gate) & (df["screen_delta"] < 0)
        & (df["L1_delta"] < 0))
    return df


# ---------------------------------------------------------------------------
# Rung 3 -- clinical actionability + precedent
# ---------------------------------------------------------------------------
_PHASE_RANK = {
    "Launched": 5, "Approved": 5, "Phase 4": 4, "Phase 3": 3,
    "Phase 2/Phase 3": 3, "Phase 2": 2, "Phase 1/Phase 2": 2,
    "Phase 1": 1, "Preclinical": 0, "Withdrawn": 0, "": -1, None: -1,
}


def phase_rank(phase: str) -> int:
    return _PHASE_RANK.get(str(phase).strip(), 0)


def build_hub_target_index(hub: pd.DataFrame) -> dict:
    """Map gene symbol -> list of (drug, clinical_phase, moa) from repurposing hub."""
    idx: dict[str, list] = {}
    for _, r in hub.iterrows():
        tgt = r.get("target")
        if pd.isna(tgt):
            continue
        for g in str(tgt).split("|"):
            g = g.strip()
            if not g:
                continue
            idx.setdefault(g, []).append(
                (r.get("pert_iname"), r.get("clinical_phase"), r.get("moa")))
    return idx


def actionability_for(gene: str, hub_idx: dict, top_k: int = 5) -> dict:
    """Return druggability summary for a gene node from the repurposing hub."""
    hits = hub_idx.get(gene, [])
    if not hits:
        return dict(druggable=False, n_drugs=0, top_clinical_phase="none",
                    top_phase_rank=-1, clinical_drugs="")
    ranked = sorted(hits, key=lambda h: phase_rank(h[1]), reverse=True)
    top_rank = phase_rank(ranked[0][1])
    top_phase = str(ranked[0][1])
    # drugs at the top phase
    top_drugs = [str(h[0]) for h in ranked if phase_rank(h[1]) == top_rank and pd.notna(h[0])]
    return dict(druggable=True, n_drugs=len(hits),
                top_clinical_phase=top_phase, top_phase_rank=top_rank,
                clinical_drugs="; ".join(dict.fromkeys(top_drugs))[:300])


def build_clinpgx_index(rel_tsv: str) -> dict:
    """Gene symbol -> set of associated chemical names from clinpgx relationships."""
    rel = pd.read_csv(rel_tsv, sep="\t",
                      usecols=["Entity1_name", "Entity1_type", "Entity2_name",
                               "Entity2_type", "Association", "PMIDs"])
    idx: dict[str, dict] = {}
    for _, r in rel.iterrows():
        gene = None
        if r["Entity1_type"] == "Gene":
            gene = r["Entity1_name"]
        elif r["Entity2_type"] == "Gene":
            gene = r["Entity2_name"]
        if gene is None:
            continue
        d = idx.setdefault(gene, {"n_annotations": 0, "pmids": set(),
                                  "associated": set()})
        d["n_annotations"] += 1
        if pd.notna(r.get("PMIDs")):
            for p in str(r["PMIDs"]).split(";"):
                p = p.strip()
                if p:
                    d["pmids"].add(p)
        # record the non-gene partner name
        other = r["Entity2_name"] if r["Entity1_type"] == "Gene" else r["Entity1_name"]
        if pd.notna(other):
            d["associated"].add(str(other))
    return idx


def clinpgx_for(gene: str, clin_idx: dict, max_pmids: int = 5) -> dict:
    d = clin_idx.get(gene)
    if not d:
        return dict(clinpgx_annotated=False, clinpgx_n=0, clinpgx_pmids="")
    pmids = list(d["pmids"])[:max_pmids]
    return dict(clinpgx_annotated=True, clinpgx_n=int(d["n_annotations"]),
                clinpgx_pmids="; ".join(pmids))


def build_druglabel_index(bygene_tsv: str) -> set:
    """Set of gene symbols with an FDA/EMA drug label annotation."""
    dl = pd.read_csv(bygene_tsv, sep="\t")
    col = "Gene Symbol" if "Gene Symbol" in dl.columns else dl.columns[1]
    return set(dl[col].dropna().astype(str))


# ---------------------------------------------------------------------------
# Trust-tier assignment (combines the three rungs)
# ---------------------------------------------------------------------------
def assign_trust_tier(row: pd.Series) -> str:
    """T1..T4 per the doctrine. Requires columns:
    screen_replicated (bool), druggable (bool), top_phase_rank (int),
    clinpgx_annotated (bool), has_label (bool), trial_precedent (str/na).
    """
    replicated = bool(row.get("screen_replicated", False))
    druggable = bool(row.get("druggable", False))
    top_rank = int(row.get("top_phase_rank", -1)) if pd.notna(row.get("top_phase_rank")) else -1
    clin = bool(row.get("clinpgx_annotated", False)) or bool(row.get("has_label", False))
    precedent = str(row.get("trial_precedent", "")).strip()
    has_precedent = precedent not in ("", "nan", "none", "no precedent found", "not_searched")

    if not replicated:
        return "T4_discovery_only"
    # replicated from here
    late_phase = top_rank >= 3  # Phase 3+ / approved
    if late_phase and (clin or has_precedent):
        return "T1_clinically_anchored"
    if druggable:
        return "T2_reproduced_druggable"
    return "T3_reproduced_not_druggable"


# ---------------------------------------------------------------------------
# v2.0 clinical-grade layer -- cross-platform independent replication + evidence ladder
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "v2.0-clinical-grade"

# Evidence-grade rubric (ordered, conservative). See CLINICAL_GRADE_REPORT.md.
#   E1 approved-in-context    : approved drug for the node in the matched genomic context
#   E2 clinical-trial precedent: node in active trial (ideally matched context)
#   E3 preclinical in vivo    : SL shown in animal/xenograft/PDX
#   E4 preclinical in vitro / paralog-established
#   E5 computational-only     : pipeline + reproducibility only; no external biology
EVIDENCE_GRADES = ("E1", "E2", "E3", "E4", "E5")


def load_sanger_score(gene_effect_csv: str) -> pd.DataFrame:
    """Load the Sanger Project Score CERES matrix (Figshare 9116732) and return
    a cell-line x gene DataFrame with clean HUGO gene columns (Entrez suffix
    stripped) and Broad ACH- row IDs preserved.

    This is an ORTHOGONAL platform to the Broad Avana/Chronos discovery matrix:
    different library (KY1.0/1.1), wet-lab (Sanger), and algorithm (CERES).
    Note: this Figshare version is pre-filtered to lines with Broad CN data, so
    its cell panel overlaps the Broad discovery panel almost entirely -- it
    supports CROSS-PLATFORM replication, not cross-cohort (held-out) replication.
    """
    df = pd.read_csv(gene_effect_csv, index_col=0)
    df.columns = [c.split(" (")[0] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def independent_replication(sanger: pd.DataFrame, muts: pd.DataFrame,
                            axes: pd.DataFrame) -> pd.DataFrame:
    """Run cross-platform replication for a set of axes on the Sanger matrix.

    `axes` must have columns driver, partner, mode. Reuses replicate_axis (the
    same one-sided MWU + Cohen's d as discovery) and applies BH-FDR across the
    TESTED axes only. Returns a per-axis DataFrame with sanger_* columns,
    sign_concordant (needs L1_cohens_d in `axes`), and an independent_replication_status.
    """
    rows = []
    for _, r in axes.iterrows():
        res = replicate_axis(sanger, muts, r["driver"], r["partner"], r["mode"])
        if "L1_cohens_d" in axes.columns:
            res["L1_cohens_d"] = r["L1_cohens_d"]
        rows.append(res)
    out = pd.DataFrame(rows)
    tested = out["screen_status"] == "tested"
    out["sanger_fdr"] = np.nan
    if tested.any():
        out.loc[tested, "sanger_fdr"] = bh_fdr(out.loc[tested, "screen_p"].values)
    if "L1_cohens_d" in out.columns:
        out["sign_concordant"] = np.sign(out["screen_d"]) == np.sign(out["L1_cohens_d"])
    else:
        out["sign_concordant"] = np.nan
    out["sanger_replicated"] = tested & out["sign_concordant"].fillna(False) & (out["sanger_fdr"] < 0.05)

    def _status(row):
        st = row["screen_status"]
        if st == "na_cn":
            return "not_testable_cn"
        if st == "insufficient_n":
            return "insufficient_n"
        if st == "partner_absent":
            return "partner_absent"
        if row["sanger_replicated"]:
            return "replicated_crossplatform"
        if bool(row.get("sign_concordant", False)):
            return "concordant_not_sig"
        return "discordant"

    out["independent_replication_status"] = out.apply(_status, axis=1)
    return out


def assign_evidence_grade(precedent_class: str) -> str:
    """Map a precedent_class to an ordered evidence grade (E1..E5).

    precedent_class taxonomy:
      clinical_target_in_trials -> E1/E2 (caller refines E1 if approved-in-context)
      established_SL_drug_dev    -> E3
      preclinical_in_vivo        -> E3
      preclinical_in_vitro       -> E4
      mechanistic_no_drug        -> E4 (pathway-plausible) [conservative floor]
      none_found                 -> E5

    This is the conservative default map; approved-in-context E1 promotions are
    assigned explicitly by the curator with a citation (context-match required).
    Promiscuous/hub-only hits cannot raise a grade.
    """
    mapping = {
        "clinical_target_in_trials": "E2",
        "established_SL_drug_dev": "E3",
        "preclinical_in_vivo": "E3",
        "preclinical_in_vitro": "E4",
        "mechanistic_no_drug": "E4",
        "none_found": "E5",
    }
    return mapping.get(precedent_class, "E5")


def assign_trust_tier_v2(evidence_grade: str, replication_status: str) -> str:
    """Two-dimensional trust tier from (evidence_grade, independent replication).

    Keeps external precedent and internal replication as orthogonal axes.
    """
    strong = evidence_grade in ("E1", "E2")
    mid = evidence_grade in ("E3", "E4")
    replicated = replication_status == "replicated_crossplatform"
    if strong and replicated:
        return "T1_anchored_and_replicated"
    if strong:
        return "T2_anchored_precedent"
    if mid and replicated:
        return "T2_anchored_precedent"
    if mid:
        return "T3_precedent_only"
    if replicated:
        return "T3_replicated_no_precedent"
    return "T4_discovery_only"
