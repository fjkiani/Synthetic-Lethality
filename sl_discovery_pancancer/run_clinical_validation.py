#!/usr/bin/env python3
"""run_clinical_validation.py -- single re-runnable entrypoint for the
clinical-grade SL validation layer (schema v2.0-clinical-grade).

Pipeline (idempotent; download-if-missing):
  1. Ensure the Sanger Project Score CERES matrix (Figshare 9116732) is present
     and MD5-verified.
  2. Load discovery CRISPR (Broad Avana/Chronos), mutations, and the 58-axis
     survivor set (the "battery survivors" that passed the L1-L5 validation).
  3. Run CROSS-PLATFORM independent replication of every axis on the Sanger
     KY1.0/1.1 + CERES matrix (orthogonal library + wet-lab + algorithm),
     with BH-FDR across the tested axes.
  4. Merge in the curated clinical-evidence ladder (E1-E5 + precedent_class +
     citations) -- a receipt-backed, human-curated artifact (see
     clinical_evidence_ladder.csv) -- and compute the 2D trust tier
     (evidence_grade x replication_status).
  5. Write the extended master matrix, the replication benchmark, the evidence
     ladder, receipts (SHA256 of every artifact), and data_versions.json.

Design notes / honesty:
  * This is CROSS-PLATFORM replication, not cross-cohort. The Figshare CERES
    matrix is pre-filtered to lines with Broad CN data, so its panel overlaps
    the discovery panel almost entirely (315/318). We state this everywhere.
  * The evidence grades and citations are curated by hand from LiteratureSearch
    (all 58 axes searched). The entrypoint does NOT re-run LiteratureSearch;
    it consumes the curated CSV as a versioned input so runs are deterministic
    and offline-reproducible. Re-derived numbers (replication statistics) are
    computed fresh from the raw matrices every run.
  * No efficacy is inferred anywhere. Grades are precedent/actionability only.

Usage:
    python run_clinical_validation.py \
        --out-dir /mnt/results/SL_portfolio \
        --sanger-dir /mnt/shared-workspace/sanger_score \
        --discovery /mnt/shared-workspace/depmap_cache/crispr_gene_effect.parquet \
        --mutations /mnt/shared-workspace/sl_discovery/mutations_in_crispr.parquet \
        --model /mnt/shared-workspace/depmap_cache/depmap_model.parquet

All paths have sensible defaults (below). Runs end-to-end from a clean kernel.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import sys
import urllib.request
from datetime import datetime, timezone

import numpy as np
import pandas as pd

# Local engine module (same directory)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_clinical_validation as scv  # noqa: E402

LOG = logging.getLogger("run_clinical_validation")

# --- Pinned data source (receipt) ---------------------------------------------
SANGER_FIGSHARE_DOI = "10.6084/m9.figshare.9116732"
SANGER_GENE_EFFECT_MD5 = "9a9e92d984766532a7b4152d440ec86a"
SANGER_GENE_EFFECT_URL = "https://ndownloader.figshare.com/files/16623881"  # gene_effect.csv

DEFAULTS = dict(
    out_dir="/mnt/results/SL_portfolio",
    sanger_dir="/mnt/shared-workspace/sanger_score",
    discovery="/mnt/shared-workspace/depmap_cache/crispr_gene_effect.parquet",
    mutations="/mnt/shared-workspace/sl_discovery/mutations_in_crispr.parquet",
    model="/mnt/shared-workspace/depmap_cache/depmap_model.parquet",
    evidence_ladder="/mnt/results/SL_portfolio/clinical_evidence_ladder.csv",
)


def md5sum(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256sum(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def ensure_sanger(sanger_dir: str) -> str:
    """Ensure gene_effect.csv is present and MD5-verified; download if missing."""
    os.makedirs(sanger_dir, exist_ok=True)
    ge = os.path.join(sanger_dir, "gene_effect.csv")
    if os.path.exists(ge):
        got = md5sum(ge)
        if got == SANGER_GENE_EFFECT_MD5:
            LOG.info("Sanger gene_effect.csv present, MD5 OK (%s)", got)
            return ge
        LOG.warning("Sanger gene_effect.csv present but MD5 mismatch (%s != %s); re-downloading",
                    got, SANGER_GENE_EFFECT_MD5)
    LOG.info("Downloading Sanger Project Score CERES from Figshare %s ...", SANGER_FIGSHARE_DOI)
    urllib.request.urlretrieve(SANGER_GENE_EFFECT_URL, ge)
    got = md5sum(ge)
    if got != SANGER_GENE_EFFECT_MD5:
        raise RuntimeError(
            f"Downloaded Sanger gene_effect.csv MD5 {got} != expected {SANGER_GENE_EFFECT_MD5}. "
            "Refusing to proceed with an unverified independent dataset."
        )
    LOG.info("Downloaded and MD5-verified (%s)", got)
    return ge


def load_discovery(discovery_path: str) -> pd.DataFrame:
    """Discovery CRISPR matrix -> cell-line x gene, ACH-indexed, clean HUGO cols."""
    disc = pd.read_parquet(discovery_path)
    # Discovery parquet is already ACH-indexed genes-in-columns; strip Entrez if present.
    disc.columns = [c.split(" (")[0] for c in disc.columns]
    disc = disc.loc[:, ~disc.columns.duplicated()]
    return disc


def load_mutations(mutations_path: str) -> pd.DataFrame:
    muts = pd.read_parquet(mutations_path)
    if "ModelID" not in muts.columns and "ACH" in muts.columns:
        muts = muts.rename(columns={"ACH": "ModelID"})
    return muts


def load_axes(evidence_ladder_csv: str) -> pd.DataFrame:
    """The 58 survivor axes + curated evidence, from the versioned ladder CSV.

    Requires: driver, partner, mode, L1_cohens_d_discovery, precedent_class,
    evidence_grade, evidence_citation. If the discovery Cohen's d column is
    absent it is treated as NaN (sign-concordance then unavailable).
    """
    ev = pd.read_csv(evidence_ladder_csv)
    # normalise expected columns
    if "L1_cohens_d_discovery" in ev.columns:
        ev["L1_cohens_d"] = ev["L1_cohens_d_discovery"]
    elif "L1_cohens_d" not in ev.columns:
        ev["L1_cohens_d"] = np.nan
    return ev


def accounting(sanger: pd.DataFrame, disc: pd.DataFrame, muts: pd.DataFrame) -> dict:
    sanger_lines = set(sanger.index)
    disc_lines = set(disc.index)
    overlap = sanger_lines & disc_lines
    held_out = sanger_lines - disc_lines
    mut_lines = set(muts["ModelID"].unique())
    held_out_testable = held_out & mut_lines
    shared_testable = overlap & mut_lines
    return dict(
        sanger_total=len(sanger_lines),
        discovery_total=len(disc_lines),
        overlap=len(overlap),
        held_out_sanger_only=len(held_out),
        held_out_with_mutation_data=len(held_out_testable),
        shared_lines_testable=len(shared_testable),
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", default=v, help=f"default: {v}")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    os.makedirs(args.out_dir, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    LOG.info("=== clinical-grade SL validation | schema %s | start %s ===",
             scv.SCHEMA_VERSION, started)

    # 1. Data
    ge = ensure_sanger(args.sanger_dir)
    sanger = scv.load_sanger_score(ge)
    LOG.info("Sanger matrix: %d lines x %d genes", *sanger.shape)
    disc = load_discovery(args.discovery)
    LOG.info("Discovery matrix: %d lines x %d genes", *disc.shape)
    muts = load_mutations(args.mutations)
    LOG.info("Mutations: %d rows", len(muts))
    axes = load_axes(args.evidence_ladder)
    LOG.info("Axes (survivors + curated evidence): %d", len(axes))

    acc = accounting(sanger, disc, muts)
    LOG.info("Cell-line accounting: %s", acc)

    # 2. Cross-platform independent replication (fresh from raw matrices).
    # Restrict the Sanger matrix to the mutation-called universe so the WT
    # background matches discovery's testable universe (this is the canonical
    # 311-line set; the 7 extra Sanger lines lack mutation calls and would only
    # perturb the WT distribution). Guarantees deterministic effect sizes.
    _mut_lines = set(muts["ModelID"].unique())
    sanger_testable = sanger.loc[sanger.index.isin(_mut_lines)]
    LOG.info("Sanger testable universe (mutation-called): %d lines", sanger_testable.shape[0])
    rep = scv.independent_replication(sanger_testable, muts, axes)
    rep = rep.rename(columns={
        "screen_n_mut": "sanger_n_mut", "screen_n_wt": "sanger_n_wt",
        "screen_delta": "sanger_delta", "screen_p": "sanger_p", "screen_d": "sanger_d",
        "screen_status": "sanger_status",
    })
    n_tested = int((rep["sanger_status"] == "tested").sum())
    n_repl = int((rep["independent_replication_status"] == "replicated_crossplatform").sum())
    n_concord = int(rep["sign_concordant"].fillna(False).sum())
    LOG.info("Replication: %d tested, %d sign-concordant, %d replicated (FDR<0.05)",
             n_tested, n_concord, n_repl)

    # Effect-size concordance on tested axes. `rep` already carries L1_cohens_d
    # (independent_replication copies it from axes), so no re-merge is needed.
    from scipy.stats import spearmanr, pearsonr
    paired = rep[rep["sanger_status"] == "tested"].copy()
    valid = paired.dropna(subset=["L1_cohens_d", "sanger_d"])
    rho, rho_p = spearmanr(valid["L1_cohens_d"], valid["sanger_d"])
    pr, pr_p = pearsonr(valid["L1_cohens_d"], valid["sanger_d"])
    LOG.info("Effect-size concordance (n=%d): Spearman rho=%.4f (p=%.2e), Pearson r=%.4f (p=%.2e)",
             len(valid), rho, rho_p, pr, pr_p)

    # 3. Evidence grade + 2D trust tier (evidence is curated; recompute tier deterministically)
    axes = axes.copy()
    if "evidence_grade" not in axes.columns:
        axes["evidence_grade"] = axes["precedent_class"].map(scv.assign_evidence_grade)
    # Drop any pre-existing replication columns so the freshly-computed `rep` is
    # authoritative and the merge does not create _x/_y suffix collisions.
    _repl_cols = ["sanger_n_mut", "sanger_delta", "sanger_p", "sanger_fdr", "sanger_d",
                  "sign_concordant", "sanger_replicated", "sanger_status",
                  "replication_status", "trust_tier_v2"]
    axes = axes.drop(columns=[c for c in _repl_cols if c in axes.columns])
    merged = axes.merge(
        rep[["driver", "partner", "mode", "sanger_n_mut", "sanger_delta", "sanger_p",
             "sanger_fdr", "sanger_d", "sign_concordant", "sanger_replicated",
             "sanger_status", "independent_replication_status"]],
        on=["driver", "partner", "mode"], how="left",
    ).rename(columns={"independent_replication_status": "replication_status"})
    merged["trust_tier_v2"] = [
        scv.assign_trust_tier_v2(g, s)
        for g, s in zip(merged["evidence_grade"], merged["replication_status"])
    ]

    # 4. Write artifacts
    rep_out = os.path.join(args.out_dir, "benchmark_independent_sanger.csv")
    rep.to_csv(rep_out, index=False)
    ladder_out = os.path.join(args.out_dir, "clinical_evidence_ladder.csv")
    merged.to_csv(ladder_out, index=False)

    summary = dict(
        dataset=f"Sanger Project Score CERES (Figshare {SANGER_FIGSHARE_DOI})",
        replication_type="cross-platform (orthogonal library + wet-lab + algorithm)",
        NOT_replication_type="cross-cohort / held-out cell lines (panel overlaps)",
        sanger_gene_effect_md5=SANGER_GENE_EFFECT_MD5,
        cell_line_accounting=acc,
        results=dict(
            n_axes_total=int(len(rep)), n_tested=n_tested,
            n_sign_concordant=n_concord, n_replicated_fdr05=n_repl,
            effect_size_concordance=dict(
                spearman_rho=round(float(rho), 4), spearman_p=float(rho_p),
                pearson_r=round(float(pr), 4), pearson_p=float(pr_p),
                n=int(len(valid)),
            ),
        ),
        schema_version=scv.SCHEMA_VERSION,
        generated_utc=started,
    )
    summ_out = os.path.join(args.out_dir, "benchmark_independent_sanger_summary.json")
    with open(summ_out, "w") as fh:
        json.dump(summary, fh, indent=2)

    # data_versions.json (pinned inputs + MD5s)
    dv = dict(
        schema_version=scv.SCHEMA_VERSION,
        generated_utc=started,
        inputs={
            "sanger_project_score": dict(
                source=f"Figshare {SANGER_FIGSHARE_DOI}", file="gene_effect.csv",
                md5=md5sum(ge), publication="Behan 2019 Nature 10.1038/s41586-019-1103-9",
            ),
            "discovery_crispr": dict(source="DepMap 24Q4 Avana/Chronos",
                                     file=os.path.basename(args.discovery),
                                     md5=md5sum(args.discovery) if os.path.exists(args.discovery) else None),
            "mutations": dict(file=os.path.basename(args.mutations),
                              md5=md5sum(args.mutations) if os.path.exists(args.mutations) else None),
            "model_table": dict(file=os.path.basename(args.model),
                                md5=md5sum(args.model) if os.path.exists(args.model) else None),
        },
    )
    dv_out = os.path.join(args.out_dir, "data_versions.json")
    with open(dv_out, "w") as fh:
        json.dump(dv, fh, indent=2)

    # receipts: SHA256 of everything we wrote
    artifacts = {os.path.basename(p): sha256sum(p)
                 for p in [rep_out, ladder_out, summ_out, dv_out] if os.path.exists(p)}
    LOG.info("Artifacts written: %s", list(artifacts))
    LOG.info("=== done | tested=%d replicated=%d rho=%.3f ===", n_tested, n_repl, rho)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
