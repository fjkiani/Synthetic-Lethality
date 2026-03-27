"""
Drug mapping layer — Step 2 of the SL pipeline.

For each candidate SL partner gene:
  A. DepMap PRISM / GDSC — quantify mutant vs WT drug sensitivity differential
  B. ChEMBL — retrieve known inhibitors, MoA, max clinical phase
  C. OncoKB  — retrieve therapeutic levels of evidence (if token available)

Returns: List[GeneDrugPair] ranked by rank_score.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .config import get_settings
from .models import (
    DepmapDrugEvidence,
    DrugClass,
    ExternalDrugEvidence,
    GeneDrugPair,
    MutationType,
)

logger = logging.getLogger(__name__)
cfg = get_settings()


# ── Drug class inference ──────────────────────────────────────────────────────

_MoA_CLASS_MAP = {
    "parp": DrugClass.PARP_INHIBITOR,
    "atr": DrugClass.ATR_INHIBITOR,
    "atm": DrugClass.ATM_INHIBITOR,
    "wee1": DrugClass.WEE1_INHIBITOR,
    "chk1": DrugClass.CHK1_INHIBITOR,
    "chk2": DrugClass.CHK1_INHIBITOR,
    "kinase": DrugClass.KINASE_INHIBITOR,
    "histone": DrugClass.EPIGENETIC,
    "hdac": DrugClass.EPIGENETIC,
    "bromodomain": DrugClass.EPIGENETIC,
    "proteasome": DrugClass.PROTEASOME,
    "topoisomerase": DrugClass.TOPOISOMERASE,
    "antimetabolite": DrugClass.ANTIMETABOLITE,
    "alkylating": DrugClass.ALKYLATING,
    "immune checkpoint": DrugClass.IMMUNOTHERAPY,
    "pd-1": DrugClass.IMMUNOTHERAPY,
    "pd-l1": DrugClass.IMMUNOTHERAPY,
    "ctla-4": DrugClass.IMMUNOTHERAPY,
}


def _infer_drug_class(text: str) -> DrugClass:
    if not text:
        return DrugClass.UNKNOWN
    t = text.lower()
    for keyword, dc in _MoA_CLASS_MAP.items():
        if keyword in t:
            return dc
    return DrugClass.OTHER


# ── A. PRISM drug-response differential ──────────────────────────────────────

def prism_drug_sensitivity(
    gene: str,
    mutant_ids: List[str],
    wt_ids: List[str],
    prism_df: pd.DataFrame,
    prism_meta: pd.DataFrame,
    top_n: int = 20,
) -> List[DepmapDrugEvidence]:
    """
    Identify PRISM compounds where mutant lines are significantly more sensitive.

    Returns a list of DepmapDrugEvidence sorted by delta_viability (most negative first).
    """
    available = set(prism_df.index)
    mut_avail = [m for m in mutant_ids if m in available]
    wt_avail = [w for w in wt_ids if w in available]

    if len(mut_avail) < 3 or len(wt_avail) < 3:
        logger.warning("Insufficient PRISM coverage for gene %s drug mapping", gene)
        return []

    mut_viab = prism_df.loc[mut_avail]
    wt_viab = prism_df.loc[wt_avail]

    results = []
    for compound in prism_df.columns:
        m = mut_viab[compound].dropna().values
        w = wt_viab[compound].dropna().values
        if len(m) < 3 or len(w) < 3:
            continue
        delta = float(m.mean() - w.mean())
        if delta >= 0:   # we want mutant MORE sensitive (lower viability LFC)
            continue
        stat, pval = stats.mannwhitneyu(m, w, alternative="less")
        results.append({
            "compound": compound,
            "delta_viability": round(delta, 4),
            "p_value": float(pval),
            "n_mut": len(m),
            "n_wt": len(w),
        })

    if not results:
        return []

    p_vals = [r["p_value"] for r in results]
    _, fdr_vals, _, _ = multipletests(p_vals, method="fdr_bh")
    for r, fdr in zip(results, fdr_vals):
        r["fdr"] = round(float(fdr), 6)

    # Filter & sort
    results = [r for r in results if r["fdr"] <= 0.25]
    results.sort(key=lambda x: x["delta_viability"])
    results = results[:top_n]

    # Resolve compound names from metadata
    name_col = next(
        (c for c in prism_meta.columns if c.lower() in ("name", "compound_name", "drug_name")),
        None,
    )
    evidence_list = []
    for r in results:
        drug_name = r["compound"]
        if name_col and r["compound"] in prism_meta.index:
            drug_name = str(prism_meta.loc[r["compound"], name_col])
        evidence_list.append(
            DepmapDrugEvidence(
                dataset="PRISM_repurposing",
                drug_name=drug_name,
                delta_viability=r["delta_viability"],
                p_value=r["p_value"],
                fdr=r["fdr"],
                n_mut=r["n_mut"],
                n_wt=r["n_wt"],
            )
        )
    return evidence_list


# ── B. ChEMBL gene → inhibitor lookup ────────────────────────────────────────

def chembl_drugs_for_gene(gene: str, max_phase_min: int = 0) -> List[ExternalDrugEvidence]:
    """
    Use ChEMBL REST API to find approved/clinical-stage inhibitors for a gene.
    Falls back to empty list on network errors.
    """
    try:
        from chembl_webresource_client.new_client import new_client

        target_res = new_client.target
        results = target_res.filter(
            target_synonym__icontains=gene,
            target_type="SINGLE PROTEIN",
        ).only(["target_chembl_id", "pref_name", "target_type"])

        hits: List[ExternalDrugEvidence] = []
        seen_drugs: set = set()

        for tgt in results[:5]:  # cap at 5 targets to avoid explosion
            tid = tgt.get("target_chembl_id")
            if not tid:
                continue

            activity_res = new_client.activity
            activities = activity_res.filter(
                target_chembl_id=tid,
                standard_type__in=["IC50", "Ki", "Kd", "EC50"],
            ).only(["molecule_chembl_id", "pchembl_value"])

            for act in activities[:30]:
                mol_id = act.get("molecule_chembl_id")
                if not mol_id or mol_id in seen_drugs:
                    continue
                seen_drugs.add(mol_id)

                mol = new_client.molecule.get(mol_id)
                if not mol:
                    continue

                phase = mol.get("max_phase") or 0
                if isinstance(phase, str):
                    try:
                        phase = int(phase)
                    except ValueError:
                        phase = 0
                if phase < max_phase_min:
                    continue

                pref_name = mol.get("pref_name") or mol_id
                moa_data = mol.get("molecule_mechanisms", [])
                moa_text = ""
                if moa_data:
                    moa_text = moa_data[0].get("mechanism_of_action", "") or ""

                hits.append(
                    ExternalDrugEvidence(
                        source="ChEMBL",
                        drug_name=pref_name,
                        drug_chembl_id=mol_id,
                        mechanism_of_action=moa_text or None,
                        drug_class=_infer_drug_class(moa_text),
                        primary_targets=[gene],
                        max_phase=phase,
                    )
                )

        return hits[:20]  # cap

    except Exception as exc:
        logger.warning("ChEMBL lookup failed for %s: %s", gene, exc)
        return []


# ── C. Open KB stack (replaces OncoKB — no token required) ──────────────────
#
# Sources: CIViC (GraphQL/nightly TSV) + CGI (bulk TSV) + JAX-CKB (curated seed)
# Fully open, CC0/public domain, zero licensing friction.

def open_kb_drugs(gene: str, variant: Optional[str] = None) -> List[ExternalDrugEvidence]:
    """
    Query the open KB stack (CIViC + CGI + JAX) for clinical evidence.
    Drop-in replacement for the old oncokb_gene_drugs().
    """
    try:
        from ..kb.kb_engine import oncokb_drop_in
        open_hits = oncokb_drop_in(gene, variant)
        results = []
        for h in open_hits:
            results.append(
                ExternalDrugEvidence(
                    source=h.source,
                    drug_name=h.drug_name,
                    mechanism_of_action=h.mechanism_of_action,
                    drug_class=h.drug_class or DrugClass.UNKNOWN,
                    primary_targets=h.primary_targets or [gene],
                    max_phase=h.max_phase,
                    oncokb_level=None,   # No OncoKB — mars_tier carried in notes
                    notes=(
                        f"Mars tier: {h.mars_tier} | "
                        f"Confidence: {h.confidence} | "
                        f"KBs: {', '.join(h.supporting_kbs)}"
                    ),
                )
            )
        return results[:15]
    except Exception as exc:
        logger.warning("Open KB query failed for %s: %s", gene, exc)
        return []


# Backward-compat alias kept so nothing else breaks
def oncokb_gene_drugs(gene: str) -> List[ExternalDrugEvidence]:
    """Deprecated alias → routes to open KB stack."""
    return open_kb_drugs(gene)


# ── D. GDSC gene → compound lookup (via compound annotation) ─────────────────

def gdsc_drugs_for_gene(
    gene: str,
    gdsc_meta: Optional[pd.DataFrame],
) -> List[ExternalDrugEvidence]:
    """
    Look up GDSC compound metadata for drugs targeting `gene`.
    gdsc_meta should have columns: drug_name, targets, target_pathway, max_conc_um.
    """
    if gdsc_meta is None:
        return []

    target_col = next(
        (c for c in gdsc_meta.columns if c.lower() in ("target", "targets", "gene_targets")),
        None,
    )
    name_col = next(
        (c for c in gdsc_meta.columns if c.lower() in ("drug_name", "name", "compound_name")),
        None,
    )
    pathway_col = next(
        (c for c in gdsc_meta.columns if c.lower() in ("target_pathway", "pathway")),
        None,
    )
    if target_col is None or name_col is None:
        return []

    mask = gdsc_meta[target_col].str.contains(gene, case=False, na=False)
    hits = []
    for _, row in gdsc_meta[mask].iterrows():
        drug_name = str(row[name_col])
        pathway = str(row[pathway_col]) if pathway_col else ""
        hits.append(
            ExternalDrugEvidence(
                source="GDSC_annotation",
                drug_name=drug_name,
                mechanism_of_action=pathway or None,
                drug_class=_infer_drug_class(pathway),
                primary_targets=[gene],
            )
        )
    return hits[:10]


# ── E. Rank score ─────────────────────────────────────────────────────────────

def compute_rank_score(
    sl_delta: float,
    sl_fdr: float,
    drug_delta_viability: Optional[float],
    drug_fdr: Optional[float],
    max_phase: Optional[int],
    oncokb_level: Optional[str],
) -> Tuple[float, Dict[str, float]]:
    """
    Composite rank score [0–1 scale before normalisation].

    Components:
      - sl_score      : (|delta_dep| / 2) * (1 - fdr)
      - drug_score    : |delta_viability| * (1 - drug_fdr) if available
      - druggability  : scaled by max_phase and OncoKB level
    """
    # SL signal
    sl_score = min(abs(sl_delta) / 2.0, 1.0) * (1.0 - min(sl_fdr, 1.0))

    # Drug response differential
    drug_score = 0.0
    if drug_delta_viability is not None and drug_fdr is not None:
        drug_score = min(abs(drug_delta_viability) / 3.0, 1.0) * (1.0 - min(drug_fdr, 1.0))

    # Druggability
    phase_score = 0.0
    if max_phase is not None:
        phase_score = {0: 0.1, 1: 0.3, 2: 0.5, 3: 0.7, 4: 1.0}.get(max_phase, 0.1)

    oncokb_score = 0.0
    if oncokb_level:
        lvl_map = {
            "LEVEL_1": 1.0, "LEVEL_2": 0.8, "LEVEL_3A": 0.6,
            "LEVEL_3B": 0.4, "LEVEL_4": 0.2,
        }
        oncokb_score = lvl_map.get(oncokb_level, 0.0)

    druggability = max(phase_score, oncokb_score)

    # Weighted composite
    rank = 0.45 * sl_score + 0.35 * drug_score + 0.20 * druggability
    components = {
        "sl_signal": round(sl_score, 4),
        "drug_response": round(drug_score, 4),
        "druggability": round(druggability, 4),
    }
    return round(rank, 4), components


# ── F. Main entry point ───────────────────────────────────────────────────────

def map_gene_to_drugs(
    gene: str,
    sl_delta: float,
    sl_fdr: float,
    mutant_ids: List[str],
    wt_ids: List[str],
    prism_df: Optional[pd.DataFrame],
    prism_meta: Optional[pd.DataFrame],
    gdsc_meta: Optional[pd.DataFrame],
    top_n_drugs: int = 10,
) -> List[GeneDrugPair]:
    """
    Full drug-mapping pipeline for one SL partner gene.
    """
    # --- PRISM ---
    prism_evidence: List[DepmapDrugEvidence] = []
    if prism_df is not None and prism_meta is not None:
        prism_evidence = prism_drug_sensitivity(
            gene, mutant_ids, wt_ids, prism_df, prism_meta, top_n=top_n_drugs * 2
        )

    # --- ChEMBL ---
    chembl_evidence = chembl_drugs_for_gene(gene)

    # --- OncoKB ---
    oncokb_evidence = oncokb_gene_drugs(gene)

    # --- GDSC meta ---
    gdsc_evidence = gdsc_drugs_for_gene(gene, gdsc_meta)

    # Merge external evidence: prefer OncoKB > ChEMBL > GDSC
    ext_by_name: Dict[str, ExternalDrugEvidence] = {}
    for ev_list in [gdsc_evidence, chembl_evidence, oncokb_evidence]:
        for ev in ev_list:
            key = ev.drug_name.lower()
            ext_by_name[key] = ev  # later sources (higher quality) overwrite

    # Build GeneDrugPair for each PRISM hit, augmented by external evidence
    pairs: List[GeneDrugPair] = []
    seen_drugs: set = set()

    for prism_ev in prism_evidence:
        key = prism_ev.drug_name.lower()
        ext_ev = ext_by_name.get(key)
        max_phase = ext_ev.max_phase if ext_ev else None
        oncokb_lvl = ext_ev.oncokb_level if ext_ev else None
        drug_class = ext_ev.drug_class if ext_ev else DrugClass.UNKNOWN

        rank, components = compute_rank_score(
            sl_delta=sl_delta,
            sl_fdr=sl_fdr,
            drug_delta_viability=prism_ev.delta_viability,
            drug_fdr=prism_ev.fdr,
            max_phase=max_phase,
            oncokb_level=oncokb_lvl,
        )
        pairs.append(
            GeneDrugPair(
                partner_gene=gene,
                drug_name=prism_ev.drug_name,
                drug_class=drug_class,
                depmap_evidence=prism_ev,
                external_evidence=ext_ev,
                rank_score=rank,
                rank_components=components,
            )
        )
        seen_drugs.add(key)

    # Add external-only entries (no PRISM hit, but known inhibitor with high precedence)
    for name_key, ext_ev in ext_by_name.items():
        if name_key in seen_drugs:
            continue
        if (ext_ev.max_phase or 0) < 2 and not ext_ev.oncokb_level:
            continue   # only include external-only if clinical stage or OncoKB evidence
        rank, components = compute_rank_score(
            sl_delta=sl_delta,
            sl_fdr=sl_fdr,
            drug_delta_viability=None,
            drug_fdr=None,
            max_phase=ext_ev.max_phase,
            oncokb_level=ext_ev.oncokb_level,
        )
        pairs.append(
            GeneDrugPair(
                partner_gene=gene,
                drug_name=ext_ev.drug_name,
                drug_class=ext_ev.drug_class,
                depmap_evidence=None,
                external_evidence=ext_ev,
                rank_score=rank,
                rank_components=components,
            )
        )

    pairs.sort(key=lambda x: -x.rank_score)
    return pairs[:top_n_drugs]
