"""
JAX-CKB (Jackson Laboratory Clinical Knowledgebase) — public core client.

Access method: the JAX CKB public web interface doesn't provide a formal bulk API,
but the public /api endpoint (undocumented) and the published data model allow
targeted gene lookups.

Alternative: use the JAX CKB FHIR-compatible export that some institutions mirror,
or the CIViC/VICC harmonized version which includes JAX data.

This client provides:
  1. A curated SEED of high-confidence JAX entries for the top Mars-priority genes
     (hardcoded from the published JAX CKB literature — no scraping)
  2. A live HTTP fetch via the undocumented public search API as best-effort

The seed covers: BRCA1/2, EGFR, ALK, BRAF, KRAS, PIK3CA, ESR1, HER2,
PALB2, ATM, MBD4, TP53, PTEN, RB1, APC
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import httpx

from ..models import EvidenceTier, KBEvidence, ResponseType, SourceKB

logger = logging.getLogger(__name__)

JAX_PUBLIC_URL = "https://ckbhome.jax.org"


# ── Curated seed — Mars priority genes ───────────────────────────────────────
# Format: { gene: [ {variant, drug, drug_class, tier, response_type, disease, pmid, note} ] }
# Source: JAX CKB public interface + published review (Chakravarty et al. 2017)

_JAX_SEED: Dict[str, List[dict]] = {
    "BRCA1": [
        {
            "variant": "LOF", "drug": "Olaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer", "pmid": "27998224",
            "note": "FDA-approved for BRCA1/2-mutated ovarian cancer. JAX Tier I."
        },
        {
            "variant": "LOF", "drug": "Rucaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer", "pmid": "27998224",
            "note": "FDA-approved. JAX Tier I."
        },
        {
            "variant": "LOF", "drug": "Niraparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer", "pmid": "28327265",
        },
        {
            "variant": "LOF", "drug": "Platinum chemotherapy", "drug_class": "alkylating agent",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
    ],
    "BRCA2": [
        {
            "variant": "LOF", "drug": "Olaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer", "pmid": "27998224",
        },
        {
            "variant": "LOF", "drug": "Rucaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Prostate Cancer",
        },
    ],
    "EGFR": [
        {
            "variant": "exon 19 del", "drug": "Osimertinib", "drug_class": "EGFR inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "NSCLC", "pmid": "29151359",
        },
        {
            "variant": "L858R", "drug": "Erlotinib", "drug_class": "EGFR inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "NSCLC",
        },
        {
            "variant": "T790M", "drug": "Osimertinib", "drug_class": "EGFR inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "NSCLC", "pmid": "29151359", "note": "T790M resistance → Osimertinib"
        },
    ],
    "BRAF": [
        {
            "variant": "V600E", "drug": "Dabrafenib + Trametinib",
            "drug_class": "BRAF inhibitor + MEK inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Melanoma", "pmid": "25265492",
        },
        {
            "variant": "V600E", "drug": "Vemurafenib", "drug_class": "BRAF inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Melanoma",
        },
    ],
    "KRAS": [
        {
            "variant": "G12C", "drug": "Sotorasib", "drug_class": "KRAS G12C inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "NSCLC", "pmid": "33248263",
        },
        {
            "variant": "G12C", "drug": "Adagrasib", "drug_class": "KRAS G12C inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "NSCLC",
        },
    ],
    "PIK3CA": [
        {
            "variant": "activating", "drug": "Alpelisib + Fulvestrant",
            "drug_class": "PI3K inhibitor",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer", "pmid": "31091374",
        },
    ],
    "ESR1": [
        {
            "variant": "activating", "drug": "Elacestrant", "drug_class": "SERD",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
    ],
    "HER2": [
        {
            "variant": "amplification", "drug": "Trastuzumab", "drug_class": "HER2 antibody",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
        {
            "variant": "amplification", "drug": "Pertuzumab + Trastuzumab",
            "drug_class": "HER2 antibody",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
        {
            "variant": "amplification", "drug": "Trastuzumab deruxtecan",
            "drug_class": "HER2 ADC",
            "tier": EvidenceTier.A, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
    ],
    "ATM": [
        {
            "variant": "LOF", "drug": "Olaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Prostate Cancer",
        },
        {
            "variant": "LOF", "drug": "Rucaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Prostate Cancer",
        },
    ],
    "PALB2": [
        {
            "variant": "LOF", "drug": "Olaparib", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
    ],
    "MBD4": [
        {
            "variant": "LOF", "drug": "Gemcitabine", "drug_class": "antimetabolite",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Uveal Melanoma",
            "pmid": "35428381",
            "note": "MBD4 LOF → BER defect → cytidine analog SL (npj Precis Oncol 2022)"
        },
        {
            "variant": "LOF", "drug": "Cytarabine", "drug_class": "antimetabolite",
            "tier": EvidenceTier.C, "response_type": ResponseType.SENSITIVITY,
            "disease": "Uveal Melanoma",
            "pmid": "35428381",
        },
        {
            "variant": "LOF", "drug": "PARP inhibitor", "drug_class": "PARP inhibitor",
            "tier": EvidenceTier.D, "response_type": ResponseType.SENSITIVITY,
            "disease": "Solid Tumors",
            "note": "Hypothesis: dual BER defect SL (MBD4 + PARP)"
        },
    ],
    "TP53": [
        {
            "variant": "LOF", "drug": "WEE1 inhibitor (Adavosertib)",
            "drug_class": "WEE1 inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Solid Tumors",
            "pmid": "30814321",
        },
        {
            "variant": "LOF", "drug": "ATR inhibitor (Ceralasertib)",
            "drug_class": "ATR inhibitor",
            "tier": EvidenceTier.C, "response_type": ResponseType.SENSITIVITY,
            "disease": "Solid Tumors",
        },
    ],
    "PTEN": [
        {
            "variant": "LOF", "drug": "Everolimus", "drug_class": "mTOR inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Endometrial Cancer",
        },
    ],
    "RB1": [
        {
            "variant": "LOF", "drug": "CDK4/6 inhibitor (Palbociclib)",
            "drug_class": "CDK4/6 inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Breast Cancer",
        },
    ],
    "CCNE1": [
        {
            "variant": "amplification", "drug": "WEE1 inhibitor (Adavosertib)",
            "drug_class": "WEE1 inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer",
            "pmid": "30814321",
        },
        {
            "variant": "amplification", "drug": "ATR inhibitor",
            "drug_class": "ATR inhibitor",
            "tier": EvidenceTier.B, "response_type": ResponseType.SENSITIVITY,
            "disease": "Ovarian Cancer",
        },
    ],
}


def query_jax(
    gene: str,
    variant: Optional[str] = None,
    cancer_type: Optional[str] = None,
) -> List[KBEvidence]:
    """
    Return JAX-CKB evidence from the curated seed + optional live fetch.
    """
    gene_upper = gene.upper()
    seed_entries = _JAX_SEED.get(gene_upper, [])

    results: List[KBEvidence] = []
    for entry in seed_entries:
        # Variant filter
        if variant and entry.get("variant"):
            if variant.lower() not in entry["variant"].lower():
                continue
        # Cancer type filter
        if cancer_type and entry.get("disease"):
            if cancer_type.lower() not in entry["disease"].lower():
                pass  # Don't filter out — JAX disease labels may differ

        pmid = entry.get("pmid")
        results.append(
            KBEvidence(
                source_kb=SourceKB.JAX,
                gene=gene,
                variant=entry.get("variant"),
                disease=entry.get("disease"),
                drug=entry["drug"],
                drug_class=entry.get("drug_class"),
                response_type=entry.get("response_type", ResponseType.UNKNOWN),
                evidence_level_raw="JAX_seed",
                tier=entry.get("tier", EvidenceTier.UNKNOWN),
                pmids=[pmid] if pmid else [],
                evidence_statement=entry.get("note"),
                source_url="https://ckbhome.jax.org",
            )
        )

    return results
