"""
ClinVar client — somatic oncogenicity classifications.

Used as a FILTER layer (not a therapy KB):
  - Identifies variants as pathogenic / likely pathogenic / VUS / benign
  - Somatic oncogenicity annotations (since ClinVar 2023)
  - Source: NCBI FTP + Entrez eutils API

We use the ClinVar VCF/TSV summary for bulk access and eutils for per-variant queries.
This gives us:
  - Variant oncogenicity confidence (to DOWN-WEIGHT benign hits in SL analysis)
  - Germline vs somatic classification
  - Cross-lab consensus (agreement_count → confidence)
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import httpx

logger = logging.getLogger(__name__)

CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CLINVAR_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


# Pathogenicity → filtering weight (0 = exclude, 1 = include with full weight)
CLINVAR_PATHOGENICITY_WEIGHT: Dict[str, float] = {
    "pathogenic": 1.0,
    "likely pathogenic": 0.9,
    "pathogenic/likely pathogenic": 0.95,
    "uncertain significance": 0.5,
    "likely benign": 0.1,
    "benign": 0.0,
    "benign/likely benign": 0.05,
    "conflicting interpretations": 0.6,
    "oncogenic": 1.0,
    "likely oncogenic": 0.85,
    "uncertain oncogenic potential": 0.5,
}


def get_variant_oncogenicity(gene: str, variant: Optional[str] = None) -> Dict[str, float]:
    """
    Query ClinVar eutils for oncogenicity/pathogenicity of a gene/variant.
    Returns a dict: {variant_description → pathogenicity_weight}

    Used to filter SL hits: if a variant is 'Benign' in ClinVar,
    we flag it as lower confidence.
    """
    query_term = f"{gene}[gene] AND (cancer[title] OR oncogenic[title])"
    if variant:
        query_term = f"{gene}[gene] AND {variant}[variant_name]"

    try:
        # Step 1: eSearch
        search_resp = httpx.get(
            CLINVAR_ESEARCH_URL,
            params={
                "db": "clinvar",
                "term": query_term,
                "retmax": 20,
                "retmode": "json",
            },
            timeout=15,
        )
        search_resp.raise_for_status()
        ids = search_resp.json().get("esearchresult", {}).get("idlist", [])

        if not ids:
            return {}

        # Step 2: eSummary
        summary_resp = httpx.get(
            CLINVAR_ESUMMARY_URL,
            params={
                "db": "clinvar",
                "id": ",".join(ids),
                "retmode": "json",
            },
            timeout=15,
        )
        summary_resp.raise_for_status()
        result = summary_resp.json().get("result", {})

        weights: Dict[str, float] = {}
        for uid in ids:
            rec = result.get(uid, {})
            clinical_sig = rec.get("clinical_significance", {})
            desc = clinical_sig.get("description", "").lower().strip()
            title = rec.get("title", uid)
            weight = CLINVAR_PATHOGENICITY_WEIGHT.get(desc, 0.5)
            weights[title] = weight

        return weights

    except Exception as exc:
        logger.warning("ClinVar query failed for %s: %s", gene, exc)
        return {}


def is_likely_benign(gene: str, variant: Optional[str] = None) -> bool:
    """
    Returns True if ClinVar consensus calls this variant Benign/Likely Benign.
    Used to down-rank SL partner suggestions with benign variant context.
    """
    weights = get_variant_oncogenicity(gene, variant)
    if not weights:
        return False
    avg = sum(weights.values()) / len(weights)
    return avg < 0.15
