"""
Session-cached live literature queries for BrM genes.

Results populate BrMGeneRow.live_literature_note ONLY.
Not included in frozen artifact — frozen receipts are authoritative.
Returns None on network failure (non-blocking).
"""

from __future__ import annotations

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_CACHE: Dict[str, Optional[str]] = {}


def query_brm_literature(
    gene: str,
    cancer_type: str = "NSCLC",
    max_results: int = 3,
) -> Optional[str]:
    """
    Query PubMed/Semantic Scholar for BrM-relevant papers on a gene.

    Returns a one-paragraph summary or None on failure.
    Session-cached: repeated calls for the same gene return cached result.
    """
    cache_key = f"{gene.upper()}::{cancer_type}"
    if cache_key in _CACHE:
        return _CACHE[cache_key]

    try:
        result = _fetch_literature(gene, cancer_type, max_results)
        _CACHE[cache_key] = result
        return result
    except Exception as exc:
        logger.warning("live_literature query failed for %s: %s", gene, exc)
        _CACHE[cache_key] = None
        return None


def _fetch_literature(gene: str, cancer_type: str, max_results: int) -> Optional[str]:
    """Internal fetch — tries Semantic Scholar API."""
    import urllib.request
    import urllib.parse
    import json

    query = f"{gene} brain metastasis {cancer_type} therapeutic target"
    encoded = urllib.parse.quote(query)
    url = (
        f"https://api.semanticscholar.org/graph/v1/paper/search"
        f"?query={encoded}&limit={max_results}"
        f"&fields=title,year,abstract"
    )

    req = urllib.request.Request(url, headers={"User-Agent": "sl-agent/2.0"})
    with urllib.request.urlopen(req, timeout=8) as resp:
        data = json.loads(resp.read())

    papers = data.get("data", [])
    if not papers:
        return None

    lines = []
    for p in papers[:max_results]:
        title = p.get("title", "")
        year  = p.get("year", "")
        abstr = (p.get("abstract") or "")[:200]
        lines.append(f"[{year}] {title}. {abstr}...")

    return " | ".join(lines)


def clear_cache() -> None:
    """Clear the session cache (for testing)."""
    _CACHE.clear()
