"""
Tests for the open KB stack — fully offline with mocked HTTP.
"""
from __future__ import annotations

from typing import List
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from sl_agent.kb.clients.cgi_client import _map_tier as cgi_map_tier, _map_response
from sl_agent.kb.clients.civic_client import _map_tier as civic_map_tier, _map_significance
from sl_agent.kb.clients.jax_client import query_jax
from sl_agent.kb.evidence_fuser import fuse_evidence
from sl_agent.kb.models import (
    EvidenceTier,
    KBEvidence,
    KBQueryInput,
    ResponseType,
    SourceKB,
)


# ── Tier mapping ──────────────────────────────────────────────────────────────

def test_civic_tier_A():
    assert civic_map_tier("A") == EvidenceTier.A

def test_civic_tier_B():
    assert civic_map_tier("B") == EvidenceTier.B

def test_civic_tier_C():
    assert civic_map_tier("C") == EvidenceTier.C

def test_civic_tier_D():
    assert civic_map_tier("D") == EvidenceTier.C   # D = case study → Tier C

def test_civic_tier_E():
    assert civic_map_tier("E") == EvidenceTier.D

def test_cgi_tier_fda():
    assert cgi_map_tier("FDA approved") == EvidenceTier.A

def test_cgi_tier_phase3():
    assert cgi_map_tier("Clinical trials Phase III") == EvidenceTier.B

def test_cgi_tier_phase1():
    assert cgi_map_tier("Clinical trials Phase I") == EvidenceTier.C

def test_cgi_tier_preclinical():
    assert cgi_map_tier("Pre-clinical") == EvidenceTier.D


# ── Response type mapping ─────────────────────────────────────────────────────

def test_civic_significance_sensitivity():
    assert _map_significance("Sensitivity/Response") == ResponseType.SENSITIVITY

def test_civic_significance_resistance():
    assert _map_significance("Resistance") == ResponseType.RESISTANCE

def test_cgi_response_sensitive():
    assert _map_response("Responsive") == ResponseType.SENSITIVITY

def test_cgi_response_resistant():
    assert _map_response("Resistant") == ResponseType.RESISTANCE


# ── JAX seed ─────────────────────────────────────────────────────────────────

def test_jax_brca1_returns_entries():
    results = query_jax("BRCA1")
    assert len(results) >= 3
    drugs = [r.drug for r in results]
    assert any("Olaparib" in d for d in drugs)

def test_jax_mbd4_has_gemcitabine():
    results = query_jax("MBD4")
    drugs = [r.drug for r in results]
    assert any("Gemcitabine" in d or "gemcitabine" in d.lower() for d in drugs)

def test_jax_mbd4_tier():
    results = query_jax("MBD4")
    gem = next((r for r in results if "Gemcitabine" in r.drug), None)
    assert gem is not None
    assert gem.tier == EvidenceTier.B

def test_jax_unknown_gene_returns_empty():
    results = query_jax("FAKEGENE123")
    assert results == []

def test_jax_cancer_type_doesnt_exclude():
    # JAX doesn't hard-filter disease labels — should still return entries
    results = query_jax("BRCA1", cancer_type="Ovarian Cancer")
    assert len(results) >= 1

def test_jax_source_is_jax():
    results = query_jax("BRCA1")
    assert all(r.source_kb == SourceKB.JAX for r in results)


# ── Evidence fuser ────────────────────────────────────────────────────────────

def _make_evidence(gene, drug, tier, kb, response=ResponseType.SENSITIVITY) -> KBEvidence:
    return KBEvidence(
        source_kb=kb,
        gene=gene,
        drug=drug,
        tier=tier,
        response_type=response,
        evidence_level_raw="test",
    )


def test_fuser_deduplicates_same_drug():
    """Same drug from 2 KBs → 1 recommendation with 2 supporting KBs."""
    ev = [
        _make_evidence("BRCA1", "Olaparib", EvidenceTier.A, SourceKB.CIVIC),
        _make_evidence("BRCA1", "Olaparib", EvidenceTier.B, SourceKB.JAX),
    ]
    q = KBQueryInput(gene="BRCA1")
    recs = fuse_evidence(ev, q)
    assert len(recs) == 1
    assert len(recs[0].supporting_kbs) == 2
    assert recs[0].tier == EvidenceTier.A   # takes best tier


def test_fuser_consensus_upgrades_tier():
    """Two KBs with Tier C → fuser upgrades to Tier B."""
    ev = [
        _make_evidence("TP53", "AZD1775", EvidenceTier.C, SourceKB.CIVIC),
        _make_evidence("TP53", "AZD1775", EvidenceTier.C, SourceKB.CGI),
    ]
    q = KBQueryInput(gene="TP53")
    recs = fuse_evidence(ev, q)
    assert recs[0].tier == EvidenceTier.B


def test_fuser_confidence_high_for_tier_a():
    ev = [_make_evidence("BRCA1", "Olaparib", EvidenceTier.A, SourceKB.CIVIC)]
    q = KBQueryInput(gene="BRCA1")
    recs = fuse_evidence(ev, q)
    assert recs[0].confidence == "high"


def test_fuser_sorts_by_tier():
    ev = [
        _make_evidence("BRCA1", "Preclinical_Drug", EvidenceTier.D, SourceKB.CGI),
        _make_evidence("BRCA1", "Olaparib", EvidenceTier.A, SourceKB.CIVIC),
        _make_evidence("BRCA1", "Investigational_Drug", EvidenceTier.C, SourceKB.JAX),
    ]
    q = KBQueryInput(gene="BRCA1")
    recs = fuse_evidence(ev, q)
    tiers = [r.tier for r in recs]
    assert tiers[0] == EvidenceTier.A   # highest tier first


def test_fuser_filters_resistance_when_excluded():
    ev = [
        _make_evidence("BRCA1", "Olaparib", EvidenceTier.A, SourceKB.CIVIC, ResponseType.SENSITIVITY),
        _make_evidence("BRCA1", "Rucaparib", EvidenceTier.B, SourceKB.JAX, ResponseType.RESISTANCE),
    ]
    q = KBQueryInput(gene="BRCA1", include_resistance=False)
    recs = fuse_evidence(ev, q)
    assert all(r.response_type != ResponseType.RESISTANCE for r in recs)


def test_fuser_min_tier_filter():
    ev = [
        _make_evidence("BRCA1", "Olaparib", EvidenceTier.A, SourceKB.CIVIC),
        _make_evidence("BRCA1", "Preclinical_X", EvidenceTier.D, SourceKB.CGI),
    ]
    q = KBQueryInput(gene="BRCA1", min_tier=EvidenceTier.C)
    recs = fuse_evidence(ev, q)
    assert all(r.tier not in (EvidenceTier.D, EvidenceTier.UNKNOWN) for r in recs)


def test_fuser_empty_input():
    q = KBQueryInput(gene="BRCA1")
    recs = fuse_evidence([], q)
    assert recs == []


# ── KB engine drop-in (with mocked clients) ───────────────────────────────────

def test_oncokb_drop_in_returns_open_evidence():
    """oncokb_drop_in() should return OpenKBEvidence objects for BRCA1."""
    # Mock both CIViC and CGI to return empty (JAX seed still works)
    with patch("sl_agent.kb.kb_engine.query_civic", return_value=[]), \
         patch("sl_agent.kb.kb_engine.query_cgi", return_value=[]):
        from sl_agent.kb.kb_engine import oncokb_drop_in
        hits = oncokb_drop_in("BRCA1")

    # JAX seed should still return results
    assert len(hits) >= 1
    assert all(hasattr(h, "drug_name") for h in hits)
    assert all(hasattr(h, "max_phase") for h in hits)
    assert all(hasattr(h, "mars_tier") for h in hits)
    # No OncoKB level — this is the open stack
    assert all(h.oncokb_level is None for h in hits)


def test_open_kb_drugs_in_drug_mapper():
    """Verify open_kb_drugs() in drug_mapper.py returns ExternalDrugEvidence objects."""
    with patch("sl_agent.kb.kb_engine.query_civic", return_value=[]), \
         patch("sl_agent.kb.kb_engine.query_cgi", return_value=[]):
        from sl_agent.core.drug_mapper import open_kb_drugs
        results = open_kb_drugs("BRCA1")

    assert len(results) >= 1
    from sl_agent.core.models import ExternalDrugEvidence
    assert all(isinstance(r, ExternalDrugEvidence) for r in results)
    # mars_tier should appear in notes
    assert all("Mars tier" in (r.notes or "") for r in results)
