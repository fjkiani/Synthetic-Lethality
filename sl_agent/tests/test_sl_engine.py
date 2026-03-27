"""
Unit tests for the SL computation engine — fully offline (synthetic data).
"""
import numpy as np
import pandas as pd
import pytest

from sl_agent.core.models import MutationType, SLQueryInput
from sl_agent.core.sl_engine import MIN_GROUP_SIZE, SLEngine


# ── Fixtures ──────────────────────────────────────────────────────────────────

N_LINES = 120
N_GENES = 50
RNG = np.random.default_rng(42)


def _make_crispr(n_lines=N_LINES, n_genes=N_GENES, seed=42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    model_ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    genes = [f"GENE{i}" for i in range(n_genes)]
    data = rng.normal(-0.4, 0.4, (n_lines, n_genes))   # tighter WT baseline
    # Inject a synthetic SL signal for GENE1: mutant lines only
    # WT lines (20-119) stay at baseline → overall median stays above pan-essential threshold
    data[:20, 1] -= 1.5   # first 20 = mutant group: strongly more dependent
    return pd.DataFrame(data, index=model_ids, columns=genes)


def _make_sample_info(n_lines=N_LINES) -> pd.DataFrame:
    model_ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    lineages = ["Ovarian Cancer"] * 60 + ["Lung Cancer"] * 60
    return pd.DataFrame(
        {"OncotreeLineage": lineages, "OncotreePrimaryDisease": lineages},
        index=model_ids,
    )


@pytest.fixture
def engine():
    crispr = _make_crispr()
    sample_info = _make_sample_info()
    return SLEngine(
        crispr_df=crispr,
        sample_info_df=sample_info,
        depmap_release="test_v1",
    )


@pytest.fixture
def basic_query():
    return SLQueryInput(
        gene="GENE0",
        mutation_type=MutationType.ANY,
        cancer_type=None,
        top_n_partners=10,
        top_n_drugs=5,
        fdr_cutoff=0.25,
        delta_dep_cutoff=0.10,
    )


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_engine_initialises(engine):
    assert engine.crispr.shape[0] == N_LINES
    assert engine.crispr.shape[1] == N_GENES


def test_pan_essential_computed(engine):
    pan_ess = engine._compute_pan_essential()
    assert isinstance(pan_ess, set)
    # In synthetic random data no gene has 90%+ dependency frequency — that's expected.
    # With real DepMap data this will be non-empty (core fitness genes).
    # Just validate the return type and that it doesn't blow up.
    assert pan_ess is not None


def test_compute_sl_partners_returns_results(engine, basic_query):
    mutant_ids = [f"ACH-{i:06d}" for i in range(20)]
    wt_ids = [f"ACH-{i:06d}" for i in range(20, 100)]

    ctx, partners, warnings = engine.compute_sl_partners(
        query=basic_query,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )
    assert ctx.n_mut_lines == 20
    assert ctx.n_wt_lines == 80
    assert len(partners) > 0


def test_synthetic_sl_signal_detected(engine, basic_query):
    """GENE1 should be ranked highly as an SL partner of GENE0."""
    mutant_ids = [f"ACH-{i:06d}" for i in range(20)]
    wt_ids = [f"ACH-{i:06d}" for i in range(20, 100)]

    _, partners, _ = engine.compute_sl_partners(
        query=basic_query,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )
    partner_genes = [p.gene for p in partners]
    assert "GENE1" in partner_genes, "Expected synthetic SL signal for GENE1"


def test_delta_dependency_negative(engine, basic_query):
    """All returned partners should have negative delta (mutant more essential)."""
    mutant_ids = [f"ACH-{i:06d}" for i in range(20)]
    wt_ids = [f"ACH-{i:06d}" for i in range(20, 100)]

    _, partners, _ = engine.compute_sl_partners(
        query=basic_query,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )
    for p in partners:
        assert p.delta_dependency < 0, f"{p.gene} delta={p.delta_dependency} should be < 0"


def test_fdr_within_cutoff(engine, basic_query):
    mutant_ids = [f"ACH-{i:06d}" for i in range(20)]
    wt_ids = [f"ACH-{i:06d}" for i in range(20, 100)]

    _, partners, _ = engine.compute_sl_partners(
        query=basic_query,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )
    for p in partners:
        assert p.fdr <= basic_query.fdr_cutoff


def test_too_few_mutants_raises(engine, basic_query):
    mutant_ids = [f"ACH-{i:06d}" for i in range(2)]   # < MIN_GROUP_SIZE
    wt_ids = [f"ACH-{i:06d}" for i in range(2, 50)]
    with pytest.raises(ValueError):
        engine.compute_sl_partners(
            query=basic_query,
            mutant_ids=mutant_ids,
            wt_ids=wt_ids,
        )


def test_codependency_computed(engine, basic_query):
    basic_query_with_codep = basic_query.model_copy(update={"include_codependency": True})
    mutant_ids = [f"ACH-{i:06d}" for i in range(20)]
    wt_ids = [f"ACH-{i:06d}" for i in range(20, 100)]

    _, partners, _ = engine.compute_sl_partners(
        query=basic_query_with_codep,
        mutant_ids=mutant_ids,
        wt_ids=wt_ids,
    )
    # At least some partners should have codependency populated
    codep_set = [p.codependency_r for p in partners if p.codependency_r is not None]
    assert len(codep_set) > 0
