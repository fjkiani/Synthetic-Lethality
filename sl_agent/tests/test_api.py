"""
Integration tests for the FastAPI endpoints (offline — mocks data layer).
"""
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock

from sl_agent.api.app import create_app
from sl_agent.core.models import (
    AnalysisScope, InputContext, MutationType,
    SLMapResult, SLPartner, CrossValidation, FrameworkSupport,
)


# ── Shared mock result ────────────────────────────────────────────────────────

def _mock_sl_result() -> SLMapResult:
    return SLMapResult(
        input_context=InputContext(
            query_gene="BRCA1",
            mutation_type=MutationType.LOF,
            cancer_type="Ovarian Cancer",
            scope=AnalysisScope.CANCER_SPECIFIC,
            n_total_lines=80,
            n_mut_lines=20,
            n_wt_lines=60,
            depmap_release="24Q4",
        ),
        sl_partners=[
            SLPartner(
                gene="PARP1",
                delta_dependency=-1.2,
                p_value=1e-6,
                fdr=0.001,
                n_mut=20,
                n_wt=60,
                test_type="wilcoxon_ranksum_one_sided",
                effect_size_cohend=-2.1,
                pathway=["DNA_Damage_Response"],
                supporting_frameworks=FrameworkSupport(SLIdR=True, SLAYER=True, literature=True),
            )
        ],
        gene_drug_pairs=[],
        cross_validation=CrossValidation(
            frameworks_checked=["SLIdR", "SLAYER", "SL_RFM_literature"],
            confirmed_pairs=["PARP1"],
            speculative_pairs=[],
        ),
    )


@pytest.fixture
def client():
    app = create_app()
    # Patch DataStore and GDSC loaders so no network I/O happens at startup.
    # Also patch the GDSC warmup loaders so the lifespan 60-s wait is instant.
    with patch("sl_agent.api.app.DataStore.ensure_loaded"), \
         patch("sl_agent.api.routes.DataStore.ensure_loaded"), \
         patch("sl_agent.api.routes.DataStore._crispr", new=pd.DataFrame()), \
         patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc1",
               return_value=pd.DataFrame()), \
         patch("sl_agent.data.gdsc_biomarker_loader.load_gdsc2",
               return_value=pd.DataFrame()):
        with TestClient(app) as c:
            yield c


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_health_endpoint(client):
    resp = client.get("/api/v1/health")
    assert resp.status_code == 200
    data = resp.json()
    assert "status" in data
    assert "version" in data


def test_analyze_returns_200_with_mock(client):
    with patch("sl_agent.api.routes.run_sl_analysis", return_value=_mock_sl_result()):
        resp = client.post(
            "/api/v1/analyze",
            json={
                "gene": "BRCA1",
                "mutation_type": "loss_of_function",
                "cancer_type": "Ovarian Cancer",
                "top_n_partners": 5,
                "top_n_drugs": 3,
            },
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "success"
    assert "result" in data
    assert data["result"]["sl_partners"][0]["gene"] == "PARP1"


def test_analyze_invalid_gene_raises_400(client):
    # Patch run_sl_analysis to raise ValueError
    with patch(
        "sl_agent.api.routes.run_sl_analysis",
        side_effect=ValueError("Gene FAKEGENE not found"),
    ):
        resp = client.post(
            "/api/v1/analyze",
            json={"gene": "FAKEGENE"},
        )
    assert resp.status_code == 400


def test_analyze_result_schema(client):
    """Verify all required output fields are present."""
    with patch("sl_agent.api.routes.run_sl_analysis", return_value=_mock_sl_result()):
        resp = client.post("/api/v1/analyze", json={"gene": "BRCA1"})
    data = resp.json()["result"]

    assert "input_context" in data
    assert "sl_partners" in data
    assert "gene_drug_pairs" in data
    assert "cross_validation" in data

    partner = data["sl_partners"][0]
    for field in ["gene", "delta_dependency", "p_value", "fdr", "n_mut", "n_wt", "test_type"]:
        assert field in partner, f"Missing field: {field}"


def test_async_job_lifecycle(client):
    """Submit async job, poll for queued status."""
    with patch("sl_agent.api.routes.run_sl_analysis", return_value=_mock_sl_result()):
        submit = client.post("/api/v1/analyze/async", json={"gene": "BRCA1"})
    assert submit.status_code == 202
    job_id = submit.json()["job_id"]
    assert job_id

    poll = client.get(f"/api/v1/result/{job_id}")
    assert poll.status_code == 200
    assert poll.json()["status"] in ("queued", "running", "done")
