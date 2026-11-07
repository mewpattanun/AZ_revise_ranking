"""
API tests against the real Parquet data (skipped when CLF_DATA_DIR has none).
"""

import pytest

from conftest import real_data

pytestmark = real_data


@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient
    from app.main import app

    return TestClient(app)


def _rank(client, **body):
    res = client.post("/api/rank", json=body)
    assert res.status_code == 200, res.text
    return res.json()


def test_rank_egfr_drops_geo_with_warning(client):
    body = _rank(client, genes=[{"hugo": "EGFR"}], sources=["depmap", "hpa", "geo", "protein"], top_n=10)
    assert body["query"]["sources"] == ["depmap", "hpa"]
    assert body["query"]["scoring_version"]
    assert body["query"]["scoring_method"] == "robust_z"
    assert any("GEO" in w for w in body["warnings"])
    assert len(body["results"]) == 10
    scores = [r["score"] for r in body["results"]]
    assert scores == sorted(scores, reverse=True)


def test_no_warning_without_geo(client):
    assert _rank(client, genes=[{"hugo": "EGFR"}], top_n=1)["warnings"] == []


def test_core_only_with_one_source_returns_results(client):
    both = _rank(client, genes=[{"hugo": "EGFR"}], core_only=True, top_n=50)
    one = _rank(client, genes=[{"hugo": "EGFR"}], sources=["depmap"], core_only=True, top_n=50)
    assert len(both["results"]) == 50
    assert len(one["results"]) == 50
    assert all(r["is_core"] for r in one["results"])


@pytest.mark.parametrize("genes", [
    [{"hugo": "EGFR", "direction": "low"}],
    [{"hugo": "IGHD1-14"}],
    [{"hugo": "ERBB2"}, {"hugo": "EGFR", "direction": "low"}],
])
def test_repeated_requests_give_same_order(client, genes):
    first = [r["ach_id"] for r in _rank(client, genes=genes, top_n=100)["results"]]
    assert first, "expected results"
    assert [r["ach_id"] for r in _rank(client, genes=genes, top_n=100)["results"]] == first


def test_constant_gene_ties_are_ordered_by_ach_id(client):
    # IGHD1-14 is constant, so every score is exactly tied. (Other genes can
    # only be checked for repeatability: the API rounds scores to 4 decimals,
    # so two slightly different scores may look tied.)
    ids = [r["ach_id"] for r in _rank(client, genes=[{"hugo": "IGHD1-14"}], top_n=100)["results"]]
    assert ids == sorted(ids)


def test_gene_without_variance(client):
    # IGHD1-14: HPA all zero, no DepMap, no protein -> constant gene, all tied.
    alone = _rank(client, genes=[{"hugo": "IGHD1-14"}], top_n=5)
    assert len(alone["results"]) == 5
    paired = _rank(client, genes=[{"hugo": "EGFR"}, {"hugo": "IGHD1-14"}], top_n=5)
    assert len(paired["results"]) == 5


def test_evidence_has_per_study_geo_and_no_geo_score(client):
    top = _rank(client, genes=[{"hugo": "EGFR"}], top_n=1)["results"][0]
    res = client.get(f"/api/celllines/{top['ach_id']}/evidence", params={"genes": "EGFR", "directions": "high"})
    assert res.status_code == 200, res.text
    ev = res.json()
    assert "GEO" not in ev.get("source_scores", {})
    egfr = ev["genes"]["EGFR"]
    assert "rna_consensus" in egfr
    if "geo" in egfr:
        assert egfr["geo"]["status"] in {"supportive", "contradictory", "mixed", "unavailable"}
        for s in egfr["geo"]["studies"]:
            assert s["verdict"] in {"supportive", "contradictory", "neutral"}
