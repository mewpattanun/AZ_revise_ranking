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


# ── Low-expression mode (priority 5) ─────────────────────────

OFF = {"below_detection", "negligible"}


def test_low_esr1_off_mode_returns_off_lines_first(client):
    body = _rank(client, genes=[{"hugo": "ESR1", "direction": "low"}], top_n=50)
    assert body["query"]["low_expression"]["mode"] == "off"
    assert body["query"]["low_expression"]["off_tpm"] == 1.0
    for r in body["results"]:
        info = r["expression"]["ESR1"]
        assert info["category"] in OFF
        assert all(v < 1.0 for v in info["tpm_by_source"].values())


def test_low_esr1_measurable_mode_returns_1_to_3_tpm_first(client):
    body = _rank(
        client, genes=[{"hugo": "ESR1", "direction": "low"}], low_expression_mode="measurable", top_n=10
    )
    for r in body["results"]:
        info = r["expression"]["ESR1"]
        assert info["category"] == "low_measurable"
        assert not info["sources_disagree"]
        assert all(1.0 <= v <= 3.0 for v in info["tpm_by_source"].values())


@pytest.mark.parametrize("gene, ach_id, expected", [
    ("ESR1", "ACH-000768", "below_detection"),   # MDA-MB-231, ER-negative
    ("PTEN", "ACH-000090", "negligible"),        # PC3, PTEN deleted
    ("CDKN2A", "ACH-000681", "negligible"),      # A549, CDKN2A deleted
    ("ESR1", "ACH-000019", "moderate"),          # MCF7, ER-positive
    # Known limitation: HCC1954 is ER-negative but has ~3 TPM ESR1 RNA, so a
    # fixed TPM cutoff cannot call it off (needs a gene-relative rule).
    ("ESR1", "ACH-000859", "moderate"),
])
def test_known_cell_lines_get_the_expected_category(gene, ach_id, expected):
    from app.services import data_service

    ensg = data_service.resolve_gene(gene)["ensembl_id"]
    cats = data_service.classify_expression(data_service.get_rna_consensus(ensg))
    assert cats.loc[ach_id, "category"] == expected


def test_protein_100_percent_has_no_rna_only_lines(client):
    # EGFR high with protein 100%: RNA-only lines (e.g. HSC1, A431) used to rank
    # on RNA alone; now only lines with protein data are scored.
    body = _rank(client, genes=[{"hugo": "EGFR"}], w_rna=0.0, w_protein=1.0, top_n=50)
    assert body["results"]
    assert all(r["scenario"] == "RNA+Protein" for r in body["results"])
    names = {r["cell_line_name"] for r in body["results"]}
    assert "HSC1" not in names and "A431" not in names


def test_both_weights_zero_returns_422(client):
    res = client.post("/api/rank", json={"genes": [{"hugo": "EGFR"}], "w_rna": 0, "w_protein": 0})
    assert res.status_code == 422


def test_note_for_erbb2_high_cd3e_low(client):
    # CD3E is a T-cell gene: every top ERBB2 line is off for it, so the API
    # explains why adding "CD3E low" does not change the list.
    body = _rank(client, genes=[{"hugo": "ERBB2"}, {"hugo": "CD3E", "direction": "low"}], top_n=20)
    assert len(body["notes"]) == 1
    assert body["notes"][0].startswith("CD3E (LOW) is below 1 TPM in all 20 cell lines shown")
    assert _rank(client, genes=[{"hugo": "ERBB2"}, {"hugo": "ESR1", "direction": "low"}], top_n=20)["notes"] == []


def test_erbb2_high_esr1_low_pushes_esr1_expressing_lines_below_off_lines(client):
    # HCC1954 (ESR1 ~3 TPM) stayed at #17 when the LOW gene only lowered the
    # geometric mean a little. Now every ESR1-off line comes first.
    body = _rank(client, genes=[{"hugo": "ERBB2"}, {"hugo": "ESR1", "direction": "low"}], top_n=200)
    cats = [r["expression"]["ESR1"]["category"] for r in body["results"]]
    assert all(c in OFF for c in cats)                      # 1,370 lines are off, so the top 200 all are
    assert "HCC1954" not in {r["cell_line_name"] for r in body["results"]}
