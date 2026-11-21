"""
Unit tests for the v3 ranking (synthetic data, no Parquet needed).

Covers:
  - GEO never feeds the score (AZ point 1)
  - RNA consensus robust z (points 2-3)
  - weighted geometric mean for several genes (point 4)
  - per-study GEO support
  - deterministic tie order, core_only per selected source, genes without data
"""

import math

import numpy as np
import pandas as pd
import pytest

from app.models import RankRequest
from app.services import data_service, ranking_service
from app.services.ranking_service import run_ranking, weighted_geometric_mean

GENE = {"hugo": "G1", "ensembl_id": "ENSG_G1", "direction": "high", "weight": 1.0}
GENE2 = {"hugo": "G2", "ensembl_id": "ENSG_G2", "direction": "high", "weight": 1.0}


def _order(results):
    return [r["ach_id"] for r in results]


# ── GEO is never scored ──────────────────────────────────────

def test_geo_table_is_never_read_for_scoring(fake_data):
    run_ranking([GENE], sources=["depmap", "hpa", "geo"])
    assert "fact_expression_geo" not in fake_data["calls"]


def test_extreme_geo_values_do_not_change_rank(fake_data):
    baseline = _order(run_ranking([GENE]))
    fake_data["expression"]["fact_expression_geo"] = {"A": 0.0, "B": 0.0, "C": 0.0, "D": 1e9}
    assert _order(run_ranking([GENE], sources=["depmap", "hpa", "geo"])) == baseline
    assert baseline == ["A", "B", "C", "D"]


def test_geo_only_selection_falls_back_to_scoring_sources(fake_data):
    assert _order(run_ranking([GENE], sources=["geo", "protein"])) == \
        _order(run_ranking([GENE], sources=["depmap", "hpa"]))


def test_request_defaults():
    req = RankRequest(genes=[{"hugo": "EGFR"}])
    assert req.sources == ["depmap", "hpa"]
    assert req.scoring_method == "robust_z"
    assert req.multi_gene_method == "geometric"
    assert "geo" not in ranking_service.SCORING_SOURCE_NAMES


# ── core_only ────────────────────────────────────────────────

def test_core_requires_every_selected_scoring_source(fake_data):
    fake_data["expression"]["fact_expression_hpa"].pop("B")
    by_id = {r["ach_id"]: r for r in run_ranking([GENE])}
    assert by_id["A"]["is_core"] is True
    assert by_id["B"]["is_core"] is False
    # With HPA deselected, DepMap alone is "every selected source".
    by_id = {r["ach_id"]: r for r in run_ranking([GENE], sources=["depmap"])}
    assert by_id["B"]["is_core"] is True


def test_core_only_with_one_source_returns_results(fake_data):
    # Regression: core_only compared coverage with ALL sources, so a
    # single-source search with core_only returned nothing.
    results = run_ranking([GENE], sources=["depmap"], core_only=True)
    assert _order(results) == ["A", "B", "C", "D"]


# ── Deterministic ties ───────────────────────────────────────

@pytest.mark.parametrize("direction", ["high", "low"])
def test_single_gene_ties_are_ordered_by_ach_id(fake_data, direction):
    gene = {**GENE, "direction": direction}
    # D, B, C, A tie; E breaks the tie group so the gene has variance.
    tied = {"D": 5.0, "B": 5.0, "E": 50.0, "C": 5.0, "A": 5.0}
    fake_data["expression"] = {"fact_expression_depmap": tied, "fact_expression_hpa": {}}
    first = _order(run_ranking([gene]))
    fake_data["expression"]["fact_expression_depmap"] = dict(reversed(list(tied.items())))
    assert _order(run_ranking([gene])) == first
    assert [a for a in first if a != "E"] == ["A", "B", "C", "D"]


@pytest.mark.parametrize("multi_gene_method", ["geometric", "summed_rrf"])
@pytest.mark.parametrize("direction", ["high", "low"])
def test_multi_gene_ties_are_deterministic(fake_data, direction, multi_gene_method):
    tied = {"D": 5.0, "B": 5.0, "E": 50.0, "C": 5.0, "A": 5.0}
    fake_data["expression"] = {"fact_expression_depmap": tied, "fact_expression_hpa": {}}
    genes = [{**GENE, "direction": direction}, {**GENE2, "direction": direction}]

    def run():
        return [(r["ach_id"], r["score"]) for r in run_ranking(genes, multi_gene_method=multi_gene_method)]

    first = run()
    fake_data["expression"]["fact_expression_depmap"] = dict(reversed(list(tied.items())))
    assert run() == first
    assert [a for a, _ in first if a != "E"] == ["A", "B", "C", "D"]


# ── Genes without usable data ────────────────────────────────

def test_constant_gene_does_not_crash(fake_data):
    flat = {"C": 0.0, "A": 0.0, "B": 0.0}
    fake_data["expression"] = {"fact_expression_depmap": flat, "fact_expression_hpa": flat}
    results = run_ranking([GENE])
    assert _order(results) == ["A", "B", "C"]
    assert len({r["score"] for r in results}) == 1


def test_gene_with_no_data_returns_empty(fake_data):
    fake_data["expression"] = {}
    assert run_ranking([GENE]) == []
    assert run_ranking([GENE, GENE2]) == []


def test_gene_with_no_data_is_skipped_in_multi_gene(fake_data):
    fake_data["by_gene"]["ENSG_G2"] = {}
    results = run_ranking([GENE, GENE2])
    assert _order(results) == ["A", "B", "C", "D"]
    assert all("No data for G2, skipped" in (r["score_flags"] or []) for r in results)


# ── Robust z and RNA consensus ───────────────────────────────

def test_robust_z_uses_median_and_mad():
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 100.0])
    z = data_service._robust_z(x)
    # median 3, MAD = median(|x - 3|) = 1, scale = 1.4826
    assert z.iloc[2] == 0.0
    assert z.iloc[4] == pytest.approx(97 / 1.4826)
    assert z.iloc[0] == pytest.approx(-2 / 1.4826)


def test_robust_z_falls_back_to_sd_when_mad_is_zero():
    x = pd.Series([0.0, 0.0, 0.0, 0.0, 10.0])
    z, mad_zero = data_service.robust_z_with_flag(x)
    assert mad_zero is True
    assert z.iloc[4] == pytest.approx(10.0 / x.std(ddof=0))


def test_robust_z_of_constant_is_zero():
    z = data_service._robust_z(pd.Series([3.0, 3.0, 3.0]))
    assert (z == 0.0).all()


def test_rna_consensus_is_mean_of_source_z(fake_data):
    fake_data["expression"]["fact_expression_hpa"].pop("D")
    sig = data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"])
    a = sig.loc["A"]
    assert a["consensus"] == pytest.approx((a["z_depmap"] + a["z_hpa"]) / 2)
    assert a["coverage"] == 2
    assert a["discordance"] == pytest.approx(abs(a["z_depmap"] - a["z_hpa"]))
    d = sig.loc["D"]
    assert d["coverage"] == 1
    assert d["consensus"] == pytest.approx(d["z_depmap"])
    assert math.isnan(d["discordance"])


def test_rna_consensus_uses_raw_tpm(fake_data):
    assert data_service.RNA_TRANSFORM == "tpm"
    sig = data_service.get_rna_consensus("ENSG_G1", ["depmap"])
    expected = data_service._robust_z(pd.Series({"A": 100.0, "B": 50.0, "C": 10.0, "D": 1.0}))
    pd.testing.assert_series_equal(sig["z_depmap"], expected, check_names=False)


def test_rna_consensus_log2_option_still_works(fake_data):
    sig = data_service.get_rna_consensus("ENSG_G1", ["depmap"], transform="log2")
    expected = data_service._robust_z(np.log2(pd.Series({"A": 100.0, "B": 50.0, "C": 10.0, "D": 1.0}) + 1))
    pd.testing.assert_series_equal(sig["z_depmap"], expected, check_names=False)



# ── RNA method: robust z, standard z or log2(TPM + 1) ────────

def test_standard_z_uses_mean_and_sd_of_the_expressed_cell_lines():
    s = pd.Series([0.0, 0.3, 10.0, 20.0, 30.0])
    z = data_service._standard_z(data_service._prepare(s, "tpm"), data_service._expressed_reference(s, "tpm"))
    sd = np.array([10.0, 20.0, 30.0]).std()          # ddof = 0, the two cell lines below 1 TPM not counted
    assert z.iloc[3] == pytest.approx(0.0)           # the mean of the expressed lines
    assert z.iloc[4] == pytest.approx(10.0 / sd)
    assert z.iloc[0] == pytest.approx(-20.0 / sd)    # off lines are scored from their real value, below the rest


def test_standard_z_is_zero_when_sd_is_zero_or_nothing_is_expressed():
    assert (data_service._standard_z(pd.Series([5.0, 5.0, 5.0])) == 0).all()
    s = pd.Series([0.0, 0.2, 0.5])
    assert (data_service._standard_z(s, data_service._expressed_reference(s, "tpm")) == 0).all()


def test_an_outlier_inflates_the_sd_but_not_the_mad():
    s = pd.Series([5.0, 6.0, 7.0, 8.0, 500.0])
    ref = data_service._expressed_reference(s, "tpm")
    robust = data_service._robust_z(s, ref).iloc[-1]       # median 7, MAD 1 -> 493 / 1.4826
    standard = data_service._standard_z(s, ref).iloc[-1]   # mean and SD are both dragged up by the 500
    assert robust == pytest.approx(493 / 1.4826)
    assert robust > 50 * standard


def test_rna_consensus_standard_z(fake_data):
    sig = data_service.get_rna_consensus("ENSG_G1", ["depmap"], method="standard_z")
    tpm = pd.Series({"A": 100.0, "B": 50.0, "C": 10.0, "D": 1.0})
    expected = data_service._standard_z(tpm, data_service._expressed_reference(tpm, "tpm"))
    pd.testing.assert_series_equal(sig["z_depmap"], expected, check_names=False)
    assert sig.attrs["rna_method"] == "standard_z"
    assert sig.attrs["mad_zero_sources"] == []


def test_rna_consensus_log2_tpm_does_not_depend_on_the_other_cell_lines(fake_data):
    sources = ["depmap", "hpa"]
    log2 = data_service.get_rna_consensus("ENSG_G1", sources, method="log2_tpm")
    robust = data_service.get_rna_consensus("ENSG_G1", sources)
    assert log2.attrs["rna_method"] == "log2_tpm"
    assert log2.at["B", "z_depmap"] == pytest.approx(np.log2(51.0))
    assert log2.at["B", "z_hpa"] == pytest.approx(np.log2(61.0))
    assert log2.at["B", "consensus"] == pytest.approx((np.log2(51.0) + np.log2(61.0)) / 2)

    fake_data["expression"]["fact_expression_depmap"].pop("A")
    fake_data["expression"]["fact_expression_hpa"].pop("A")
    log2_without_a = data_service.get_rna_consensus("ENSG_G1", sources, method="log2_tpm")
    robust_without_a = data_service.get_rna_consensus("ENSG_G1", sources)
    assert log2_without_a.at["B", "consensus"] == pytest.approx(log2.at["B", "consensus"])      # unchanged
    assert robust_without_a.at["B", "consensus"] != pytest.approx(robust.at["B", "consensus"])  # moved


def test_rna_methods_agree_on_the_raw_tpm_columns(fake_data):
    tpm_cols = {}
    for method in data_service.RNA_METHODS:
        sig = data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"], method=method)
        tpm_cols[method] = sig[["tpm_depmap", "tpm_hpa"]]
    pd.testing.assert_frame_equal(tpm_cols["robust_z"], tpm_cols["standard_z"])
    pd.testing.assert_frame_equal(tpm_cols["robust_z"], tpm_cols["log2_tpm"])


def test_log2_magnitude_score_is_only_min_maxed():
    s = pd.Series([0.0, 2.0, 4.0])
    assert ranking_service.magnitude_score(s, "log2_tpm").tolist() == pytest.approx([0.0, 0.5, 1.0])
    # z-scores still go through asinh: asinh(2) / asinh(4) is about 0.689, not 0.5
    assert ranking_service.magnitude_score(s).iloc[1] == pytest.approx(np.arcsinh(2.0) / np.arcsinh(4.0))


@pytest.mark.parametrize("rna_method", data_service.RNA_METHODS)
def test_run_ranking_works_with_every_rna_method(fake_data, rna_method):
    results = run_ranking([GENE], rna_method=rna_method)
    assert _order(results) == ["A", "B", "C", "D"]
    scores = [r["score"] for r in results]
    assert scores == sorted(scores, reverse=True) and scores[0] > scores[-1]


@pytest.mark.parametrize("rna_method", data_service.RNA_METHODS)
def test_results_show_the_rna_value_of_the_chosen_method(fake_data, rna_method):
    results = {r["ach_id"]: r for r in run_ranking([GENE], rna_method=rna_method)}
    sig = data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"], method=rna_method)
    for ach, r in results.items():
        assert r["expression"]["G1"]["rna_value"] == pytest.approx(sig.at[ach, "consensus"], abs=1e-3)
    if rna_method == "log2_tpm":     # B: DepMap 50, HPA 60
        assert results["B"]["expression"]["G1"]["rna_value"] == pytest.approx((np.log2(51) + np.log2(61)) / 2, abs=1e-3)


def test_default_rna_method_is_robust_z(fake_data):
    assert run_ranking([GENE]) == run_ranking([GENE], rna_method="robust_z")
    assert RankRequest(genes=[{"hugo": "EGFR"}]).rna_method == "robust_z"


def test_unknown_rna_method_is_rejected(fake_data):
    from pydantic import ValidationError

    with pytest.raises(ValueError):
        data_service.get_rna_consensus("ENSG_G1", ["depmap"], method="mean")
    with pytest.raises(ValueError):
        run_ranking([GENE], rna_method="mean")
    with pytest.raises(ValidationError):
        RankRequest(genes=[{"hugo": "EGFR"}], rna_method="mean")


# ── Expressed-only reference: DepMap and HPA agree for sparse genes ──

def _sparse_gene_sources():
    """A gene that is off in 90% of cell lines: DepMap has tiny continuous values,
    HPA reports 0.1 steps (so exact zeros)."""
    rng = np.random.default_rng(1)
    base = np.concatenate([rng.uniform(0, 0.05, 900), rng.lognormal(np.log(80), 0.8, 100)])
    base[-1] = 163.19
    depmap = pd.Series(base)
    hpa = pd.Series(np.round(base * 0.77, 1))
    hpa.iloc[-1] = 126.0
    return depmap, hpa


def _rna_z(s):
    """z as the RNA consensus computes it: median/MAD from the cell lines at or above 1 TPM."""
    return data_service._robust_z(data_service._prepare(s, "tpm"), data_service._expressed_reference(s, "tpm"))


def test_expressed_reference_drops_cell_lines_below_1_tpm():
    s = pd.Series({"a": 0.0, "b": 0.4, "c": 0.99, "d": 1.0, "e": 7.0, "f": 40.0})
    ref = data_service._expressed_reference(s, "tpm")
    assert list(ref.index) == ["d", "e", "f"]
    assert list(ref) == [1.0, 7.0, 40.0]


def test_median_and_mad_come_from_the_expressed_cell_lines_only():
    # 10 off cell lines and 5 expressed ones: median 12, MAD = median(|x - 12|) = 7
    s = pd.Series([0.0] * 5 + [0.3] * 5 + [5.0, 8.0, 12.0, 20.0, 35.0])
    z = _rna_z(s)
    assert z.iloc[12] == 0.0                                  # the median of the expressed lines
    assert z.iloc[-1] == pytest.approx((35.0 - 12.0) / (1.4826 * 7.0))
    # off cell lines are still scored, from their real value, below every expressed one
    assert z.iloc[0] == pytest.approx((0.0 - 12.0) / (1.4826 * 7.0))
    assert z.iloc[:10].max() < z.iloc[10:].min()


def test_expressed_reference_makes_depmap_and_hpa_agree_for_a_sparse_gene():
    depmap, hpa = _sparse_gene_sources()
    z_dep, z_hpa = _rna_z(depmap).iloc[-1], _rna_z(hpa).iloc[-1]
    assert 0.8 < z_dep / z_hpa < 1.25
    # Without removing the cell lines below 1 TPM DepMap's tiny non-zero MAD gives a z thousands of times larger
    plain = data_service._robust_z(data_service._prepare(depmap, "tpm")).iloc[-1]
    assert plain > 100 * z_dep


def test_expressed_reference_keeps_the_order_of_cell_lines():
    depmap, _ = _sparse_gene_sources()
    z = _rna_z(depmap)
    assert z.rank(method="first").equals(depmap.rank(method="first"))


def test_expressed_reference_does_not_change_a_gene_expressed_everywhere():
    s = pd.Series({"A": 5.0, "B": 12.0, "C": 30.0, "D": 80.0, "E": 300.0})
    plain = data_service._robust_z(data_service._prepare(s, "tpm"))
    pd.testing.assert_series_equal(plain, _rna_z(s))


def test_rna_scale_falls_back_to_sd_of_the_expressed_cell_lines_when_mad_is_zero():
    # expressed lines 5, 5, 5, 5, 25: median 5, MAD 0 -> SD of those five = 8 (off lines not counted)
    s = pd.Series([0.0, 0.0, 5.0, 5.0, 5.0, 5.0, 25.0])
    z = _rna_z(s)
    assert z.iloc[-1] == pytest.approx(20.0 / 8.0)
    assert data_service._mad_is_zero(data_service._expressed_reference(s, "tpm"))


def test_gene_off_in_every_cell_line_gets_z_zero_everywhere():
    s = pd.Series({"A": 0.0, "B": 0.2, "C": 0.9})
    assert (_rna_z(s) == 0.0).all()


def test_magnitude_score_keeps_order_and_range():
    z = pd.Series({"a": -2.0, "b": 0.0, "c": 1.0, "d": 60.0})
    s = ranking_service.magnitude_score(z)
    assert list(s.sort_values().index) == ["a", "b", "c", "d"]
    assert s.min() == 0.0 and s.max() == 1.0
    # asinh compression: an extreme outlier does not squash the rest to ~0
    assert s["c"] > 0.2


# ── Weighted geometric mean (multi-gene) ─────────────────────

def test_geometric_mean_one_failed_gene_pulls_score_down():
    scores = pd.DataFrame({"G1": [0.9, 0.6], "G2": [0.0, 0.6]}, index=["one_fails", "both_ok"])
    gm = weighted_geometric_mean(scores)
    assert gm["one_fails"] < 0.01
    assert gm["both_ok"] == pytest.approx(0.6, abs=1e-5)
    # An arithmetic mean would have preferred one_fails (0.45 vs 0.6 is close);
    # the geometric mean makes the failed gene decisive.
    assert gm["both_ok"] > gm["one_fails"]


def test_geometric_mean_skips_missing_gene():
    scores = pd.DataFrame({"G1": [0.8], "G2": [np.nan]}, index=["x"])
    assert weighted_geometric_mean(scores)["x"] == pytest.approx(0.8, abs=1e-5)


def test_geometric_mean_weights():
    scores = pd.DataFrame({"G1": [0.9, 0.3], "G2": [0.3, 0.9]}, index=["x", "y"])
    equal = weighted_geometric_mean(scores)
    assert equal["x"] == pytest.approx(equal["y"])
    weighted = weighted_geometric_mean(scores, {"G1": 3.0, "G2": 1.0})
    assert weighted["x"] == pytest.approx(math.exp((3 * math.log(0.9) + math.log(0.3)) / 4), abs=1e-5)
    assert weighted["x"] > weighted["y"]


def test_gene_weight_changes_multi_gene_order(fake_data):
    # A is best for G1, B is best for G2, symmetric; C is worst for both.
    g1 = {"A": 100.0, "B": 50.0, "C": 1.0}
    g2 = {"A": 50.0, "B": 100.0, "C": 1.0}
    fake_data["by_gene"]["ENSG_G1"] = {"fact_expression_depmap": g1, "fact_expression_hpa": g1}
    fake_data["by_gene"]["ENSG_G2"] = {"fact_expression_depmap": g2, "fact_expression_hpa": g2}
    equal = run_ranking([GENE, GENE2])
    assert equal[0]["score"] == equal[1]["score"]       # symmetric tie, broken by ach_id
    assert _order(run_ranking([{**GENE, "weight": 5.0}, GENE2]))[:2] == ["A", "B"]
    assert _order(run_ranking([GENE, {**GENE2, "weight": 5.0}]))[:2] == ["B", "A"]


def test_low_gene_is_inverted_in_geometric_mean(fake_data):
    # Same data for both genes: HIGH G1 + LOW G2 asks for opposite things,
    # so no cell line can be good at both and the extremes are pulled down.
    results = run_ranking([GENE, {**GENE2, "direction": "low"}])
    by_id = {r["ach_id"]: r for r in results}
    assert by_id["A"]["per_gene_rank"] == {"G1": 1, "G2": 4}
    assert by_id["D"]["per_gene_rank"] == {"G1": 4, "G2": 1}
    assert by_id["A"]["score"] < 0.01 and by_id["D"]["score"] < 0.01


# ── Per-study GEO support ────────────────────────────────────

@pytest.mark.parametrize("z, direction, expected", [
    (1.5, "high", "supportive"),
    (1.0, "high", "supportive"),
    (0.5, "high", "neutral"),
    (-1.0, "high", "contradictory"),
    (-1.5, "low", "supportive"),
    (1.5, "low", "contradictory"),
])
def test_geo_verdict(z, direction, expected):
    assert data_service.geo_verdict(z, direction) == expected


def test_classify_geo_counts_studies():
    study_z = pd.DataFrame({
        "ach_id": ["A", "A", "A", "B", "B", "C", "C"],
        "gse_id": ["S1", "S2", "S3", "S1", "S2", "S1", "S2"],
        "z":      [2.0, 1.5, 0.0, -2.0, -1.2, 2.0, -2.0],
        "n_samples": [10] * 7,
    })
    out = data_service.classify_geo(study_z, "high")
    assert out.loc["A", "geo_status"] == "supportive"
    assert (out.loc["A", "geo_supportive"], out.loc["A", "geo_total"]) == (2, 3)
    assert out.loc["B", "geo_status"] == "contradictory"
    assert out.loc["C", "geo_status"] == "mixed"
    assert data_service.classify_geo(study_z, "low").loc["B", "geo_status"] == "supportive"


def _study(gse, n, high_ach=None, prefix=""):
    """n samples in one study; every sample is its own cell line except one high one."""
    rows = [(f"{prefix}{gse}_cl{i}", f"{gse}_gsm{i}", gse, 10.0 + i % 3) for i in range(n - 1)]
    rows.append((high_ach or f"{gse}_hi", f"{gse}_gsm_hi", gse, 1000.0))
    return rows


def test_geo_study_z_is_within_study_and_skips_small_studies(fake_data):
    fake_data["geo_samples"] = (
        _study("GSE_BIG", data_service.GEO_MIN_STUDY_SAMPLES, high_ach="A")
        + _study("GSE_SMALL", data_service.GEO_MIN_STUDY_SAMPLES - 1, high_ach="B")
    )
    sz = data_service.get_geo_study_z("ENSG_G1")
    assert set(sz["gse_id"]) == {"GSE_BIG"}
    a = sz[sz["ach_id"] == "A"].iloc[0]
    assert a["z"] > data_service.GEO_SUPPORT_Z
    assert a["n_samples"] == data_service.GEO_MIN_STUDY_SAMPLES


def test_geo_support_is_attached_without_changing_order(fake_data):
    fake_data["geo_samples"] = _study("GSE_BIG", 12, high_ach="A")
    without = run_ranking([GENE])
    with_geo = run_ranking([GENE], include_geo=True)
    assert _order(with_geo) == _order(without)
    assert [r["score"] for r in with_geo] == [r["score"] for r in without]
    by_id = {r["ach_id"]: r["geo"]["G1"] for r in with_geo}
    assert by_id["A"]["status"] == "supportive"
    assert by_id["B"]["status"] == "unavailable"
    assert all("geo" not in r or r["geo"] is None for r in without)
