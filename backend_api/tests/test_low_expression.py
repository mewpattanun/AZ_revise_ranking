"""
Unit tests for the low-expression mode (AZ priority 5), synthetic data.

"Classify first, then rank": every cell line gets an expression category from
raw TPM, and for a LOW gene the category decides the order before the score.
Thresholds (config): below detection < 0.1 TPM, off < 1 TPM, moderate > 3 TPM.
"""

import pandas as pd
import pytest
from pydantic import ValidationError

from app import config
from app.models import RankRequest
from app.services import data_service
from app.services.ranking_service import run_ranking

LOW = {"hugo": "G1", "ensembl_id": "ENSG_G1", "direction": "low", "weight": 1.0}
HIGH_G2 = {"hugo": "G2", "ensembl_id": "ENSG_G2", "direction": "high", "weight": 1.0}

# ach_id order is the OPPOSITE of the TPM order inside the off group, so a
# noise-driven order (D, C, B, A) and a tied order (A, B, C, D) look different.
DEPMAP = {
    "A": 0.9, "B": 0.6, "C": 0.06, "D": 0.0,   # off (below 1 TPM in both sources)
    "E": 2.0, "F": 1.2,                         # low measurable (1-3 TPM)
    "G": 50.0, "H": 10.0,                       # moderate (> 3 TPM)
}
HPA = {"A": 0.8, "B": 0.3, "C": 0.0, "D": 0.0, "E": 1.5, "F": 1.0, "G": 40.0, "H": 8.0}


def _order(results):
    return [r["ach_id"] for r in results]


@pytest.fixture
def low_data(fake_data):
    fake_data["expression"] = {"fact_expression_depmap": dict(DEPMAP), "fact_expression_hpa": dict(HPA)}
    return fake_data


# ── Categories ───────────────────────────────────────────────

def test_thresholds_are_the_agreed_values():
    assert config.EXPR_OFF_TPM == 1.0
    assert config.EXPR_MODERATE_TPM == 3.0
    assert config.EXPR_DETECTION_TPM == {"depmap": 0.1, "hpa": 0.1}


def test_classify_expression_bands(low_data):
    cats = data_service.classify_expression(data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"]))
    assert cats["category"].to_dict() == {
        "A": "negligible", "B": "negligible", "C": "below_detection", "D": "below_detection",
        "E": "low_measurable", "F": "low_measurable", "G": "moderate", "H": "moderate",
    }
    assert cats.loc["A", "tpm"] == pytest.approx(0.85)
    assert cats.loc["A", "tpm_depmap"] == 0.9 and cats.loc["A", "tpm_hpa"] == 0.8


def test_category_uses_the_highest_source(fake_data):
    # Off only when EVERY source agrees. X: DepMap says off, HPA says expressed
    # -> low measurable, flagged. Y: DepMap 1.17 vs HPA 11.5 (BT549 CDH1) ->
    # moderate; both say expressed, so it is not an off/on disagreement.
    fake_data["expression"] = {
        "fact_expression_depmap": {"X": 0.5, "Y": 1.17, "Z": 20.0},
        "fact_expression_hpa": {"X": 2.0, "Y": 11.5, "Z": 30.0},
    }
    cats = data_service.classify_expression(data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"]))
    assert cats.loc["X", "category"] == "low_measurable"
    assert cats.loc["Y", "category"] == "moderate"
    assert cats["sources_disagree"].astype(bool).to_dict() == {"X": True, "Y": False, "Z": False}


def test_one_source_only_is_classified_on_that_source(low_data):
    low_data["expression"]["fact_expression_hpa"].pop("E")
    cats = data_service.classify_expression(data_service.get_rna_consensus("ENSG_G1", ["depmap", "hpa"]))
    assert cats.loc["E", "category"] == "low_measurable"
    assert cats.loc["E", "tpm"] == 2.0
    assert pd.isna(cats.loc["E", "tpm_hpa"])


# ── Single LOW gene: the three modes ─────────────────────────

def test_off_mode_puts_off_lines_first_and_ties_them(low_data):
    # Off lines tie on RNA, so they follow the ach_id tie rule, not TPM noise.
    assert _order(run_ranking([LOW], low_expression_mode="off")) == ["A", "B", "C", "D", "F", "E", "H", "G"]


def test_off_mode_is_the_default(low_data):
    assert _order(run_ranking([LOW])) == _order(run_ranking([LOW], low_expression_mode="off"))
    assert RankRequest(genes=[{"hugo": "EGFR"}]).low_expression_mode == "off"


def test_noise_below_1_tpm_does_not_change_the_order(low_data):
    first = _order(run_ranking([LOW]))
    # Shuffle the TPM values inside the off group: the order must not change.
    low_data["expression"]["fact_expression_depmap"].update({"A": 0.0, "B": 0.06, "C": 0.9, "D": 0.6})
    low_data["expression"]["fact_expression_hpa"].update({"A": 0.0, "B": 0.0, "C": 0.8, "D": 0.3})
    assert _order(run_ranking([LOW])) == first


def test_measurable_mode_ranks_1_to_3_tpm_first_and_off_last(low_data):
    assert _order(run_ranking([LOW], low_expression_mode="measurable")) == ["F", "E", "H", "G", "A", "B", "C", "D"]


def test_measurable_mode_puts_source_disagreement_after_agreement(low_data):
    # X: DepMap 1.5 but HPA 0 -> low measurable, but the sources disagree, so it
    # comes after E and F even though its consensus score is the lowest.
    low_data["expression"]["fact_expression_depmap"]["X"] = 1.5
    low_data["expression"]["fact_expression_hpa"]["X"] = 0.0
    order = _order(run_ranking([LOW], low_expression_mode="measurable"))
    assert order[:3] == ["F", "E", "X"]
    # In "off" mode X is simply low measurable (tier after the off lines).
    # (only the tier is checked; the order of near-equal 1-2 TPM lines inside it is noise
    # in this 9-line fixture, where the median itself falls on one of the lines)
    assert set(_order(run_ranking([LOW]))[4:7]) == {"X", "F", "E"}


def test_score_mode_keeps_the_old_order(low_data):
    # Legacy: plain ascending score, so noise below 1 TPM decides the order.
    assert _order(run_ranking([LOW], low_expression_mode="score")) == ["D", "C", "B", "A", "F", "E", "H", "G"]


def test_protein_still_orders_lines_inside_the_off_group(low_data):
    low_data["protein"] = {"A": 9.0, "B": 5.0, "C": 3.0, "D": 1.0, "E": 2.0, "F": 2.0, "G": 8.0, "H": 6.0}
    assert _order(run_ranking([LOW]))[:4] == ["D", "C", "B", "A"]


def test_cell_line_without_rna_comes_last(low_data):
    low_data["protein"] = {"P": 0.0, "A": 1.0, "G": 5.0}   # P has protein but no RNA
    results = run_ranking([LOW], top_n=50)
    assert _order(results)[-1] == "P"
    assert results[-1]["expression"]["G1"]["category"] == "missing"


def test_post_filter_resort_keeps_the_tiers(low_data):
    # core_only triggers the post-filter re-sort; the tiers must survive it.
    assert _order(run_ranking([LOW], core_only=True)) == ["A", "B", "C", "D", "F", "E", "H", "G"]


def test_q6_boost_keeps_the_tiers(low_data, monkeypatch):
    # Q6 inflates low scores of non-matching lines; it must not lift a
    # measurable line above an off line.
    from app.services import q6_lineage_scorer

    monkeypatch.setattr(data_service, "get_full_dim", lambda: pd.DataFrame())
    monkeypatch.setattr(
        q6_lineage_scorer, "compute_q6_score",
        lambda *a, **k: pd.DataFrame({
            "ach_id": list(DEPMAP),
            "q6_score": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            "match_level": ["none"] * 4 + ["same_lineage"] * 4,
        }),
    )
    order = _order(run_ranking([LOW], target_lineage="Breast", q6_boost_enabled=True))
    assert set(order[:4]) == {"A", "B", "C", "D"}


# ── HIGH genes are not changed ───────────────────────────────

@pytest.mark.parametrize("mode", ["off", "measurable", "score"])
def test_high_gene_order_does_not_depend_on_the_mode(low_data, mode):
    high = {**LOW, "direction": "high"}
    assert _order(run_ranking([high], low_expression_mode=mode)) == ["G", "H", "E", "F", "A", "B", "C", "D"]


# ── Several genes (geometric mean) ───────────────────────────

def _multi(low_data, mode):
    flat = {k: 5.0 for k in DEPMAP}    # G2 is the same for every line: only G1 decides
    low_data["by_gene"]["ENSG_G2"] = {"fact_expression_depmap": flat, "fact_expression_hpa": flat}
    return _order(run_ranking([LOW, HIGH_G2], low_expression_mode=mode))


def test_multi_gene_off_mode_gives_off_lines_full_desirability(low_data):
    assert _multi(low_data, "off")[:4] == ["A", "B", "C", "D"]


def test_multi_gene_measurable_mode_treats_off_as_a_failed_gene(low_data):
    assert _multi(low_data, "measurable")[-4:] == ["A", "B", "C", "D"]


# ── What the API returns ─────────────────────────────────────

def test_results_show_tpm_and_category(low_data):
    by_id = {r["ach_id"]: r for r in run_ranking([LOW])}
    info = dict(by_id["C"]["expression"]["G1"])
    assert isinstance(info.pop("rna_value"), float)      # the RNA z of the chosen method, checked in test_ranking.py
    assert info == {
        "category": "below_detection",
        "tpm": 0.03,
        "tpm_by_source": {"depmap": 0.06, "hpa": 0.0},
        "sources_disagree": False,
    }
    assert by_id["G"]["expression"]["G1"]["category"] == "moderate"
    assert "_tier" not in by_id["G"]


def test_disagreement_is_flagged(fake_data):
    fake_data["expression"] = {
        "fact_expression_depmap": {"X": 0.5, "Y": 20.0, "Z": 0.1},
        "fact_expression_hpa": {"X": 2.5, "Y": 30.0, "Z": 0.2},
    }
    by_id = {r["ach_id"]: r for r in run_ranking([LOW])}
    def disagree(row):   # DepMap is almost all below 1 TPM in this fixture, so a MAD=0 note may also appear
        return [f for f in (row["score_flags"] or []) if f.startswith("Sources disagree")]
    assert disagree(by_id["X"]) == ["Sources disagree on expressed vs off (DepMap 0.5, HPA 2.5 TPM)"]
    assert disagree(by_id["Y"]) == []


def test_unknown_mode_is_rejected(low_data):
    with pytest.raises(ValidationError):
        RankRequest(genes=[{"hugo": "EGFR"}], low_expression_mode="lowest")
    with pytest.raises(ValueError):
        run_ranking([LOW], low_expression_mode="lowest")
