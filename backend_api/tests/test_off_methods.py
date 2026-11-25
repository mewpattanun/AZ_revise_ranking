"""
Unit tests for the "off cutoff" methods (gene-relative off), synthetic data.

Below 1 TPM a gene is always off. The other methods can also call a cell line
off at 1 TPM or more when it is low for THAT gene (e.g. VIM ~2 TPM in epithelial
lines vs ~1,000 TPM in mesenchymal lines). Default "fixed" = below 1 TPM only.
"""

import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from app.models import RankRequest
from app.services import data_service
from app.services.data_service import classify_expression, off_by_method
from app.services.ranking_service import run_ranking

METHODS = ["fixed", "percent_of_reference", "mixture", "distance_from_on", "strict_mixture", "custom"]


def _ids(prefix, n):
    return [f"{prefix}{i:03d}" for i in range(n)]


def vim_like():
    """Clear on/off gene: 10 lines below 1 TPM, 60 'low' lines ~2 TPM, 60 'on' lines ~1,000 TPM."""
    rng = np.random.default_rng(0)
    off = pd.Series(rng.uniform(0.0, 0.8, 10), index=_ids("Z", 10))
    low = pd.Series(2 ** rng.normal(np.log2(3), 0.3, 60) - 1, index=_ids("L", 60))
    on = pd.Series(2 ** rng.normal(np.log2(1000), 0.5, 60) - 1, index=_ids("H", 60))
    return pd.concat([off, low, on])


def ptens_like():
    """No on/off split: every line around 20 TPM, plus 3 lines below 1 TPM."""
    rng = np.random.default_rng(1)
    on = pd.Series(2 ** rng.normal(np.log2(21), 0.4, 120) - 1, index=_ids("H", 120))
    off = pd.Series([0.1, 0.5, 0.9], index=_ids("Z", 3))
    return pd.concat([on, off])


# ── off_by_method ────────────────────────────────────────────

@pytest.mark.parametrize("method", METHODS)
def test_below_1_tpm_is_always_off(method):
    level = vim_like()
    off, _ = off_by_method(level, method, cutoff_tpm=10)
    assert off[level < 1].all()


def test_fixed_is_below_1_tpm_only():
    level = vim_like()
    off, reason = off_by_method(level, "fixed")
    assert off.equals(level < 1)
    assert reason is None


@pytest.mark.parametrize("method", ["mixture", "distance_from_on", "strict_mixture", "percent_of_reference"])
def test_clear_split_gene_low_group_is_off(method):
    level = vim_like()
    off, reason = off_by_method(level, method)
    low_group = level.index.str.startswith("L")
    on_group = level.index.str.startswith("H")
    assert off[low_group].mean() > 0.9          # the ~2 TPM lines are off for this gene
    assert not off[on_group].any()              # the ~1,000 TPM lines are not
    assert isinstance(reason, str) and reason


@pytest.mark.parametrize("method", ["mixture", "distance_from_on", "strict_mixture"])
def test_gene_without_a_split_falls_back_to_1_tpm(method):
    level = ptens_like()
    off, reason = off_by_method(level, method)
    assert off.equals(level < 1)
    assert reason is None


def test_percent_of_reference_threshold():
    # 10 lines at 100 TPM and 10 at 1..10 TPM: 90th percentile of lines >= 1 TPM = 100,
    # so off below 10 TPM.
    level = pd.Series([100.0] * 10 + list(range(1, 11)) + [0.5], index=_ids("X", 21))
    off, reason = off_by_method(level, "percent_of_reference")
    assert off.equals(level < 10)
    assert "below 10.0 TPM" in reason


def test_custom_cutoff():
    level = pd.Series([0.5, 3.0, 9.9, 10.0, 50.0], index=list("abcde"))
    off, reason = off_by_method(level, "custom", cutoff_tpm=10)
    assert off.tolist() == [True, True, True, False, False]
    assert reason == "below the chosen cutoff of 10 TPM"
    # a cutoff below 1 TPM still keeps "below 1 TPM is off"
    off, _ = off_by_method(level, "custom", cutoff_tpm=0.2)
    assert off.tolist() == [True, False, False, False, False]


def test_custom_needs_a_cutoff_and_unknown_method_fails():
    with pytest.raises(ValueError):
        off_by_method(vim_like(), "custom")
    with pytest.raises(ValueError):
        off_by_method(vim_like(), "lowest")


def test_too_few_lines_falls_back_to_1_tpm():
    level = pd.Series([0.5, 2.0, 2.5, 900.0, 1100.0], index=list("abcde"))
    off, reason = off_by_method(level, "strict_mixture")
    assert off.equals(level < 1) and reason is None


def test_group_fit_is_cached():
    level = vim_like()
    off_by_method(level, "strict_mixture")
    size = len(data_service._GROUP_FIT_CACHE)
    off_by_method(level, "mixture")              # same gene data: no new fit
    assert len(data_service._GROUP_FIT_CACHE) == size


# ── Categories and ranking ───────────────────────────────────

def _signal(level):
    return pd.DataFrame({"tpm_depmap": level, "tpm_hpa": level})


def test_classify_marks_low_relative_only_at_1_tpm_or_more():
    level = vim_like()
    cats = classify_expression(_signal(level), "strict_mixture")
    rel = cats["category"] == "low_relative"
    assert rel.sum() > 50
    assert (level[rel[rel].index] >= 1).all()
    assert cats.loc[rel, "off_reason"].str.contains("clear on/off split").all()
    assert cats.loc[~rel, "off_reason"].isna().all()
    # default: nothing is low_relative, below 1 TPM keeps its old categories
    default = classify_expression(_signal(level))
    assert not (default["category"] == "low_relative").any()
    assert set(default.loc[level < 1, "category"]) <= {"below_detection", "negligible"}


@pytest.fixture
def vim_data(fake_data):
    level = vim_like().round(3).to_dict()
    fake_data["expression"] = {"fact_expression_depmap": dict(level), "fact_expression_hpa": dict(level)}
    return fake_data


LOW = {"hugo": "G1", "ensembl_id": "ENSG_G1", "direction": "low", "weight": 1.0}


def test_low_search_puts_gene_relative_off_lines_in_the_off_group(vim_data):
    results = run_ranking([LOW], top_n=200, off_method="strict_mixture")
    cats = [r["expression"]["G1"]["category"] for r in results]
    first_on = next(i for i, c in enumerate(cats) if c not in data_service.OFF_CATEGORIES)
    assert all(c not in data_service.OFF_CATEGORIES for c in cats[first_on:])   # off group comes first
    assert "low_relative" in cats[:first_on]
    rel = next(r for r in results if r["expression"]["G1"]["category"] == "low_relative")
    assert "clear on/off split" in rel["expression"]["G1"]["off_reason"]


def test_below_1_tpm_comes_before_gene_relative_off_which_is_ordered_by_tpm(vim_data):
    # Inside the off group: lines below 1 TPM first (tied), then "Low for GENE"
    # lines from the lowest TPM up, because 2 vs 20 TPM is a real difference.
    results = run_ranking([LOW], top_n=200, off_method="strict_mixture")
    cats = [r["expression"]["G1"]["category"] for r in results]
    n_below = sum(c in ("below_detection", "negligible") for c in cats)
    assert all(c in ("below_detection", "negligible") for c in cats[:n_below])
    rel = [r["expression"]["G1"]["tpm"] for r in results if r["expression"]["G1"]["category"] == "low_relative"]
    assert cats[n_below:n_below + len(rel)] == ["low_relative"] * len(rel)
    assert rel == sorted(rel)


def test_tpm_greater_than_1_puts_gene_relative_off_lines_last(vim_data):
    results = run_ranking([LOW], top_n=200, off_method="strict_mixture", low_expression_mode="measurable")
    cats = [r["expression"]["G1"]["category"] for r in results]
    assert cats[0] not in data_service.OFF_CATEGORIES
    assert cats[-1] in data_service.OFF_CATEGORIES


def test_default_ranking_is_unchanged_by_the_new_option(vim_data):
    assert [r["ach_id"] for r in run_ranking([LOW], top_n=200)] == \
        [r["ach_id"] for r in run_ranking([LOW], top_n=200, off_method="fixed")]


def test_run_ranking_validates_off_method(vim_data):
    with pytest.raises(ValueError):
        run_ranking([LOW], off_method="lowest")
    with pytest.raises(ValueError):
        run_ranking([LOW], off_method="custom")


# ── API request model ────────────────────────────────────────

def test_request_off_method_defaults_and_validation():
    req = RankRequest(genes=[{"hugo": "VIM"}])
    assert req.off_method == "fixed" and req.off_cutoff_tpm is None
    assert RankRequest(genes=[{"hugo": "VIM"}], off_method="custom", off_cutoff_tpm=10).off_cutoff_tpm == 10
    with pytest.raises(ValidationError):
        RankRequest(genes=[{"hugo": "VIM"}], off_method="custom")          # needs a cutoff
    with pytest.raises(ValidationError):
        RankRequest(genes=[{"hugo": "VIM"}], off_method="custom", off_cutoff_tpm=0.5)   # at least 1
    with pytest.raises(ValidationError):
        RankRequest(genes=[{"hugo": "VIM"}], off_method="lowest")
