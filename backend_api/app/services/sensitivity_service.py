"""
Weight-sensitivity analysis (AZ priority 7).

Scientists often ask: "If I change the RNA:Protein weight, does my top
candidate change?" Priority 7 answers this by running the ranking at
several weight settings and reporting (a) how much the top-K overlap
changes, (b) a rank-correlation across settings, and (c) which
candidates are "weight-sensitive" and therefore need a second look.

This service does NOT change the ranking algorithm or the default
weight — it only makes the dependence visible.
"""

from typing import Any
from scipy.stats import spearmanr
from app.services.ranking_service import run_ranking


# 5 weight settings tested in the sensitivity sweep 
# Default (70:30) is marked so overlap/correlation can be reported
# relative to it.
WEIGHT_SETTINGS: list[tuple[str, float, float]] = [
    ("rna_only", 1.0, 0.0),
    ("90_10",    0.9, 0.1),
    ("70_30",    0.7, 0.3),   # default
    ("50_50",    0.5, 0.5),
    ("30_70",    0.3, 0.7),
]

DEFAULT_SETTING = "70_30"

# Thresholds for flagging a "weight-sensitive" candidate.
# A cell is flagged only when BOTH conditions hold, so candidates far
# down the list (where position changes don't matter to the scientist)
# are not reported as noise.
SENSITIVE_MIN_RANK_CHANGE = 10   # rank must move by at least this many positions
SENSITIVE_MAX_DEFAULT_RANK = 50  # only cells inside this default-rank window count


# Core sensitivity computation 

def compute_sensitivity(
    resolved_genes: list[dict],
    base_request_kwargs: dict[str, Any],
) -> dict[str, Any]:
    """
    Run `run_ranking` at every weight setting in WEIGHT_SETTINGS using
    the same filters and options as the user's original request, and
    return a structured report.

    Parameters
    ----------
    resolved_genes : list[dict]
        The gene list already resolved to Ensembl IDs (as the /api/rank
        router builds it). Passed straight to run_ranking.
    base_request_kwargs : dict
        Everything else from the RankRequest that should be held
        constant across the sweep (filters, mode, sources, etc.).
        `w_rna` and `w_protein` are overwritten per setting, so the
        caller can safely include or omit them.

    Returns
    -------
    dict with:
      - settings_tested           list of labels
      - default_setting           label of the default (70_30)
      - rankings_per_setting      {label: [{ach_id, rank}, ...]}
      - top_10_overlap            {label: float in [0, 1]}
      - rank_correlation          {label: Spearman rho vs default}
      - sensitive_candidates      list[dict], see _find_sensitive()
    """
    sweep_kwargs = {k: v for k, v in base_request_kwargs.items()
                    if k not in ("w_rna", "w_protein")}

    rankings_per_setting: dict[str, list[dict]] = {}
    for label, w_rna, w_prot in WEIGHT_SETTINGS:
        results = run_ranking(
            genes=resolved_genes,
            w_rna=w_rna,
            w_protein=w_prot,
            **sweep_kwargs,
        )
        # Keep only what sensitivity needs; downstream code should never
        # need the full result rows here.
        rankings_per_setting[label] = [
            {"ach_id": r["ach_id"], "rank": r["rank"]} for r in results
        ]

    return {
        "settings_tested": [label for label, _, _ in WEIGHT_SETTINGS],
        "default_setting": DEFAULT_SETTING,
        "rankings_per_setting": rankings_per_setting,
        "top_10_overlap": _top_k_overlap(rankings_per_setting, k=10),
        "rank_correlation": _rank_correlation(rankings_per_setting),
        "sensitive_candidates": _find_sensitive(
            rankings_per_setting,
            min_change=SENSITIVE_MIN_RANK_CHANGE,
            top_n=SENSITIVE_MAX_DEFAULT_RANK,
        ),
    }


# Metric helpers 

def _top_k_overlap(rankings: dict[str, list[dict]], k: int = 10) -> dict[str, float]:
    """
    Fraction of each setting's top-K that is also in the default's top-K.
    The default compared against itself is always 1.0 — kept so the UI
    can render a complete row without special casing.
    """
    default_top = {r["ach_id"] for r in rankings[DEFAULT_SETTING][:k]}
    if not default_top:
        return {label: 0.0 for label in rankings}
    return {
        label: len(default_top & {r["ach_id"] for r in rows[:k]}) / k
        for label, rows in rankings.items()
    }


def _rank_correlation(rankings: dict[str, list[dict]]) -> dict[str, float]:
    """
    Spearman rank correlation between each setting and the default,
    computed on the intersection of cells they both rank. Default is
    reported as 1.0 for the same UI-completeness reason as above.
    Returns 0.0 when correlation is undefined (e.g. no overlap).
    """
    default_ranks = {r["ach_id"]: r["rank"] for r in rankings[DEFAULT_SETTING]}
    out: dict[str, float] = {}
    for label, rows in rankings.items():
        if label == DEFAULT_SETTING:
            out[label] = 1.0
            continue
        setting_ranks = {r["ach_id"]: r["rank"] for r in rows}
        common = sorted(set(default_ranks) & set(setting_ranks))
        if len(common) < 2:
            out[label] = 0.0
            continue
        rho, _ = spearmanr(
            [default_ranks[a] for a in common],
            [setting_ranks[a] for a in common],
        )
        out[label] = float(rho) if rho == rho else 0.0   # guard against NaN
    return out


def _find_sensitive(
    rankings: dict[str, list[dict]],
    min_change: int,
    top_n: int,
) -> list[dict]:
    """
    Flag cells whose rank moves across settings. A cell is reported only if

      1. It sits inside the top-`top_n` at the DEFAULT weight
         (candidates the scientist will actually read).
      2. Its rank changes by at least `min_change` positions in at least
         one other setting.

    A cell missing from a setting is treated as rank = len(ranking) + 1
    for that setting, so disappearances count as big moves.
    """
    default_rows = rankings[DEFAULT_SETTING]
    default_ranks = {r["ach_id"]: r["rank"] for r in default_rows}
    worst_rank = len(default_rows) + 1

    # All cells in the default top_n
    candidates = [r["ach_id"] for r in default_rows[:top_n]]

    sensitive: list[dict] = []
    for ach_id in candidates:
        default_rank = default_ranks[ach_id]
        ranks_per_setting: dict[str, int] = {}
        for label, rows in rankings.items():
            row = next((r for r in rows if r["ach_id"] == ach_id), None)
            ranks_per_setting[label] = row["rank"] if row else worst_rank

        rank_values = list(ranks_per_setting.values())
        max_change = max(rank_values) - min(rank_values)
        if max_change >= min_change:
            sensitive.append({
                "ach_id": ach_id,
                "default_rank": default_rank,
                "max_rank_change": max_change,
                "ranks_per_setting": ranks_per_setting,
            })

    # Report the most sensitive candidates first
    return sorted(sensitive, key=lambda c: -c["max_rank_change"])
