"""
Core ranking algorithm.

Pipeline:
  1. Expression score: DepMap + HPA are standardised separately (gene-wise
     robust z of TPM) and averaged into ONE RNA consensus. The score is
     computed from the SIZE of that consensus z (asinh-compressed, then
     min-max to [0, 1]), so a cell line far above the rest scores clearly
     higher than one just above the median. The rank is only a display of the
     sorted score. (Legacy rank-based methods rrf / zscore / percentile are
     still selectable for comparison.)
  2. Protein score: same as RNA, from the size of the protein robust z
     (asinh → [0, 1]); legacy Z-Score + Percentile → RRF for the legacy methods
  3. Combined RNA + Protein (4 scenarios, configurable weights)
  4. Mutation/Fusion hard filter (include/exclude/ignore)
  5. Multi-gene combination: weighted geometric mean of the per-gene scores
     (desirability in 0-1, LOW genes inverted), so one failed gene pulls the
     score down. The old summed RRF is still selectable (multi_gene_method).
  6. Optional post-filters (disease, lineage, core-only)

Low-expression mode (AZ priority 5): every cell line gets an expression
category from raw TPM (below detection / negligible / low measurable /
moderate, see data_service.classify_expression). For LOW genes the category
decides the order first ("classify first, then rank"): below 1 TPM all cell
lines tie on RNA, so noise like 0.03 vs 0.05 TPM no longer orders them.

Note on RRF averaging:
  When fusing the per-source expression ranks, we AVERAGE the RRF contributions
  (divide by the number of rank lists that actually cover each cell line) instead
  of summing them. This prevents cell lines that are missing a data source (e.g.
  no HPA) from being penalised: a cell line with very high RNA in the source
  it has is judged on the quality of those ranks, not on how many sources exist.
  Coverage is surfaced separately (source_count / is_core) so users who want
  multi-source validation can filter for it instead of it being baked into score.
"""

import numpy as np
import pandas as pd
from scipy import stats

from app.config import DEFAULT_LOW_EXPRESSION_MODE, LOW_EXPRESSION_MODES, RRF_K
from app.services import data_service


# ── Statistical helpers ───────────────────────────────────────

def _zscore_rank(series: pd.Series) -> pd.Series:
    """Z-score normalise, then rank descending (rank 1 = highest z-score)."""
    z = stats.zscore(series, nan_policy="omit")
    return pd.Series(z, index=series.index).rank(ascending=False, method="min")


def _percentile_rank(series: pd.Series) -> pd.Series:
    """Percentile rank, then rank descending (rank 1 = highest percentile)."""
    pct = series.rank(pct=True, method="average")
    return pct.rank(ascending=False, method="min")


def _rrf(rank_dict: dict[str, pd.Series], k: int = RRF_K, average: bool = False) -> pd.Series:
    """
    Reciprocal Rank Fusion across multiple rank Series.

    average=False (default): raw sum of 1/(k+rank) — used for multi-gene fusion.
    average=True: divide each cell line's summed score by the number of rank
                  lists that actually covered it. This makes cell lines with
                  fewer data sources comparable to those with more, so a strong
                  cell line missing one source is not pushed down the ranking.
    """
    all_ids: set[str] = set()
    for s in rank_dict.values():
        all_ids.update(s.index)

    scores = pd.Series(0.0, index=list(all_ids))
    counts = pd.Series(0, index=list(all_ids))
    for ranks in rank_dict.values():
        for idx in ranks.index:
            scores[idx] += 1.0 / (k + ranks[idx])
            counts[idx] += 1

    if average:
        # Avoid division by zero; counts is >= 1 wherever scores was touched.
        counts = counts.replace(0, 1)
        scores = scores / counts

    return scores


def _sort_by_score(df: pd.DataFrame, score_col: str, ascending: bool) -> pd.DataFrame:
    """
    Sort by score with ties broken by ach_id, so the same query always gives
    the same order. Without this, tied cell lines follow DuckDB's (unordered)
    row order. Ties are common: "low" searches put many cell lines at the
    floor, and a constant gene gives every cell line the same score.
    """
    return df.sort_values(
        [score_col, "ach_id"], ascending=[ascending, True], kind="mergesort"
    ).reset_index(drop=True)


def _norm01(s: pd.Series) -> pd.Series:
    """Min-max normalise to [0, 1]."""
    if len(s) == 0:
        return s
    mn, mx = s.min(), s.max()
    if mx == mn:
        return pd.Series(0.5, index=s.index)
    return (s - mn) / (mx - mn)


# ── Expression ranking for one gene ──────────────────────────

SCORING_METHODS = ("robust_z", "rrf", "zscore", "percentile")

# RNA sources that feed the score. GEO is evidence only and never scored.
SCORING_SOURCE_NAMES = list(data_service.RNA_SOURCE_TABLES)


def magnitude_score(consensus: pd.Series, rna_method: str = "robust_z") -> pd.Series:
    """
    Score from the SIZE of the consensus (keeps magnitude, unlike a rank).

    asinh(z) is ~z near 0 and ~log(2z) for big z, so an outlier at z = 60 does
    not squash every other cell line to ~0 when we min-max to [0, 1]. It is
    monotonic, so the order is the order of the consensus z (HIGH: biggest
    first, LOW: ascending, same as before).
    log2(TPM + 1) values are already compressed by the log, so for rna_method
    "log2_tpm" they are only min-maxed.
    """
    values = consensus.astype(float)
    if rna_method == "log2_tpm":
        return _norm01(values)
    return _norm01(np.arcsinh(values))


def score_expression(signal: pd.DataFrame, scoring_method: str = "robust_z") -> pd.Series:
    """
    Expression score per cell line from the RNA consensus signal.
    scoring_method:
      "robust_z"  → magnitude of the RNA consensus (default, point 3). The consensus is the
                    robust z, standard z or log2(TPM + 1) the signal was built with
                    (signal.attrs["rna_method"], see data_service.RNA_METHODS)
      "rrf" / "zscore" / "percentile" → legacy rank-based (RRF of z-rank / pct-rank)
    """
    if scoring_method not in SCORING_METHODS:
        raise ValueError(f"unknown scoring_method: {scoring_method!r}")
    if len(signal) < 2:
        return pd.Series(dtype=float)
    consensus = signal["consensus"]

    if scoring_method == "robust_z":
        return magnitude_score(consensus, signal.attrs.get("rna_method", "robust_z"))

    rank_dict: dict[str, pd.Series] = {}
    if scoring_method in ("rrf", "zscore"):
        rank_dict["rna_zscore"] = _zscore_rank(consensus)
    if scoring_method in ("rrf", "percentile"):
        rank_dict["rna_percentile"] = _percentile_rank(consensus)
    return _rrf(rank_dict, average=True)


def rank_expression(
    ensembl_id: str,
    scoring_method: str = "robust_z",
    enabled_sources: list[str] | None = None,
    rna_method: str = "robust_z",
) -> tuple[pd.Series, pd.Series]:
    """
    Compute the expression score for one gene.
    enabled_sources: RNA sources that feed the consensus (e.g. ["depmap", "hpa"]).
                     None means all sources.
    rna_method: how each source is put on a common scale (data_service.RNA_METHODS).
    Returns (scores indexed by ach_id, source_coverage Series).
    """
    # DepMap + HPA share the same underlying data → one RNA consensus per cell line
    # (mean of per-source robust z-scores; discordance is QC only, see data_service)
    signal = data_service.get_rna_consensus(ensembl_id, enabled_sources, method=rna_method)
    if len(signal) < 2:
        return pd.Series(dtype=float), pd.Series(dtype=int)
    return score_expression(signal, scoring_method), signal["coverage"].astype(int)


def _gene_flags(signal: pd.DataFrame, prot_scores: pd.Series) -> dict[str, list[str]]:
    """
    Flags per cell line for one gene: sources (RNA and/or protein) where MAD = 0,
    so the robust z used the SD instead. Only sources that cover the cell line.
    """
    flags: dict[str, list[str]] = {}
    for src in signal.attrs.get("mad_zero_sources", []):
        col = f"z_{src}"
        if col not in signal.columns:
            continue
        label = data_service.RNA_SOURCE_LABELS.get(src, src)
        for ach in signal.index[signal[col].notna()]:
            flags.setdefault(ach, []).append(f"{label} MAD=0 (SD used)")
    if prot_scores.attrs.get("mad_zero"):
        for ach in prot_scores.index:
            flags.setdefault(ach, []).append("Protein MAD=0 (SD used)")
    return flags


# ── Low-expression mode (AZ priority 5) ──────────────────────

# Order of the categories for a LOW gene, per mode (lower tier = shown first).
# "missing" = no RNA value (protein only) always comes last.
_LOW_TIERS = {
    "off": {"below_detection": 0, "negligible": 0, "low_measurable": 1, "moderate": 2},
    "measurable": {"low_measurable": 0, "moderate": 2, "below_detection": 3, "negligible": 3},
}
_MISSING_TIER = 4


def low_expression_tier(category: str | None, mode: str, sources_disagree: bool = False) -> int:
    """
    Sort tier of a cell line for a LOW gene (0 = best). Mode "score": no tiers.
    In "measurable" mode a low measurable line where one source says off
    (e.g. DepMap 1.4, HPA 0) is less certain, so it comes after the lines
    where both sources agree.
    """
    if mode not in _LOW_TIERS:
        return 0
    tier = _LOW_TIERS[mode].get(category, _MISSING_TIER)
    if mode == "measurable" and category == "low_measurable" and sources_disagree:
        tier = 1
    return tier


def apply_low_expression_mode(expr: pd.Series, categories: pd.DataFrame, mode: str) -> pd.Series:
    """
    Set the RNA score of off cell lines (below 1 TPM in every source) for a LOW gene.

      off        → all off cell lines share the best low score, so noise below
                   1 TPM no longer orders them; protein can still separate them
      measurable → all off cell lines share the worst score (being off is not
                   what was asked for; in a multi-gene search this pulls the
                   geometric mean down like any other failed gene)
      score      → unchanged (old behaviour)

    Works on any score where LOWER = more "low" (the score is then sorted
    ascending, or inverted for the geometric mean).
    """
    if mode not in _LOW_TIERS or expr.empty or categories.empty:
        return expr
    is_off = categories["category"].isin(data_service.OFF_CATEGORIES)
    off = categories.index[is_off].intersection(expr.index)
    if len(off) == 0:
        return expr
    span = float(expr.max() - expr.min()) or 1.0
    out = expr.astype(float).copy()
    # Just outside the range of the other cell lines, so off lines are strictly
    # first (or last) after the min-max scaling in combine_rna_protein.
    if mode == "off":
        out[off] = float(expr.min()) - 1e-6 * span
    else:
        out[off] = float(expr.max()) + 1e-6 * span
    return out


def apply_low_expression_to_protein(
    prot_rrf: pd.Series,
    categories: pd.DataFrame,
    mode: str,
) -> pd.Series:
    """
    P6 LOW-mode fix: extend P5's off-tier treatment to protein.

    apply_low_expression_mode (P5) pushes RNA of off cells to the extreme
    so noise below 1 TPM does not order them. Protein was left untouched,
    which created a bug in mixed-coverage off tiers: cells without protein
    kept combined_score = 0 while cells with low protein scored > 0, so
    RNA-only cells unfairly ranked above cells with measured protein.

    This helper mirrors the push for protein, but only when coverage is
    MIXED (some off cells have protein, some don't). When every off cell
    has protein (Prab's P5 Case A), protein is left alone so it can still
    order the off group — keeping test_protein_still_orders_lines_inside_
    the_off_group valid.
    """
    if mode not in _LOW_TIERS or prot_rrf.empty or categories.empty:
        return prot_rrf
    is_off = categories["category"].isin(data_service.OFF_CATEGORIES)
    off = categories.index[is_off]
    if len(off) == 0:
        return prot_rrf

    # Case A (preserve P5): every off cell has protein → leave alone
    off_in_prot = off.intersection(prot_rrf.index)
    if len(off_in_prot) == len(off):
        return prot_rrf

    # Case B (P6 fix): mixed coverage → tie off cells on protein too,
    # mirroring the RNA push direction so sort order stays consistent.
    span = float(prot_rrf.max() - prot_rrf.min()) or 1.0
    out = prot_rrf.astype(float).copy()
    if mode == "off":
        out[off_in_prot] = float(prot_rrf.min()) - 1e-6 * span
    else:
        out[off_in_prot] = float(prot_rrf.max()) + 1e-6 * span
    return out


def _expression_info(categories: pd.DataFrame, ach_id: str, consensus: pd.Series | None = None) -> dict:
    """Category and raw TPM of one cell line for one gene (shown in the results).
    rna_value = the RNA consensus of the chosen rna_method (robust z, standard z or log2(TPM + 1))."""
    if ach_id not in categories.index:
        return {"category": "missing", "tpm": None, "tpm_by_source": {}, "sources_disagree": False, "rna_value": None}
    row = categories.loc[ach_id]
    by_source = {
        c[len("tpm_"):]: round(float(row[c]), 2)
        for c in categories.columns
        if c.startswith("tpm_") and pd.notna(row[c])
    }
    value = consensus.get(ach_id) if consensus is not None else None
    return {
        "category": row["category"],
        "tpm": round(float(row["tpm"]), 2),
        "tpm_by_source": by_source,
        "sources_disagree": bool(row["sources_disagree"]),
        "rna_value": round(float(value), 3) if value is not None and pd.notna(value) else None,
    }


def _disagree_flags(categories: pd.DataFrame) -> dict[str, list[str]]:
    """Flag cell lines where one source says off (< 1 TPM) and another says expressed."""
    flags: dict[str, list[str]] = {}
    if categories.empty:
        return flags
    for ach in categories.index[categories["sources_disagree"].astype(bool)]:
        info = _expression_info(categories, ach)
        values = ", ".join(
            f"{data_service.RNA_SOURCE_LABELS.get(s, s)} {v:g}" for s, v in info["tpm_by_source"].items()
        )
        flags.setdefault(ach, []).append(f"Sources disagree on expressed vs off ({values} TPM)")
    return flags


# ── Protein ranking for one gene ─────────────────────────────

def rank_protein(ensembl_id: str, scoring_method: str = "robust_z") -> pd.Series:
    """
    Protein score for one gene.
    robust_z (default): same as RNA — robust z of the intensities across cell
    lines → asinh → [0, 1]. scores.attrs["mad_zero"] says whether MAD was 0.
    Legacy methods: Z-Score + Percentile → RRF.
    """
    df = data_service.get_protein_for_gene(ensembl_id)
    if len(df) < 2:
        return pd.Series(dtype=float)

    intensity = df.set_index("ach_id")["protein_intensity"].astype(float)

    if scoring_method == "robust_z":
        z, mad_zero = data_service.robust_z_with_flag(intensity)
        scores = magnitude_score(z)
        scores.attrs["mad_zero"] = mad_zero
        return scores

    rank_dict = {
        "protein_zscore": _zscore_rank(intensity),
        "protein_percentile": _percentile_rank(intensity),
    }
    # Protein has a single source; averaging vs summing gives the same ordering,
    # but we average for consistency with the expression side.
    scores = _rrf(rank_dict, average=True)
    scores.attrs["mad_zero"] = False
    return scores


# ── Combine RNA + Protein (4 scenarios) ──────────────────────

def combine_rna_protein(
    expr_rrf: pd.Series,
    prot_rrf: pd.Series,
    w_rna: float,
    w_prot: float,
    invert: bool = False,
) -> pd.DataFrame:
    """
    Combine expression and protein scores (P6: score and coverage are separate).

    The score is the weighted mean over the data types the cell line HAS, so a
    missing data type is not counted as zero:
      RNA+Protein  : (w_rna * RNA + w_prot * protein) / (w_rna + w_prot)
      RNA only     : RNA score
      Protein only : protein score
    A data type with weight 0 carries no evidence: a cell line whose only data
    type has weight 0 (e.g. RNA only with protein at 100%) is left out, because
    under those weights there is nothing to score it on.
    modality_coverage = share of the total weight the cell line has data for.

    invert=True flips both to "low is good" (1 - score), used for LOW genes in
    multi-gene searches so that a higher score is always better.
    Returns DataFrame with: ach_id, combined_score, scenario, modality_coverage.
    """
    w_total = w_rna + w_prot
    if w_total <= 0:
        raise ValueError("w_rna and w_protein cannot both be 0")
    expr_norm = _norm01(expr_rrf)
    prot_norm = _norm01(prot_rrf)
    if invert:
        expr_norm = 1.0 - expr_norm
        prot_norm = 1.0 - prot_norm
    all_ids = set(expr_norm.index) | set(prot_norm.index)

    rows = []
    for ach in all_ids:
        if ach is None or (isinstance(ach, float) and np.isnan(ach)):
            continue
        # What data the cell line HAS (shown as the scenario) ...
        has_rna = ach in expr_norm.index and not np.isnan(expr_norm.get(ach, np.nan))
        has_prot = ach in prot_norm.index and not np.isnan(prot_norm.get(ach, np.nan))
        if has_rna and has_prot:
            scenario = "RNA+Protein"
        elif has_rna:
            scenario = "RNA only"
        elif has_prot:
            scenario = "Protein only"
        else:
            continue

        # ... and what counts in the score: only data types with weight > 0.
        use_rna = has_rna and w_rna > 0
        use_prot = has_prot and w_prot > 0
        if use_rna and use_prot:
            score = (w_rna * expr_norm[ach] + w_prot * prot_norm[ach]) / w_total
        elif use_rna:
            score = float(expr_norm[ach])
        elif use_prot:
            score = float(prot_norm[ach])
        else:
            continue   # its only data type has weight 0: nothing to score it on
        modality_coverage = ((w_rna if use_rna else 0.0) + (w_prot if use_prot else 0.0)) / w_total

        rows.append({
            "ach_id": ach,
            "combined_score": score,
            "scenario": scenario,
            "modality_coverage": modality_coverage,
        })

    # Explicit columns: a gene with no usable data must give an empty table
    # (the caller returns no results), not a column-less frame that crashes.
    return pd.DataFrame(rows, columns=["ach_id", "combined_score", "scenario", "modality_coverage"])


# ── Multi-gene fusion ────────────────────────────────────────

MULTI_GENE_METHODS = ("geometric", "summed_rrf")
GEOMETRIC_EPS = 1e-6   # numerical stability only (log of 0); not a tuning knob


def weighted_geometric_mean(
    scores: pd.DataFrame,
    weights: dict[str, float] | None = None,
    eps: float = GEOMETRIC_EPS,
) -> pd.Series:
    """
    Weighted geometric mean of per-gene desirability, one row per cell line:

        exp[ sum_g w_g * log(D_g + eps) / sum_g w_g ]

    scores: DataFrame (rows = cell lines, columns = genes) with D_g in (0, 1],
            higher = better for every gene (LOW genes already inverted).
            NaN = the cell line has no score for that gene.
    Missing-gene policy: a gene without a score is SKIPPED for that cell line
    (the weights of the remaining genes are renormalised), so a cell line is
    not penalised for data it does not have. Cell lines with no score at all get NaN.
    """
    w = pd.Series({g: 1.0 for g in scores.columns})
    if weights:
        w.update(pd.Series(weights, dtype=float).reindex(scores.columns).dropna())
    have = scores.notna()
    wsum = have.mul(w, axis=1).sum(axis=1)
    logs = np.log(scores.astype(float).clip(lower=0.0) + eps).fillna(0.0)
    num = logs.mul(w, axis=1).sum(axis=1)
    return np.exp(num / wsum.where(wsum > 0))


# ── Mutation / Fusion hard filters ───────────────────────────

def apply_mutation_filter(
    df: pd.DataFrame, mode: str, ensembl_id: str
) -> pd.DataFrame:
    """Apply hard mutation filter: include / exclude / ignore."""
    if mode == "ignore":
        return df
    mutated = data_service.get_mutated_cell_lines(ensembl_id)
    if not mutated:
        return df if mode == "exclude" else df.iloc[0:0]  # include with no mutations = empty

    if mode == "include":
        return df[df["ach_id"].isin(mutated)]
    else:  # exclude
        return df[~df["ach_id"].isin(mutated)]


def apply_fusion_filter(
    df: pd.DataFrame, mode: str, hugo: str
) -> pd.DataFrame:
    """Apply hard fusion filter: include / exclude / ignore.
    Matches fusion rows by HUGO symbol (see get_fused_cell_lines)."""
    if mode == "ignore":
        return df
    fused = data_service.get_fused_cell_lines(hugo)
    if not fused:
        return df if mode == "exclude" else df.iloc[0:0]

    if mode == "include":
        return df[df["ach_id"].isin(fused)]
    else:
        return df[~df["ach_id"].isin(fused)]

# ── Instability / Metabolite / miRNA filters ───────────────────────────

def apply_instability_filter(
    df: pd.DataFrame, msi_max: float | None, cin_max: float | None
) -> pd.DataFrame:
    """Exclude cell lines above the given MSI / CIN thresholds."""
    if msi_max is None and cin_max is None:
        return df
    unstable = data_service.get_unstable_cell_lines(msi_max, cin_max)
    return df[~df["ach_id"].isin(unstable)] if unstable else df


def apply_metabolite_filter(
    df: pd.DataFrame, metabolite: str | None, threshold: float | None
) -> pd.DataFrame:
    """Exclude cell lines whose level of a metabolite exceeds the threshold."""
    if metabolite is None or threshold is None:
        return df
    high = data_service.get_high_metabolite_cell_lines(metabolite, threshold)
    return df[~df["ach_id"].isin(high)] if high else df


def apply_mirna_filter(
    df: pd.DataFrame, mirna_id: str | None, threshold: float | None
) -> pd.DataFrame:
    """Exclude cell lines whose expression of a miRNA exceeds the threshold."""
    if mirna_id is None or threshold is None:
        return df
    high = data_service.get_high_mirna_cell_lines(mirna_id, threshold)
    return df[~df["ach_id"].isin(high)] if high else df


def _score_flags(
    flags_per_gene: dict[str, dict[str, list[str]]],
    ach_id: str,
    genes_skipped: list[str] | None = None,
) -> list[str] | None:
    """
    Collect one cell line's flags over all genes (gene-prefixed if several):
    MAD=0 notes, plus genes skipped because the cell line has no score for them.
    """
    multi = len(flags_per_gene) > 1
    out = [
        f"{hugo}: {text}" if multi else text
        for hugo, per_cell in flags_per_gene.items()
        for text in per_cell.get(ach_id, [])
    ]
    if isinstance(genes_skipped, list):
        out += [f"No data for {hugo}, skipped" for hugo in genes_skipped]
    return out or None


# ── Full ranking pipeline ────────────────────────────────────

def run_ranking(
    genes: list[dict],       # [{"hugo": "EGFR", "direction": "high", "ensembl_id": "ENSG..."}]
    w_rna: float = 0.7,
    w_protein: float = 0.3,
    mutation_mode: str = "ignore",
    fusion_mode: str = "ignore",
    msi_max: float | None = None,
    cin_max: float | None = None,
    exclude_metabolite: str | None = None,
    metabolite_threshold: float | None = None,
    exclude_mirna: str | None = None,
    mirna_threshold: float | None = None,
    disease_filter: str | None = None,
    lineage_filter: str | None = None,
    subtype_filter: str | None = None,
    target_disease: str | None = None,
    target_lineage: str | None = None,
    target_subtype: str | None = None,
    q6_boost_enabled: bool = True,
    core_only: bool = False,
    require_protein: bool = False,
    top_n: int = 20,
    scoring_method: str = "robust_z",
    sources: list[str] | None = None,
    assay_type: str | None = None,
    include_geo: bool = False,
    multi_gene_method: str = "geometric",
    low_expression_mode: str = DEFAULT_LOW_EXPRESSION_MODE,
    rna_method: str = "robust_z",
) -> list[dict]:
    """
    Execute the full ranking pipeline. Returns list of ranked cell line dicts.
    rna_method: robust_z (default), standard_z or log2_tpm, see data_service.RNA_METHODS.
    Only the RNA scale changes; protein, weights and the rest of the pipeline are the same.
    """
    # Per-gene scoring
    combined_per_gene: dict[str, pd.DataFrame] = {}
    coverage_per_gene: dict[str, pd.Series] = {}
    flags_per_gene: dict[str, dict[str, list[str]]] = {}
    categories_per_gene: dict[str, pd.DataFrame] = {}
    consensus_per_gene: dict[str, pd.Series] = {}

    if multi_gene_method not in MULTI_GENE_METHODS:
        raise ValueError(f"unknown multi_gene_method: {multi_gene_method!r}")
    if low_expression_mode not in LOW_EXPRESSION_MODES:
        raise ValueError(f"unknown low_expression_mode: {low_expression_mode!r}")
    if rna_method not in data_service.RNA_METHODS:
        raise ValueError(f"unknown rna_method: {rna_method!r}")
    # Geometric mean needs "higher = better" for every gene, so LOW genes are
    # inverted (1 - score). Single-gene search and summed RRF keep the old scores.
    invert_low = len(genes) > 1 and multi_gene_method == "geometric"

    # Only scoring RNA sources are used; "geo" and "protein" are dropped.
    # Selecting no scoring source (e.g. only "geo") means all of them.
    rna_sources = [s for s in (sources or []) if s in SCORING_SOURCE_NAMES] or SCORING_SOURCE_NAMES

    for gene in genes:
        ensg = gene["ensembl_id"]

        # 1. Expression score from the RNA consensus (only selected RNA sources)
        signal = data_service.get_rna_consensus(ensg, rna_sources, method=rna_method)
        expr_rrf = score_expression(signal, scoring_method)
        coverage = signal["coverage"].astype(int) if len(signal) >= 2 else pd.Series(dtype=int)
        coverage_per_gene[gene["hugo"]] = coverage

        # 1b. Expression category from raw TPM (shown for every gene). For a
        #     LOW gene it also sets the RNA score of off cell lines (priority 5).
        categories = data_service.classify_expression(signal)
        categories_per_gene[gene["hugo"]] = categories
        consensus_per_gene[gene["hugo"]] = signal["consensus"]
        if gene["direction"] == "low":
            expr_rrf = apply_low_expression_mode(expr_rrf, categories, low_expression_mode)

        # 2. Protein (always scored, independent of the RNA source filter)
        prot_rrf = rank_protein(ensg, scoring_method)
        # P6 LOW-mode fix: in mixed-coverage off tier, tie protein too so
        # RNA-only cells do not get a spurious advantage from missing protein.
        if gene["direction"] == "low":
            prot_rrf = apply_low_expression_to_protein(prot_rrf, categories, low_expression_mode)
        flags = _gene_flags(signal, prot_rrf)
        for ach, texts in _disagree_flags(categories).items():
            flags.setdefault(ach, []).extend(texts)
        flags_per_gene[gene["hugo"]] = flags

        # 3. Combine
        combined = combine_rna_protein(
            expr_rrf, prot_rrf, w_rna, w_protein,
            invert=invert_low and gene["direction"] == "low",
        )

        # 4. Mutation/Fusion filter
        combined = apply_mutation_filter(combined, mutation_mode, ensg)
        combined = apply_fusion_filter(combined, fusion_mode, gene["hugo"])
        combined = apply_instability_filter(combined, msi_max, cin_max)
        combined = apply_metabolite_filter(combined, exclude_metabolite, metabolite_threshold)
        combined = apply_mirna_filter(combined, exclude_mirna, mirna_threshold)

        combined_per_gene[gene["hugo"]] = combined

    # 5. Multi-gene combination
    # Sort tier per cell line: for a single LOW gene the expression category
    # comes first ("classify first, then rank"); otherwise every tier is 0.
    single_low = len(genes) == 1 and genes[0]["direction"] == "low"

    def tier_of(ach: str) -> int:
        if not single_low:
            return 0
        cats = categories_per_gene[genes[0]["hugo"]]
        if ach not in cats.index:
            return low_expression_tier(None, low_expression_mode)
        return low_expression_tier(
            cats.at[ach, "category"], low_expression_mode, bool(cats.at[ach, "sources_disagree"])
        )

    if len(genes) == 1:
        gene = genes[0]
        final_df = combined_per_gene[gene["hugo"]].copy()
        ascending = gene["direction"] == "low"
        final_df["tier"] = [tier_of(a) for a in final_df["ach_id"]]
        final_df = final_df.sort_values(
            ["tier", "combined_score", "ach_id"], ascending=[True, ascending, True], kind="mergesort"
        ).reset_index(drop=True)
        final_df["final_rank"] = range(1, len(final_df) + 1)
        final_df["final_score"] = final_df["combined_score"]
    else:
        gene_ranks: dict[str, pd.Series] = {}
        desirability: dict[str, pd.Series] = {}
        for gene in genes:
            df_g = combined_per_gene[gene["hugo"]].copy()
            # rank 1 = best for the requested direction (display + summed RRF)
            ascending = gene["direction"] == "low" and not invert_low
            df_g = _sort_by_score(df_g, "combined_score", ascending)
            df_g["gene_rank"] = range(1, len(df_g) + 1)
            gene_ranks[gene["hugo"]] = df_g.set_index("ach_id")["gene_rank"]
            desirability[gene["hugo"]] = df_g.set_index("ach_id")["combined_score"]

        if multi_gene_method == "geometric":
            # Weighted geometric mean of per-gene desirability: one failed gene
            # meaningfully lowers the score. Missing gene = skipped (see helper).
            weights = {g["hugo"]: float(g.get("weight", 1.0)) for g in genes}
            multi_rrf = weighted_geometric_mean(pd.DataFrame(desirability), weights)
        else:
            # Legacy: raw-sum RRF of per-gene ranks, divided by the theoretical
            # maximum (each gene at rank 1 -> n_genes x 1/(k+1)).
            multi_rrf = _rrf(gene_ranks, average=False) / (len(genes) / (RRF_K + 1))

        rows = []
        for ach in multi_rrf.index:
            coverages = []
            scenarios = []
            for gene in genes:
                gdf = combined_per_gene[gene["hugo"]]
                match = gdf[gdf["ach_id"] == ach]
                if len(match) > 0:
                    coverages.append(match["modality_coverage"].iloc[0])
                    scenarios.append(match["scenario"].iloc[0])

            rows.append({
                "ach_id": ach,
                "final_score": float(multi_rrf[ach]),
                "modality_coverage": float(np.mean(coverages)) if coverages else 0.0,
                "scenario": " | ".join(scenarios),
                # genes this cell line has no score for (skipped in the geometric mean)
                "genes_skipped": [
                    g["hugo"] for g in genes if ach not in desirability[g["hugo"]].index
                ],
            })

        final_df = pd.DataFrame(
            rows, columns=["ach_id", "final_score", "modality_coverage", "scenario", "genes_skipped"]
        )
        final_df = _sort_by_score(final_df, "final_score", ascending=False)
        final_df["final_rank"] = range(1, len(final_df) + 1)

    if final_df.empty:
        return []

    # 6. Merge metadata
    ach_ids = final_df["ach_id"].tolist()
    meta = data_service.get_cell_lines_metadata(ach_ids)

    # 7. Apply post-filters
    results = []
    for _, row in final_df.iterrows():
        cl = meta.get(row["ach_id"], {})

        # Disease filter
        if disease_filter and (cl.get("primary_disease") or "").lower() != disease_filter.lower():
            continue

        # Lineage filter
        if lineage_filter and (cl.get("lineage") or "").lower() != lineage_filter.lower():
            continue

        # Subtype filter
        if subtype_filter and (cl.get("subtype") or "").lower() != subtype_filter.lower():
            continue

        # Core only: require data in every SELECTED RNA source (for first gene).
        # Comparing with all sources made core_only return nothing when the
        # user picked a single source, since coverage can then be at most 1.
        is_core = False
        first_hugo = genes[0]["hugo"]
        if first_hugo in coverage_per_gene:
            cov = coverage_per_gene[first_hugo]
            is_core = bool(row["ach_id"] in cov.index and cov[row["ach_id"]] >= len(rna_sources))

        if core_only and not is_core:
            continue

        # Require protein: at least one searched gene has protein data
        # (how this should work for several genes is still to be decided).
        if require_protein and "Protein" not in str(row.get("scenario", "")):
            continue

        # Collect mutations for include mode
        mutations = None
        if mutation_mode == "include":
            all_muts = []
            for gene in genes:
                gene_muts = data_service.get_mutations_for_gene(gene["ensembl_id"])
                all_muts.extend([m for m in gene_muts if m["ach_id"] == row["ach_id"]])
            mutations = all_muts if all_muts else None

        # Collect fusions for include mode
        fusions = None
        if fusion_mode == "include":
            all_fusions = []
            for gene in genes:
                gene_fusions = data_service.get_fusions_for_gene(gene["hugo"])
                all_fusions.extend([f for f in gene_fusions if f["ach_id"] == row["ach_id"]])
            fusions = all_fusions if all_fusions else None

        score_col = "final_score" if "final_score" in row.index else "combined_score"

        # Per-gene ranks for multi-gene searches
        per_gene_rank = None
        if len(genes) > 1:
            per_gene_rank = {}
            for gene in genes:
                gr = gene_ranks.get(gene["hugo"])
                if gr is not None and row["ach_id"] in gr.index:
                    per_gene_rank[gene["hugo"]] = int(gr[row["ach_id"]])

        results.append({
            "ach_id": row["ach_id"],
            "cell_line_name": cl.get("cell_line_name"),
            "primary_disease": cl.get("primary_disease"),
            "lineage": cl.get("lineage"),
            "growth_pattern": cl.get("growth_pattern"),
            "score": round(float(row[score_col]), 4),
            "modality_coverage": round(float(row["modality_coverage"]), 2),
            "scenario": row.get("scenario", ""),
            "is_core": is_core,
            "mutations": mutations,
            "fusions": fusions,
            "per_gene_rank": per_gene_rank,
            "score_flags": _score_flags(flags_per_gene, row["ach_id"], row.get("genes_skipped")),
            # Raw TPM + expression category per gene (priority 5)
            "expression": {
                g["hugo"]: _expression_info(
                    categories_per_gene[g["hugo"]], row["ach_id"], consensus_per_gene[g["hugo"]]
                )
                for g in genes
            },
            "_tier": tier_of(row["ach_id"]),
        })

    # Determine sort direction for post-filter re-ranking.
    # For single gene with "low" direction, lower score = better → sort ascending.
    # For multi-gene, direction is already baked into RRF ranks → sort descending.
    _sort_descending = not single_low

    def sort_key(r: dict):
        # Tier first (expression category for a single LOW gene), then score.
        return (r["_tier"], -r["score"] if _sort_descending else r["score"])

    # Re-rank after filters (stable sort: equal keys keep the ach_id tie order)
    results.sort(key=sort_key)
    for i, r in enumerate(results):
        r["rank"] = i + 1

    # Q6 soft boost
    if q6_boost_enabled and (target_disease or target_lineage or target_subtype):
        from app.services.q6_lineage_scorer import compute_q6_score

        dim_full = data_service.get_full_dim()
        q6 = compute_q6_score(
            dim_full,
            query_lineage=target_lineage,
            query_disease=target_disease,
            query_subtype=target_subtype
        )
        q6_map = q6.set_index("ach_id")[["q6_score", "match_level"]].to_dict("index")

        for r in results:
            info = q6_map.get(r["ach_id"], {})
            q6_score = info.get("q6_score", 0.0)
            r["q6_score"] = q6_score
            r["match_level"] = info.get("match_level", "none")
            r["base_score"] = r["score"]
            if _sort_descending:
                # "high" direction: reward matching → higher score = better
                # Formula: final = base × (0.7 + 0.3 × q6)   (α = 0.3)
                r["score"] = round(r["score"] * (0.7 + 0.3 * q6_score), 4)
            else:
                # "low" direction: penalise non-matching → inflate their score
                # Matching (q6=1) → score × 1.0 (unchanged, stays near top)
                # No match (q6=0) → score × 1.3 (pushed down in ascending sort)
                r["score"] = round(r["score"] * (1.0 + 0.3 * (1.0 - q6_score)), 4)

        # Re-sort by boosted score (respecting direction and tier)
        results.sort(key=sort_key)
        for i, r in enumerate(results):
            r["rank"] = i + 1

    results = results[:top_n]
    for r in results:
        r.pop("_tier", None)

    # GEO support for the returned cell lines (extra evidence only:
    # it is added after ranking and never changes score or rank)
    if include_geo and results:
        geo_per_gene = {}
        for gene in genes:
            geo = data_service.get_geo_support(gene["ensembl_id"], gene["direction"])
            if geo is not None:
                geo_per_gene[gene["hugo"]] = geo[0]
        if geo_per_gene:
            for r in results:
                r["geo"] = {
                    hugo: data_service.geo_summary_for_cell_line(summary, r["ach_id"])
                    for hugo, summary in geo_per_gene.items()
                }

    # Q7 assay compatibility warnings (does NOT change ranking order)
    if assay_type and assay_type != "none":
        from app.services.q7_warning_checker import check_q7_warnings

        dim = pd.DataFrame(results)[["ach_id", "cell_line_name", "lineage", "growth_pattern"]]
        q7 = check_q7_warnings(dim, assay_type)
        q7_map = q7.set_index("ach_id")[["q7_status", "q7_warning"]].to_dict("index")
        for r in results:
            info = q7_map.get(r["ach_id"], {})
            r["q7_status"] = info.get("q7_status", "OK")
            r["q7_warning"] = info.get("q7_warning", "")

    return results