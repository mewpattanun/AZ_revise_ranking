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

from app.config import RRF_K
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


def magnitude_score(consensus: pd.Series) -> pd.Series:
    """
    Score from the SIZE of the consensus robust z (keeps magnitude, unlike a rank).

    asinh(z) is ~z near 0 and ~log(2z) for big z, so an outlier at z = 60 does
    not squash every other cell line to ~0 when we min-max to [0, 1]. It is
    monotonic, so the order is the order of the consensus z (HIGH: biggest
    first, LOW: ascending, same as before).
    """
    return _norm01(np.arcsinh(consensus.astype(float)))


def score_expression(signal: pd.DataFrame, scoring_method: str = "robust_z") -> pd.Series:
    """
    Expression score per cell line from the RNA consensus signal.
    scoring_method:
      "robust_z"  → magnitude of the consensus robust z (default, point 3)
      "rrf" / "zscore" / "percentile" → legacy rank-based (RRF of z-rank / pct-rank)
    """
    if scoring_method not in SCORING_METHODS:
        raise ValueError(f"unknown scoring_method: {scoring_method!r}")
    if len(signal) < 2:
        return pd.Series(dtype=float)
    consensus = signal["consensus"]

    if scoring_method == "robust_z":
        return magnitude_score(consensus)

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
) -> tuple[pd.Series, pd.Series]:
    """
    Compute the expression score for one gene.
    enabled_sources: RNA sources that feed the consensus (e.g. ["depmap", "hpa"]).
                     None means all sources.
    Returns (scores indexed by ach_id, source_coverage Series).
    """
    # DepMap + HPA share the same underlying data → one RNA consensus per cell line
    # (mean of per-source robust z-scores; discordance is QC only, see data_service)
    signal = data_service.get_rna_consensus(ensembl_id, enabled_sources)
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
    Combine expression and protein scores.
    invert=True flips both to "low is good" (1 - score), used for LOW genes in
    multi-gene searches so that a higher score is always better.
    Returns DataFrame with: ach_id, combined_score, scenario, coverage_confidence.
    """
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
        has_rna = ach in expr_norm.index and not np.isnan(expr_norm.get(ach, np.nan))
        has_prot = ach in prot_norm.index and not np.isnan(prot_norm.get(ach, np.nan))

        if has_rna and has_prot:
            score = w_rna * expr_norm[ach] + w_prot * prot_norm[ach]
            scenario = "RNA+Protein"
            confidence = 1.0
        elif has_rna:
            score = w_rna * float(expr_norm[ach])
            scenario = "RNA only"
            confidence = w_rna
        elif has_prot:
            score = w_prot * float(prot_norm[ach])
            scenario = "Protein only"
            confidence = w_prot
        else:
            continue

        rows.append({
            "ach_id": ach,
            "combined_score": score,
            "scenario": scenario,
            "coverage_confidence": confidence,
        })

    return pd.DataFrame(rows)


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
    top_n: int = 20,
    scoring_method: str = "robust_z",
    sources: list[str] | None = None,
    assay_type: str | None = None,
    include_geo: bool = False,
    multi_gene_method: str = "geometric",
) -> list[dict]:
    """
    Execute the full ranking pipeline. Returns list of ranked cell line dicts.
    """
    # Per-gene scoring
    combined_per_gene: dict[str, pd.DataFrame] = {}
    coverage_per_gene: dict[str, pd.Series] = {}
    flags_per_gene: dict[str, dict[str, list[str]]] = {}

    if multi_gene_method not in MULTI_GENE_METHODS:
        raise ValueError(f"unknown multi_gene_method: {multi_gene_method!r}")
    # Geometric mean needs "higher = better" for every gene, so LOW genes are
    # inverted (1 - score). Single-gene search and summed RRF keep the old scores.
    invert_low = len(genes) > 1 and multi_gene_method == "geometric"

    for gene in genes:
        ensg = gene["ensembl_id"]

        # 1. Expression score from the RNA consensus (only selected RNA sources)
        rna_sources = [s for s in (sources or ["depmap", "hpa"]) if s != "protein"]
        signal = data_service.get_rna_consensus(ensg, rna_sources)
        expr_rrf = score_expression(signal, scoring_method)
        coverage = signal["coverage"].astype(int) if len(signal) >= 2 else pd.Series(dtype=int)
        coverage_per_gene[gene["hugo"]] = coverage

        # 2. Protein (always scored, independent of the RNA source filter)
        prot_rrf = rank_protein(ensg, scoring_method)
        flags_per_gene[gene["hugo"]] = _gene_flags(signal, prot_rrf)

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
    if len(genes) == 1:
        gene = genes[0]
        final_df = combined_per_gene[gene["hugo"]].copy()
        ascending = gene["direction"] == "low"
        final_df = final_df.sort_values("combined_score", ascending=ascending).reset_index(drop=True)
        final_df["final_rank"] = range(1, len(final_df) + 1)
        final_df["final_score"] = final_df["combined_score"]
    else:
        gene_ranks: dict[str, pd.Series] = {}
        desirability: dict[str, pd.Series] = {}
        for gene in genes:
            df_g = combined_per_gene[gene["hugo"]].copy()
            # rank 1 = best for the requested direction (display + summed RRF)
            ascending = gene["direction"] == "low" and not invert_low
            df_g = df_g.sort_values("combined_score", ascending=ascending).reset_index(drop=True)
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
            confidences = []
            scenarios = []
            for gene in genes:
                gdf = combined_per_gene[gene["hugo"]]
                match = gdf[gdf["ach_id"] == ach]
                if len(match) > 0:
                    confidences.append(match["coverage_confidence"].iloc[0])
                    scenarios.append(match["scenario"].iloc[0])

            rows.append({
                "ach_id": ach,
                "final_score": float(multi_rrf[ach]),
                "coverage_confidence": float(np.mean(confidences)) if confidences else 0.0,
                "scenario": " | ".join(scenarios),
                # genes this cell line has no score for (skipped in the geometric mean)
                "genes_skipped": [
                    g["hugo"] for g in genes if ach not in desirability[g["hugo"]].index
                ],
            })

        final_df = pd.DataFrame(rows)
        final_df = final_df.sort_values("final_score", ascending=False).reset_index(drop=True)
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

        # Core only: require data in all expression sources (for first gene)
        is_core = False
        first_hugo = genes[0]["hugo"]
        if first_hugo in coverage_per_gene:
            cov = coverage_per_gene[first_hugo]
            is_core = row["ach_id"] in cov.index and cov[row["ach_id"]] >= len(data_service.RNA_SOURCE_TABLES)

        if core_only and not is_core:
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
            "confidence": round(float(row["coverage_confidence"]), 2),
            "scenario": row.get("scenario", ""),
            "is_core": is_core,
            "mutations": mutations,
            "fusions": fusions,
            "per_gene_rank": per_gene_rank,
            "score_flags": _score_flags(flags_per_gene, row["ach_id"], row.get("genes_skipped")),
        })

    # Determine sort direction for post-filter re-ranking.
    # For single gene with "low" direction, lower score = better → sort ascending.
    # For multi-gene, direction is already baked into RRF ranks → sort descending.
    _sort_descending = not (len(genes) == 1 and genes[0]["direction"] == "low")

    # Re-rank after filters
    results.sort(key=lambda x: x["score"], reverse=_sort_descending)
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

        # Re-sort by boosted score (respecting direction)
        results.sort(key=lambda x: x["score"], reverse=_sort_descending)
        for i, r in enumerate(results):
            r["rank"] = i + 1

    results = results[:top_n]

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