"""
Data access layer — DuckDB queries for gene lookup, expression, mutations, etc.

Note: Parquet files use categorical encoding for string columns (ach_id, ensembl_id,
hugo_symbol, etc.). DuckDB reads these as DICTIONARY type, which can cause
ConversionErrors with parameterised queries. We CAST to VARCHAR where needed.
"""

import hashlib
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy import stats

from app.config import (
    DEFAULT_OFF_METHOD, EXPR_DETECTION_TPM, EXPR_MODERATE_TPM, EXPR_OFF_TPM, OFF_DISTANCE_SD,
    OFF_METHODS, OFF_MIN_GROUP_SHARE, OFF_MIN_LINES, OFF_MIN_LOG2_GAP, OFF_MIXTURE_POSTERIOR,
    OFF_REF_FRACTION, OFF_REF_QUANTILE, OFF_STRICT_MIN_BIC_GAIN, OFF_STRICT_MIN_D,
    OFF_STRICT_MIN_GROUP_SHARE, OFF_STRICT_POSTERIOR, RRF_K,
)
from app.database import query, query_df


# ── Gene queries ──────────────────────────────────────────────

def search_genes(q: str, limit: int = 10) -> list[dict]:
    """Fuzzy search genes by HUGO symbol prefix."""
    return query(
        """
        SELECT DISTINCT
            CAST(hugo_symbol AS VARCHAR) AS hugo_symbol,
            CAST(ensembl_id AS VARCHAR) AS ensembl_id
        FROM dim_genes
        WHERE UPPER(CAST(hugo_symbol AS VARCHAR)) LIKE UPPER(? || '%')
        ORDER BY hugo_symbol
        LIMIT ?
        """,
        [q, limit],
    )


def resolve_gene(hugo: str) -> dict | None:
    """Resolve a HUGO symbol to its Ensembl ID."""
    rows = query(
        """
        SELECT
            CAST(hugo_symbol AS VARCHAR) AS hugo_symbol,
            CAST(ensembl_id AS VARCHAR) AS ensembl_id
        FROM dim_genes
        WHERE UPPER(CAST(hugo_symbol AS VARCHAR)) = UPPER(?)
        LIMIT 1
        """,
        [hugo],
    )
    return rows[0] if rows else None


# ── Expression queries ────────────────────────────────────────

def get_expression_for_gene(ensembl_id: str, source_table: str):
    """Get TPM values for a gene from one expression source. Returns DataFrame."""
    return query_df(
        f"""
        SELECT CAST(ach_id AS VARCHAR) AS ach_id, AVG(tpm) AS tpm
        FROM {source_table}
        WHERE CAST(ensembl_id AS VARCHAR) = ?
        GROUP BY ach_id
        """,
        [ensembl_id],
    )


# DepMap and HPA RNA values derive from the same underlying data, so they are
# NOT independent evidence. Each view is standardised separately and then
# combined ONCE into a single RNA consensus. (GEO is excluded: heterogeneous
# studies, see ranking_service.)
RNA_SOURCE_TABLES = {
    "depmap": "fact_expression_depmap",
    "hpa": "fact_expression_hpa",
}

_RNA_SIGNAL_COLUMNS = ["consensus", "coverage", "discordance"]
RNA_SOURCE_LABELS = {"depmap": "DepMap", "hpa": "HPA"}

# How expression values are prepared before the robust z-score is computed:
#   "tpm"  -> raw TPM  (current choice: the supervisor asked for TPM, not log TPM;
#                       asinh in the score stops a few very high cell lines from
#                       squashing the rest)
#   "log2" -> log2(TPM + 1)
# Still switched in ONE place.
RNA_TRANSFORM = "tpm"


def _prepare(tpm: pd.Series, transform: str) -> pd.Series:
    """Prepare TPM values for the robust z-score according to `transform`."""
    tpm = tpm.astype(float).clip(lower=0)
    if transform == "tpm":
        return tpm
    if transform == "log2":
        return np.log2(tpm + 1.0)
    raise ValueError(f"unknown RNA transform: {transform!r} (use 'tpm' or 'log2')")


def _expressed_reference(tpm: pd.Series, transform: str) -> pd.Series:
    """
    Values used to estimate the typical level and spread of a gene (median, MAD).
    Cell lines below EXPR_OFF_TPM (the "gene is off" cutoff, 1 TPM) are removed
    here, so noise below 1 TPM does not set the scale. Without this, DepMap's tiny
    continuous values (0.01, 0.02, ...) give a tiny non-zero MAD and z in the
    thousands, while HPA (0.1 steps, many exact zeros) has MAD = 0, so the two
    sources disagreed by orders of magnitude for the same cell line.
    The median and MAD are therefore those of the expressed cell lines. The z-score
    itself is still computed for every cell line from its real value, so the order of
    cell lines is unchanged and the off lines get a z below the expressed ones.
    """
    tpm = tpm.astype(float).clip(lower=0)
    return _prepare(tpm[tpm >= EXPR_OFF_TPM], transform)


def _robust_z(values: pd.Series, reference: pd.Series | None = None) -> pd.Series:
    """
    Gene-wise robust z-score across cell lines: (x - median) / (1.4826 * MAD).
    If MAD is 0 (most reference values tie) fall back to the standard deviation of the
    reference. If that is 0 too (the gene is constant or has no reference values),
    every z is 0.
    `reference` (any index) is used for the median, MAD and SD instead of the values
    themselves (see _expressed_reference).
    """
    ref = values if reference is None else reference
    med = ref.median()
    scale = (ref - med).abs().median() * 1.4826
    if not np.isfinite(scale) or scale == 0:
        scale = ref.std(ddof=0)
    if not np.isfinite(scale) or scale == 0:
        return pd.Series(0.0, index=values.index)
    return (values - med) / scale


# How the RNA of one source is put on a common scale (the same for every gene):
#   robust_z   (x - median) / (1.4826 * MAD)   the default
#   standard_z (x - mean) / SD
#   log2_tpm   log2(TPM + 1); no z-score, so a cell line's value does not depend on the others
RNA_METHODS = ("robust_z", "standard_z", "log2_tpm")


def _standard_z(values: pd.Series, reference: pd.Series | None = None) -> pd.Series:
    """
    Gene-wise standard z-score across cell lines: (x - mean) / SD, with the mean and SD of
    `reference` (see _expressed_reference), so cell lines below 1 TPM do not set the scale.
    If the SD is 0 or there are no reference values, every z is 0.
    """
    ref = values if reference is None else reference
    mean, sd = ref.mean(), ref.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=values.index)
    return (values - mean) / sd


def _mad_is_zero(values: pd.Series) -> bool:
    """True when MAD = 0, i.e. the robust z had to fall back to the SD."""
    vals = values.dropna()
    if len(vals) == 0:
        return False
    return bool((vals - vals.median()).abs().median() == 0)


def robust_z_with_flag(values: pd.Series) -> tuple[pd.Series, bool]:
    """Robust z of `values` plus whether MAD = 0 (SD fallback), for reuse (protein)."""
    return _robust_z(values), _mad_is_zero(values)


def get_rna_consensus(
    ensembl_id: str,
    enabled_sources: list[str] | None = None,
    transform: str | None = None,
    method: str = "robust_z",
) -> pd.DataFrame:
    """
    One provenance-aware RNA signal per cell line.

    Per source: prepare TPM (see RNA_TRANSFORM) -> gene-wise score across cell lines.
    `method` (see RNA_METHODS) says which score: robust z (default), standard z, or
    log2(TPM + 1). For both z-scores the reference (median and MAD, or mean and SD) comes
    from the cell lines at or above 1 TPM only (those below are removed from the reference,
    see _expressed_reference) and is then used for every cell line of that source.
    log2(TPM + 1) has no reference: every cell line keeps its own log2 value.
      consensus   = mean of the matched per-source scores (sources that cover the cell line)
      coverage    = number of sources used
      discordance = |z_depmap - z_hpa| when 2+ sources cover the cell line, else NaN
                    (QC only: it never changes the score or the ranking)

    Returns a DataFrame indexed by ach_id with columns
    consensus, coverage, discordance, z_<source> (the per-source score of `method`) and
    tpm_<source> (raw TPM, used for the expression category) for each source used.
    df.attrs["mad_zero_sources"] lists the sources where MAD = 0 for this gene
    (robust z then used the SD instead; surfaced to the user as a flag).
    df.attrs["rna_method"] is `method`.
    """
    if method not in RNA_METHODS:
        raise ValueError(f"unknown rna_method: {method!r} (use one of {', '.join(RNA_METHODS)})")
    transform = transform or RNA_TRANSFORM
    zcols: dict[str, pd.Series] = {}
    tpmcols: dict[str, pd.Series] = {}
    mad_zero_sources: list[str] = []
    for name, table in RNA_SOURCE_TABLES.items():
        if enabled_sources is not None and name not in enabled_sources:
            continue
        df = get_expression_for_gene(ensembl_id, table)
        if len(df) < 2:
            continue
        tpm = df.set_index("ach_id")["tpm"].astype(float)
        if method == "log2_tpm":
            zcols[f"z_{name}"] = _prepare(tpm, "log2")
        else:
            prepared = _prepare(tpm, transform)
            reference = _expressed_reference(tpm, transform)
            if method == "standard_z":
                zcols[f"z_{name}"] = _standard_z(prepared, reference)
            else:
                if _mad_is_zero(reference):
                    mad_zero_sources.append(name)
                zcols[f"z_{name}"] = _robust_z(prepared, reference)
        tpmcols[f"tpm_{name}"] = tpm.clip(lower=0)

    if not zcols:
        empty = pd.DataFrame(columns=_RNA_SIGNAL_COLUMNS, dtype=float)
        empty.attrs["mad_zero_sources"] = []
        empty.attrs["rna_method"] = method
        return empty

    zmat = pd.DataFrame(zcols)
    coverage = zmat.notna().sum(axis=1)
    spread = zmat.max(axis=1) - zmat.min(axis=1)
    out = pd.DataFrame({
        "consensus": zmat.mean(axis=1, skipna=True),
        "coverage": coverage,
        "discordance": spread.where(coverage >= 2),
    })
    result = pd.concat([out, zmat, pd.DataFrame(tpmcols)], axis=1)
    result.attrs["mad_zero_sources"] = mad_zero_sources
    result.attrs["rna_method"] = method
    return result


# ── Expression categories (AZ priority 5) ─────────────────────

# Ordered from "not expressed" to "expressed". "missing" = no RNA value at all.
# low_relative = 1 TPM or more, but off for this gene under the chosen off
# cutoff method (gene-relative; never produced by the default "fixed").
EXPRESSION_CATEGORIES = ("below_detection", "negligible", "low_relative", "low_measurable", "moderate")
OFF_CATEGORIES = ("below_detection", "negligible", "low_relative")
# Below 1 TPM in every source: differences here are noise, so these lines tie
# on RNA (and protein) in LOW searches. low_relative lines are off for their
# gene but at >= 1 TPM, so their real TPM differences (e.g. 2 vs 20) still count.
BELOW_1_CATEGORIES = ("below_detection", "negligible")


# ── Off cutoff methods (gene-relative "off") ─────────────────

_GROUP_FIT_CACHE: "OrderedDict[str, dict]" = OrderedDict()
_GROUP_FIT_CACHE_SIZE = 256


def _two_group_fit(log_level: pd.Series) -> dict:
    """
    Fit 1 and 2 Gaussian groups to log2(TPM + 1) of one gene across cell lines.
    Cached by the data itself (a ranking request, and the sensitivity sweep,
    classify the same gene several times).
    Returns the groups (low / high mean and SD in log2 units, the smallest group
    share), Ashman's D (how well the groups separate), the BIC gain of 2 groups
    over 1, and P(low group) for every cell line.
    """
    from sklearn.mixture import GaussianMixture   # only needed for these methods

    values = np.ascontiguousarray(log_level.to_numpy(dtype=float))
    key = hashlib.sha1(values.tobytes()).hexdigest()
    if key in _GROUP_FIT_CACHE:
        _GROUP_FIT_CACHE.move_to_end(key)
        fit = _GROUP_FIT_CACHE[key]
    else:
        x = values.reshape(-1, 1)
        one = GaussianMixture(1, random_state=0).fit(x)
        two = GaussianMixture(2, random_state=0).fit(x)
        mu = two.means_.ravel()
        sd = np.sqrt(two.covariances_.ravel())
        hi = int(np.argmax(mu))
        lo = 1 - hi
        fit = {
            "mu_lo": float(mu[lo]), "mu_hi": float(mu[hi]),
            "sd_lo": float(sd[lo]), "sd_hi": float(sd[hi]),
            "min_share": float(two.weights_.min()),
            "ashman_d": float(np.sqrt(2) * abs(mu[hi] - mu[lo]) / np.sqrt(sd[0] ** 2 + sd[1] ** 2)),
            "bic_gain": float(one.bic(x) - two.bic(x)),
            "p_low": two.predict_proba(x)[:, lo],
        }
        _GROUP_FIT_CACHE[key] = fit
        if len(_GROUP_FIT_CACHE) > _GROUP_FIT_CACHE_SIZE:
            _GROUP_FIT_CACHE.popitem(last=False)
    return {**fit, "p_low": pd.Series(fit["p_low"], index=log_level.index)}


def off_by_method(
    level: pd.Series, method: str = DEFAULT_OFF_METHOD, cutoff_tpm: float | None = None
) -> tuple[pd.Series, str | None]:
    """
    Which cell lines count as off for one gene.

    level: TPM per cell line (the highest source). Below EXPR_OFF_TPM a line is
    always off; `method` (see config.OFF_METHODS) can also call a line off at
    1 TPM or more when it is low for THIS gene.
    Returns (bool Series "off", reason) where reason explains, in plain words,
    why a line at >= 1 TPM was called off (None if the method added nothing).
    Data-driven methods fall back to the fixed rule when the gene has too few
    cell lines or no clear groups.
    """
    if method not in OFF_METHODS:
        raise ValueError(f"unknown off_method: {method!r} (use one of {', '.join(OFF_METHODS)})")
    level = level.astype(float)
    below_1 = level < EXPR_OFF_TPM
    if method == "fixed" or (method != "custom" and level.notna().sum() < OFF_MIN_LINES):
        return below_1, None

    if method == "custom":
        if cutoff_tpm is None:
            raise ValueError("off_method 'custom' needs off_cutoff_tpm")
        cut = max(float(cutoff_tpm), EXPR_OFF_TPM)
        return level < cut, f"below the chosen cutoff of {cut:g} TPM"

    if method == "percent_of_reference":
        expressing = level[level >= EXPR_OFF_TPM]
        if expressing.empty:
            return below_1, None
        ref = float(expressing.quantile(OFF_REF_QUANTILE))
        cut = max(OFF_REF_FRACTION * ref, EXPR_OFF_TPM)
        return level < cut, (
            f"below {cut:.1f} TPM = {OFF_REF_FRACTION:.0%} of this gene's reference level "
            f"({ref:.0f} TPM, the {OFF_REF_QUANTILE * 100:.0f}th percentile of expressing lines)"
        )

    log_level = np.log2(level.clip(lower=0) + 1.0)
    fit = _two_group_fit(log_level)
    lo_tpm, hi_tpm = 2 ** fit["mu_lo"] - 1, 2 ** fit["mu_hi"] - 1
    gap = fit["mu_hi"] - fit["mu_lo"]

    if method in ("mixture", "distance_from_on"):
        split = (
            gap >= OFF_MIN_LOG2_GAP
            and gap > max(fit["sd_lo"], fit["sd_hi"])
            and fit["min_share"] >= OFF_MIN_GROUP_SHARE
        )
        if not split:
            return below_1, None
        if method == "mixture":
            off = fit["p_low"] > OFF_MIXTURE_POSTERIOR
            reason = f"in the low group of a 2-group split (about {lo_tpm:.1f} vs {hi_tpm:.0f} TPM)"
        else:
            z = (log_level - fit["mu_hi"]) / fit["sd_hi"]
            off = z < -OFF_DISTANCE_SD
            reason = (f"more than {OFF_DISTANCE_SD:g} SD below the expressing group "
                      f"(typically {hi_tpm:.0f} TPM)")
        return below_1 | off, reason

    # strict_mixture: only genes with a clear on/off split
    clear = (
        fit["ashman_d"] > OFF_STRICT_MIN_D
        and fit["bic_gain"] > OFF_STRICT_MIN_BIC_GAIN
        and fit["min_share"] >= OFF_STRICT_MIN_GROUP_SHARE
    )
    if not clear:
        return below_1, None
    off = fit["p_low"] > OFF_STRICT_POSTERIOR
    return below_1 | off, (
        f"in the low group of a clear on/off split (about {lo_tpm:.1f} vs {hi_tpm:.0f} TPM, "
        f"separation D = {fit['ashman_d']:.1f})"
    )


def classify_expression(
    signal: pd.DataFrame, off_method: str = DEFAULT_OFF_METHOD, off_cutoff_tpm: float | None = None
) -> pd.DataFrame:
    """
    Expression category per cell line from raw TPM (get_rna_consensus output).

    A cell line gets the HIGHEST band any of its sources puts it in, so it only
    counts as off when every source agrees (conservative):
      below_detection : every source below its detection limit (EXPR_DETECTION_TPM)
      negligible      : every source below EXPR_OFF_TPM (1 TPM)
      low_relative    : 1 TPM or more, but off for this gene under `off_method`
                        (never with the default "fixed"; see off_by_method)
      low_measurable  : highest source from 1 to 3 TPM
      moderate        : highest source above EXPR_MODERATE_TPM (3 TPM)

    Returns a DataFrame indexed by ach_id with columns:
      category, tpm (mean over sources), tpm_<source>,
      sources_disagree (one source says off, another says expressed),
      off_reason (why a low_relative line is off, else None).
    Cell lines without any RNA value are not in the result ("missing").
    """
    tpm_cols = [c for c in signal.columns if c.startswith("tpm_")]
    if not tpm_cols:
        return pd.DataFrame(columns=["category", "tpm", "sources_disagree", "off_reason"])
    tpm = signal[tpm_cols].astype(float)
    have = tpm.notna()
    top = tpm.max(axis=1)

    # every available source below its own detection limit
    detect = pd.Series({c: EXPR_DETECTION_TPM.get(c[len("tpm_"):], 0.0) for c in tpm_cols})
    below_detect = (tpm.lt(detect, axis=1) | ~have).all(axis=1)

    category = np.select(
        [below_detect, top < EXPR_OFF_TPM, top <= EXPR_MODERATE_TPM],
        ["below_detection", "negligible", "low_measurable"],
        default="moderate",
    )
    off = tpm.lt(EXPR_OFF_TPM)
    on = tpm.ge(EXPR_OFF_TPM)
    out = pd.DataFrame({
        "category": category,
        "tpm": tpm.mean(axis=1),
        "sources_disagree": off.any(axis=1) & on.any(axis=1),
        "off_reason": None,
    }, index=signal.index)
    out = pd.concat([out, tpm], axis=1)
    out = out[have.any(axis=1)]

    # Gene-relative off: lines at >= 1 TPM that the chosen method calls off.
    if off_method != "fixed" and not out.empty:
        method_off, reason = off_by_method(top.loc[out.index], off_method, off_cutoff_tpm)
        relative = method_off & (top.loc[out.index] >= EXPR_OFF_TPM)
        out.loc[relative, "category"] = "low_relative"
        out.loc[relative, "off_reason"] = reason
    return out


# ── GEO support (extra evidence only, never part of the score) ──
# GEO samples come from different studies with different processing, so they are
# NOT pooled. For each study a gene-wise robust z is computed *within that study*
# (same transform and robust z as the RNA consensus). A cell line is "supportive"
# in a study when its z points the way the user asked for.
GEO_MIN_STUDY_SAMPLES = 10   # a study needs at least this many samples for the gene
GEO_SUPPORT_Z = 1.0          # z >= +1 (HIGH) / <= -1 (LOW) supports; the opposite contradicts

_GEO_MAP_CACHE: pd.DataFrame | None = None
_GEO_Z_CACHE: dict[str, pd.DataFrame] = {}   # per-gene results (ranking + detail panel ask for the same gene)
_GEO_Z_CACHE_MAX = 64
_GEO_COLUMNS = ["ach_id", "gse_id", "z", "n_samples"]


def _geo_study_map() -> pd.DataFrame | None:
    """GSM -> GSE map (dim_geo_samples.parquet). None if it has not been built."""
    global _GEO_MAP_CACHE
    if _GEO_MAP_CACHE is not None:
        return _GEO_MAP_CACHE
    try:
        df = query_df(
            """
            SELECT CAST(gsm_id AS VARCHAR) AS gsm_id, CAST(gse_id AS VARCHAR) AS gse_id
            FROM dim_geo_samples
            """
        )
    except Exception:
        return None
    _GEO_MAP_CACHE = df.dropna().drop_duplicates("gsm_id")
    return _GEO_MAP_CACHE


def _geo_samples_for_gene(ensembl_id: str) -> pd.DataFrame:
    """One TPM value per GEO sample (GSM) for a gene."""
    return query_df(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               CAST(GSM_ID AS VARCHAR) AS gsm_id,
               tpm
        FROM fact_expression_geo
        WHERE CAST(ensembl_id AS VARCHAR) = ?
        """,
        [ensembl_id],
    )


def get_geo_study_z(ensembl_id: str) -> pd.DataFrame | None:
    """
    Per (cell line, study) robust z of the gene, computed within each study.
    Columns: ach_id, gse_id, z, n_samples (samples in that study). Studies with
    fewer than GEO_MIN_STUDY_SAMPLES samples are skipped; a cell line with
    several samples in one study gets their mean z. Samples without a study are
    ignored. Returns None when the GEO data / study map is not available.
    """
    study_map = _geo_study_map()
    if study_map is None:
        return None
    if ensembl_id in _GEO_Z_CACHE:
        return _GEO_Z_CACHE[ensembl_id]
    try:
        samples = _geo_samples_for_gene(ensembl_id)
    except Exception:
        return None
    empty = pd.DataFrame(columns=_GEO_COLUMNS)
    if samples.empty:
        return empty

    df = samples.dropna(subset=["tpm"]).merge(study_map, on="gsm_id", how="inner")
    df = df.groupby(["gse_id", "gsm_id", "ach_id"], as_index=False)["tpm"].mean()
    rows = []
    for gse, g in df.groupby("gse_id"):
        g = g.drop_duplicates("gsm_id")
        if len(g) < GEO_MIN_STUDY_SAMPLES:
            continue
        tpm_g = pd.Series(g["tpm"].to_numpy(), index=g["gsm_id"].to_numpy())
        z = _robust_z(_prepare(tpm_g, RNA_TRANSFORM), _expressed_reference(tpm_g, RNA_TRANSFORM))
        per_cl = g.assign(z=z.reindex(g["gsm_id"]).to_numpy()).groupby("ach_id")["z"].mean()
        rows.append(pd.DataFrame({
            "ach_id": per_cl.index, "gse_id": gse, "z": per_cl.to_numpy(), "n_samples": len(g),
        }))
    result = pd.concat(rows, ignore_index=True) if rows else empty
    if len(_GEO_Z_CACHE) >= _GEO_Z_CACHE_MAX:
        _GEO_Z_CACHE.pop(next(iter(_GEO_Z_CACHE)))
    _GEO_Z_CACHE[ensembl_id] = result
    return result


def geo_verdict(z: float, direction: str) -> str:
    """One study: supportive / contradictory / neutral for the requested direction."""
    s = z if direction == "high" else -z
    if s >= GEO_SUPPORT_Z:
        return "supportive"
    if s <= -GEO_SUPPORT_Z:
        return "contradictory"
    return "neutral"


def classify_geo(study_z: pd.DataFrame, direction: str) -> pd.DataFrame:
    """
    One row per cell line with GEO data: geo_supportive / geo_contradictory /
    geo_total (studies) and geo_status:
      supportive     more studies support than contradict, and at least half support
      contradictory  more studies contradict than support, and at least half contradict
      mixed          anything else (including studies that are neutral)
    """
    cols = ["geo_status", "geo_supportive", "geo_contradictory", "geo_total"]
    if study_z is None or study_z.empty:
        return pd.DataFrame(columns=cols)
    s = study_z["z"] if direction == "high" else -study_z["z"]
    d = pd.DataFrame({
        "ach_id": study_z["ach_id"],
        "sup": (s >= GEO_SUPPORT_Z).astype(int),
        "con": (s <= -GEO_SUPPORT_Z).astype(int),
    })
    agg = d.groupby("ach_id").agg(
        geo_supportive=("sup", "sum"), geo_contradictory=("con", "sum"), geo_total=("sup", "size"),
    )
    sup, con, tot = agg["geo_supportive"], agg["geo_contradictory"], agg["geo_total"]
    agg["geo_status"] = np.where(
        (sup > con) & (sup * 2 >= tot), "supportive",
        np.where((con > sup) & (con * 2 >= tot), "contradictory", "mixed"),
    )
    return agg[cols]


def get_geo_support(ensembl_id: str, direction: str):
    """(per-cell-line summary, per-study z table) or None if GEO data is not available."""
    study_z = get_geo_study_z(ensembl_id)
    if study_z is None:
        return None
    return classify_geo(study_z, direction), study_z


def geo_summary_for_cell_line(summary: pd.DataFrame, ach_id: str) -> dict:
    """Plain dict for one cell line; 'unavailable' when it has no usable GEO study."""
    if ach_id not in summary.index:
        return {"status": "unavailable", "supportive": 0, "contradictory": 0, "total": 0}
    r = summary.loc[ach_id]
    return {
        "status": str(r["geo_status"]),
        "supportive": int(r["geo_supportive"]),
        "contradictory": int(r["geo_contradictory"]),
        "total": int(r["geo_total"]),
    }


# ── Protein queries ───────────────────────────────────────────

def get_protein_for_gene(ensembl_id: str):
    """Get protein intensity for a gene. Returns DataFrame."""
    return query_df(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               AVG(protein_intensity) AS protein_intensity
        FROM fact_proteomics
        WHERE CAST(ensembl_id AS VARCHAR) = ?
        GROUP BY ach_id
        """,
        [ensembl_id],
    )


# ── Mutation queries ──────────────────────────────────────────

def get_mutations_for_gene(ensembl_id: str) -> list[dict]:
    """Get all mutations for a gene."""
    return query(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               CAST(hugo_symbol AS VARCHAR) AS hugo_symbol,
               chrom, pos, ref, alt,
               variant_type, variant_info, protein_change,
               is_driver, is_hotspot, is_damaging, is_lof, allele_freq
        FROM fact_mutations
        WHERE CAST(ensembl_id AS VARCHAR) = ?
        """,
        [ensembl_id],
    )


def get_mutated_cell_lines(ensembl_id: str) -> set[str]:
    """Get set of ach_ids that have mutations in this gene."""
    rows = query(
        """
        SELECT DISTINCT CAST(ach_id AS VARCHAR) AS ach_id
        FROM fact_mutations
        WHERE CAST(ensembl_id AS VARCHAR) = ?
        """,
        [ensembl_id],
    )
    return {r["ach_id"] for r in rows}


# ── Fusion queries ────────────────────────────────────────────

def get_fusions_for_gene(hugo: str) -> list[dict]:
    """Get fusions involving a gene (as gene1 or gene2), any confidence level.
    Matches on HUGO symbol because the ensg columns in fact_fusions store
    'SYMBOL (ENSG.version)' strings, not bare Ensembl IDs."""
    return query(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               CAST(gene1_hugo AS VARCHAR) AS gene1_hugo,
               CAST(gene2_hugo AS VARCHAR) AS gene2_hugo,
               fusion_name, ffpm,
               CAST(confidence AS VARCHAR) AS confidence,
               reading_frame, supporting_reads
        FROM fact_fusions
        WHERE (CAST(gene1_hugo AS VARCHAR) = ? OR CAST(gene2_hugo AS VARCHAR) = ?)
        """,
        [hugo, hugo],
    )


def get_fused_cell_lines(hugo: str) -> set[str]:
    """Get set of ach_ids that have fusions involving this gene (any confidence).
    Matches on HUGO symbol (see note in get_fusions_for_gene)."""
    rows = query(
        """
        SELECT DISTINCT CAST(ach_id AS VARCHAR) AS ach_id
        FROM fact_fusions
        WHERE (CAST(gene1_hugo AS VARCHAR) = ? OR CAST(gene2_hugo AS VARCHAR) = ?)
        """,
        [hugo, hugo],
    )
    return {r["ach_id"] for r in rows}

# ── Non-gene-expression exclusion queries ─────────────────────

def get_unstable_cell_lines(msi_max: float | None, cin_max: float | None) -> set[str]:
    """Cell lines exceeding the given genomic instability thresholds.
    signatures_clean.parquet: index=ModelID (ACH IDs), cols: MSIScore, CIN, etc."""
    clauses = []
    params = []
    if msi_max is not None:
        clauses.append("MSIScore > ?")
        params.append(msi_max)
    if cin_max is not None:
        clauses.append("CIN > ?")
        params.append(cin_max)
    if not clauses:
        return set()

    rows = query(
        f"""
        SELECT DISTINCT CAST(ModelID AS VARCHAR) AS ach_id
        FROM fact_signatures
        WHERE {" OR ".join(clauses)}
        """,
        params,
    )
    return {r["ach_id"] for r in rows}


def get_high_metabolite_cell_lines(metabolite: str, threshold: float) -> set[str]:
    """Cell lines whose level of a metabolite exceeds the threshold.
    metabolomics_clean.parquet is wide: cols = CCLE_ID, DepMap_ID, + 225 metabolites."""
    import re
    safe_col = re.sub(r"[^a-zA-Z0-9_ ()/-]", "", metabolite)
    try:
        rows = query(
            f"""
            SELECT DISTINCT CAST(DepMap_ID AS VARCHAR) AS ach_id
            FROM fact_metabolomics
            WHERE "{safe_col}" > ?
            """,
            [threshold],
        )
    except Exception:
        return set()
    return {r["ach_id"] for r in rows}


def get_high_mirna_cell_lines(mirna_id: str, threshold: float) -> set[str]:
    """Cell lines whose expression of a miRNA exceeds the threshold.
    mirna_clean.parquet is transposed: rows = miRNAs, cols = ACH IDs.
    Read with pandas since DuckDB can't easily query this layout."""
    import pandas as pd
    from app.config import PARQUET

    path = PARQUET.get("fact_mirna")
    if path is None or not path.exists():
        return set()

    df = pd.read_parquet(path)
    # Find the string column holding miRNA names
    str_cols = df.select_dtypes(include="object").columns
    if len(str_cols) == 0:
        return set()

    mirna_col = str_cols[0]
    row = df[df[mirna_col] == mirna_id]
    if row.empty:
        return set()

    # Remaining columns are ACH IDs with expression values
    values = row.drop(columns=[mirna_col]).iloc[0]
    return set(values[values > threshold].index.astype(str))


# ── Cell line metadata ────────────────────────────────────────

def get_cell_line(ach_id: str) -> dict | None:
    """Get metadata for a single cell line."""
    rows = query(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               cell_line_name,
               CAST(primary_disease AS VARCHAR) AS primary_disease,
               CAST(lineage AS VARCHAR) AS lineage,
               CAST(lineage_subtype AS VARCHAR) AS lineage_subtype,
               CAST(sex AS VARCHAR) AS sex,
               CAST(growth_pattern AS VARCHAR) AS growth_pattern
        FROM dim_cell_lines
        WHERE CAST(ach_id AS VARCHAR) = ?
        LIMIT 1
        """,
        [ach_id],
    )
    return rows[0] if rows else None


def get_cell_lines_metadata(ach_ids: list[str]) -> dict[str, dict]:
    """Get metadata for multiple cell lines, returned as {ach_id: {...}}."""
    if not ach_ids:
        return {}
    placeholders = ", ".join(["CAST(? AS VARCHAR)"] * len(ach_ids))
    rows = query(
        f"""
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
               cell_line_name,
               CAST(primary_disease AS VARCHAR) AS primary_disease,
               CAST(lineage AS VARCHAR) AS lineage,
               CAST(lineage_subtype AS VARCHAR) AS lineage_subtype,
               CAST("Subtype" AS VARCHAR) AS subtype,
               CAST(sex AS VARCHAR) AS sex,
               CAST(growth_pattern AS VARCHAR) AS growth_pattern
        FROM dim_cell_lines
        WHERE CAST(ach_id AS VARCHAR) IN ({placeholders})
        """,
        ach_ids,
    )
    return {r["ach_id"]: r for r in rows}


# ── Filter options ────────────────────────────────────────────

def _normalize_disease(s: str) -> str:
    """Normalize disease casing: 'Bone cancer' → 'Bone Cancer'."""
    return " ".join(w.capitalize() for w in s.split())


def get_diseases() -> list[str]:
    """Get distinct disease values (case-normalized, no duplicates)."""
    rows = query(
        """
        SELECT DISTINCT CAST(primary_disease AS VARCHAR) AS primary_disease
        FROM dim_cell_lines
        WHERE primary_disease IS NOT NULL
          AND CAST(primary_disease AS VARCHAR) != ''
        """
    )
    seen = set()
    result = []
    for r in rows:
        d = r["primary_disease"]
        if not d:
            continue
        n = _normalize_disease(d)
        if n not in seen:
            seen.add(n)
            result.append(n)
    return sorted(result)


def get_lineages() -> list[str]:
    """Get distinct lineage values."""
    rows = query(
        """
        SELECT DISTINCT CAST(lineage AS VARCHAR) AS lineage
        FROM dim_cell_lines
        WHERE lineage IS NOT NULL
          AND CAST(lineage AS VARCHAR) != ''
        ORDER BY lineage
        """
    )
    return [r["lineage"] for r in rows if r["lineage"]]


def get_subtypes() -> list[str]:
    """Get distinct Subtype values."""
    rows = query(
        """
        SELECT DISTINCT CAST("Subtype" AS VARCHAR) AS subtype
        FROM dim_cell_lines
        WHERE "Subtype" IS NOT NULL
          AND CAST("Subtype" AS VARCHAR) != ''
        ORDER BY subtype
        """
    )
    return [r["subtype"] for r in rows if r["subtype"]]


def get_disease_lineage_mapping() -> list[dict]:
    """Return all distinct (disease, lineage) pairs for linked dropdown filtering."""
    rows = query(
        """
        SELECT DISTINCT
            CAST(primary_disease AS VARCHAR) AS disease,
            CAST(lineage AS VARCHAR) AS lineage
        FROM dim_cell_lines
        WHERE primary_disease IS NOT NULL
          AND CAST(primary_disease AS VARCHAR) != ''
          AND lineage IS NOT NULL
          AND CAST(lineage AS VARCHAR) != ''
        ORDER BY disease, lineage
        """
    )
    return [{"disease": _normalize_disease(r["disease"]), "lineage": r["lineage"]}
            for r in rows if r["disease"] and r["lineage"]]


def get_disease_lineage_subtype_mapping() -> list[dict]:
    """Return all distinct (disease, lineage, subtype) triples for cascading filters."""
    rows = query(
        """
        SELECT DISTINCT
            CAST(primary_disease AS VARCHAR) AS disease,
            CAST(lineage AS VARCHAR) AS lineage,
            CAST("Subtype" AS VARCHAR) AS subtype
        FROM dim_cell_lines
        WHERE primary_disease IS NOT NULL
          AND CAST(primary_disease AS VARCHAR) != ''
          AND lineage IS NOT NULL
          AND CAST(lineage AS VARCHAR) != ''
          AND "Subtype" IS NOT NULL
          AND CAST("Subtype" AS VARCHAR) != ''
        ORDER BY disease, lineage, subtype
        """
    )
    return [{"disease": _normalize_disease(r["disease"]),
             "lineage": r["lineage"],
             "subtype": r["subtype"]}
            for r in rows if r["disease"] and r["lineage"] and r["subtype"]]


# ── Full dim table (for Q6 scoring) ────────────────────────────

def get_full_dim() -> pd.DataFrame:
    """Return full dim_cell_lines including Subtype for Q6 scoring."""
    return query_df(
        """
        SELECT CAST(ach_id AS VARCHAR) AS ach_id,
            cell_line_name,
            CAST(primary_disease AS VARCHAR) AS primary_disease,
            CAST(lineage AS VARCHAR) AS lineage,
            CAST("Subtype" AS VARCHAR) AS Subtype,
            CAST(growth_pattern AS VARCHAR) AS growth_pattern
        FROM dim_cell_lines
        """
    )


# ── Evidence breakdown (per-source scoring) ──────────────────

EXPRESSION_SOURCES = [
    ("DepMap", "fact_expression_depmap"),
    ("HPA",    "fact_expression_hpa"),
]


def get_evidence_for_cell_line(
    ach_id: str,
    genes: list[dict],   # [{"hugo": "EGFR", "ensembl_id": "ENSG...", "direction": "high"}]
) -> dict:
    """
    Compute full evidence breakdown for a single cell line across all queried genes.

    Returns {
        "genes": {
            "EGFR": {
                "direction": "high",
                "sources": [
                    {"source": "DepMap", "tpm": 342.7, "z_score": 2.41,
                     "percentile": 97.8, "z_rank": 3, "pct_rank": 2, "total": 1840},
                    ...
                ],
                "source_count": 3,
                "expression_rank": 2,
                "expression_total": 1840,
            },
            ...
        },
        "protein": {"intensity": 4.82, "z_score": 1.87, "rank": 5, "total": 300},
        "mutations": [...],
        "fusions": [...],
        "data_coverage": {"rna_sources": 3, "has_protein": true, "has_mutations": true, "total_types": 5, "max_types": 5}
    }
    """
    result = {"genes": {}}

    all_mutations = []
    all_fusions = []
    has_protein = False
    total_rna_sources = 0

    for gene in genes:
        ensg = gene["ensembl_id"]
        hugo = gene["hugo"]
        direction = gene.get("direction", "high")
        ascending = direction == "low"

        gene_data = {
            "direction": direction,
            "sources": [],
            "source_count": 0,
            "expression_rank": None,
            "expression_total": 0,
        }

        # Per-source details (shown for provenance only; the overall expression
        # rank below uses the single RNA consensus)
        for src_name, src_table in EXPRESSION_SOURCES:
            df = get_expression_for_gene(ensg, src_table)
            if len(df) < 2:
                continue

            tpm_series = df.set_index("ach_id")["tpm"]
            total_cls = len(tpm_series)

            if ach_id not in tpm_series.index:
                continue

            tpm_val = float(tpm_series[ach_id])
            if pd.isna(tpm_val):
                continue

            # Z-score
            z_values = stats.zscore(tpm_series, nan_policy="omit")
            z_series = pd.Series(z_values, index=tpm_series.index) if not isinstance(z_values, pd.Series) else z_values
            z_val = float(z_series[ach_id])

            # Percentile
            pct_series = tpm_series.rank(pct=True, method="average")
            pct_val = float(pct_series[ach_id]) * 100

            # Z-Score rank (rank 1 = best for HIGH, rank 1 = lowest for LOW)
            z_rank_series = z_series.rank(ascending=ascending, method="min")
            z_rank_raw = z_rank_series[ach_id]
            z_rank = int(z_rank_raw) if pd.notna(z_rank_raw) else None

            # Percentile rank
            pct_rank_series = tpm_series.rank(ascending=ascending, method="min")
            pct_rank_raw = pct_rank_series[ach_id]
            pct_rank = int(pct_rank_raw) if pd.notna(pct_rank_raw) else None

            gene_data["sources"].append({
                "source": src_name,
                "tpm": round(tpm_val, 2),
                "z_score": round(z_val, 2),
                "percentile": round(pct_val, 1),
                "z_rank": z_rank,
                "pct_rank": pct_rank,
                "total": total_cls,
            })

        gene_data["source_count"] = len(gene_data["sources"])
        if gene_data["source_count"] > total_rna_sources:
            total_rna_sources = gene_data["source_count"]

        # Overall expression rank from the single RNA consensus (DepMap + HPA collapsed)
        signal = get_rna_consensus(ensg)
        if len(signal) >= 2 and ach_id in signal.index:
            consensus = signal["consensus"]
            # The rank only depends on the order of the consensus z, which is the
            # same for every scoring method.
            rank_c = consensus.rank(ascending=ascending, method="min")
            rrf_scores = 2.0 / (RRF_K + rank_c)
            rrf_ranked = rrf_scores.rank(ascending=False, method="min")
            rrf_val = rrf_ranked[ach_id]
            gene_data["expression_rank"] = int(rrf_val) if pd.notna(rrf_val) else None
            gene_data["expression_total"] = len(rrf_ranked)

            # Consensus value + discordance QC for this cell line
            row = signal.loc[ach_id]
            disc = row["discordance"]
            gene_data["rna_consensus"] = {
                "z": round(float(row["consensus"]), 2),
                "discordance": round(float(disc), 2) if pd.notna(disc) else None,
                "sources_used": int(row["coverage"]),
                # sources covering this cell line where MAD = 0 (SD fallback)
                "mad_zero_sources": [
                    RNA_SOURCE_LABELS.get(s, s)
                    for s in signal.attrs.get("mad_zero_sources", [])
                    if pd.notna(row.get(f"z_{s}"))
                ],
            }
            # Per-source robust z (provenance)
            for src in gene_data["sources"]:
                zc = f"z_{src['source'].lower()}"
                if zc in signal.columns and pd.notna(row.get(zc)):
                    src["robust_z"] = round(float(row[zc]), 2)

        # GEO support (extra evidence; does not change score or rank)
        geo = get_geo_support(ensg, direction)
        if geo is not None:
            geo_summary, geo_studies = geo
            gene_data["geo"] = geo_summary_for_cell_line(geo_summary, ach_id)
            mine = geo_studies[geo_studies["ach_id"] == ach_id].sort_values("gse_id")
            gene_data["geo"]["studies"] = [
                {
                    "gse_id": str(r.gse_id),
                    "z": round(float(r.z), 2),
                    "n_samples": int(r.n_samples),
                    "verdict": geo_verdict(float(r.z), direction),
                }
                for r in mine.itertuples()
            ]

        # Mutations for this gene + cell line
        gene_muts = get_mutations_for_gene(ensg)
        cl_muts = [m for m in gene_muts if m["ach_id"] == ach_id]
        all_mutations.extend(cl_muts)

        # Fusions for this gene + cell line
        gene_fusions = get_fusions_for_gene(hugo)
        cl_fusions = [f for f in gene_fusions if f["ach_id"] == ach_id]
        all_fusions.extend(cl_fusions)

        result["genes"][hugo] = gene_data

    # Protein (use first gene's ensembl_id for protein lookup)
    protein_data = None
    if genes:
        ensg = genes[0]["ensembl_id"]
        prot_df = get_protein_for_gene(ensg)
        if len(prot_df) >= 2:
            intensity_series = prot_df.set_index("ach_id")["protein_intensity"]
            if ach_id in intensity_series.index:
                intensity_val = float(intensity_series[ach_id])
                if pd.notna(intensity_val):
                    has_protein = True
                    total_cls = len(intensity_series)

                    # Z-score
                    z_values = stats.zscore(intensity_series, nan_policy="omit")
                    z_series = pd.Series(z_values, index=intensity_series.index) if not isinstance(z_values, pd.Series) else z_values
                    z_val = float(z_series[ach_id])

                    # Z-score rank
                    z_rank_series = z_series.rank(ascending=False, method="min")
                    z_rank_raw = z_rank_series[ach_id]
                    z_rank = int(z_rank_raw) if pd.notna(z_rank_raw) else None

                    # Robust z (what the protein score is computed from) + MAD=0 flag
                    rz_series, prot_mad_zero = robust_z_with_flag(intensity_series.astype(float))
                    rz_val = float(rz_series[ach_id])

                    # Percentile
                    pct_series = intensity_series.rank(pct=True, method="average")
                    pct_val = float(pct_series[ach_id]) * 100

                    # Percentile rank
                    pct_rank_series = intensity_series.rank(ascending=False, method="min")
                    pct_rank_raw = pct_rank_series[ach_id]
                    pct_rank = int(pct_rank_raw) if pd.notna(pct_rank_raw) else None

                    protein_data = {
                        "intensity": round(intensity_val, 2),
                        "z_score": round(z_val, 2) if pd.notna(z_val) else None,
                        "z_rank": z_rank,
                        "robust_z": round(rz_val, 2) if pd.notna(rz_val) else None,
                        "mad_zero": prot_mad_zero,
                        "percentile": round(pct_val, 1) if pd.notna(pct_val) else None,
                        "pct_rank": pct_rank,
                        "rank": z_rank,
                        "total": total_cls,
                    }

    result["protein"] = protein_data
    result["mutations"] = all_mutations
    result["fusions"] = all_fusions

    # Data coverage
    has_mutations = len(all_mutations) > 0
    has_fusions = len(all_fusions) > 0
    # RNA counts as ONE evidence type (DepMap + HPA are collapsed into one consensus)
    total_types = (1 if total_rna_sources > 0 else 0) + (1 if has_protein else 0) + (1 if has_mutations else 0)
    max_types = 1 + 1 + 1  # RNA consensus + protein + mutations

    result["data_coverage"] = {
        "rna_sources": total_rna_sources,
        "has_protein": has_protein,
        "has_mutations": has_mutations,
        "has_fusions": has_fusions,
        "total_types": total_types,
        "max_types": max_types,
    }

    return result
