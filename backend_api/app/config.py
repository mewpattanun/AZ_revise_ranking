"""
Application configuration.
Reads settings from environment variables with sensible defaults.
"""

import os
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()  # reads .env file if present


# ── Data directory ────────────────────────────────────────────
# Points to the folder containing Parquet files produced by csv_to_parquet.py
DATA_DIR = Path(os.getenv(
    "CLF_DATA_DIR",
    str(Path(__file__).resolve().parent.parent.parent / "data_prepare" / "data"),
))

# ── Parquet file paths ────────────────────────────────────────
PARQUET = {
    "dim_cell_lines":       DATA_DIR / "dim_cell_lines.parquet",
    "dim_genes":            DATA_DIR / "dim_genes.parquet",
    "fact_expression_depmap": DATA_DIR / "fact_expression_depmap.parquet",
    "fact_expression_hpa":  DATA_DIR / "fact_expression_hpa.parquet",
    "fact_expression_geo":  DATA_DIR / "fact_expression_geo.parquet",
    "dim_geo_samples":      DATA_DIR / "dim_geo_samples.parquet",   # GSM -> GSE, built by data_prepare/build_geo_study_map.py
    "fact_proteomics":      DATA_DIR / "fact_proteomics.parquet",
    "fact_mutations":       DATA_DIR / "fact_mutations.parquet",
    "fact_fusions":         DATA_DIR / "fact_fusions.parquet",
    "fact_metabolomics":    DATA_DIR / "metabolomics_clean.parquet",
    "fact_mirna":           DATA_DIR / "mirna_clean.parquet",
    "fact_signatures":      DATA_DIR / "signatures_clean.parquet",
    "dim_cell_line_parents": DATA_DIR / "dim_cell_line_parents.parquet",
}

# ── Ranking defaults ─────────────────────────────────────────
DEFAULT_W_RNA = 0.7
DEFAULT_W_PROTEIN = 0.3
DEFAULT_TOP_N = 20
RRF_K = 60

# Echoed in every /api/rank response so exported results can be traced back
# to the scoring method that produced them. Bump when scoring changes.
SCORING_VERSION = "3.0.0.dev3"

# ── Low-expression mode (AZ priority 5) ──────────────────────
# Each cell line gets an expression category from its raw TPM. The category is
# shown with every result and decides the order of "low" searches.
# Thresholds were checked on real data, see
# docs/benchmarks/2026-11-08_low_expression_thresholds.md
EXPR_DETECTION_TPM = {"depmap": 0.1, "hpa": 0.1}  # below = "below detection" (HPA reports in 0.1 steps, so this means 0)
EXPR_OFF_TPM = 1.0        # below = "negligible": the gene counts as off
EXPR_MODERATE_TPM = 3.0   # above = "moderate"; from 1 to 3 TPM = "low measurable"

# What "low" means in a search:
#   off        = button "TPM below 1": off cell lines first, tied on RNA (default)
#   measurable = button "TPM greater than 1": 1-3 TPM first, off cell lines last
#   score      = old behaviour: sort by score only (categories still shown)
LOW_EXPRESSION_MODES = ("off", "measurable", "score")
DEFAULT_LOW_EXPRESSION_MODE = "off"

# ── Off cutoff: how a gene is judged "off" in a cell line ─────
# Below EXPR_OFF_TPM (1 TPM) a gene is ALWAYS off. The other methods can also
# call a line off at 1 TPM or more when it is low compared with that gene's
# usual level (e.g. VIM 2.5 TPM in epithelial MCF7 vs ~1,000 in mesenchymal
# lines). The level of a cell line = its highest source TPM (as for the
# categories). Comparison on real data: docs/benchmarks/2026-11-08_low_expression_thresholds.md
#   fixed                below 1 TPM only (default, what the AZ team asked for)
#   percent_of_reference below OFF_REF_FRACTION x the OFF_REF_QUANTILE of lines >= 1 TPM
#   mixture              2-group split of log2(TPM + 1); the low group is off
#   distance_from_on     more than OFF_DISTANCE_SD below the expressing group
#   strict_mixture       mixture, only for genes with a clear on/off split
#   custom               below a cutoff the user chooses (>= 1 TPM)
OFF_METHODS = ("fixed", "percent_of_reference", "mixture", "distance_from_on", "strict_mixture", "custom")
DEFAULT_OFF_METHOD = "fixed"
OFF_METHOD_LABELS = {
    "fixed": "Fixed 1 TPM",
    "percent_of_reference": "% of reference",
    "mixture": "Mixture model",
    "distance_from_on": "Distance from on-group",
    "strict_mixture": "Strict mixture",
    "custom": "Custom cutoff",
}
OFF_REF_QUANTILE = 0.90        # percent_of_reference: reference = 90th percentile of lines >= 1 TPM ...
OFF_REF_FRACTION = 0.10        # ... and off below 10 % of it
OFF_MIN_LOG2_GAP = 2.0         # mixture / distance: the two groups must differ >= 4-fold ...
OFF_MIN_GROUP_SHARE = 0.05     # ... and each hold >= 5 % of the cell lines
OFF_MIXTURE_POSTERIOR = 0.5    # mixture: off if more likely in the low group
OFF_DISTANCE_SD = 3.0          # distance_from_on: off if > 3 SD below the expressing group
OFF_STRICT_MIN_D = 2.0         # strict_mixture: Ashman's D (separation of the 2 groups) above this ...
OFF_STRICT_MIN_BIC_GAIN = 10.0 # ... 2 groups fit clearly better than 1 (BIC) ...
OFF_STRICT_MIN_GROUP_SHARE = 0.03
OFF_STRICT_POSTERIOR = 0.9     # ... and off only if >= 90 % sure it is in the low group
OFF_MIN_LINES = 20             # fewer cell lines with data: data-driven methods fall back to fixed

# ── Server ───────────────────────────────────────────────────
API_HOST = os.getenv("CLF_HOST", "0.0.0.0")
API_PORT = int(os.getenv("CLF_PORT", "8000"))
CORS_ORIGINS = os.getenv("CLF_CORS_ORIGINS", "*").split(",")
