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

# ── Server ───────────────────────────────────────────────────
API_HOST = os.getenv("CLF_HOST", "0.0.0.0")
API_PORT = int(os.getenv("CLF_PORT", "8000"))
CORS_ORIGINS = os.getenv("CLF_CORS_ORIGINS", "*").split(",")
