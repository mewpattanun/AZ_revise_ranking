"""
Shared test setup.

Unit tests use synthetic data via monkeypatching data_service, so they need no
Parquet files. Tests marked `real_data` run against the Parquet files in
CLF_DATA_DIR and are skipped when those files are not present.

Run from the repo root:  backend_api/venv/Scripts/python -m pytest backend_api/tests
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.config import PARQUET  # noqa: E402

HAS_REAL_DATA = all(
    PARQUET[name].exists()
    for name in (
        "dim_genes", "dim_cell_lines", "fact_expression_depmap", "fact_expression_hpa",
        "fact_expression_geo", "dim_geo_samples", "fact_proteomics",
    )
)

real_data = pytest.mark.skipif(not HAS_REAL_DATA, reason="Parquet data not found in CLF_DATA_DIR")


@pytest.fixture
def fake_data(monkeypatch):
    """
    Patch data_service with a small synthetic dataset.

    Returns a dict the test can edit before calling run_ranking:
      expression[table]          -> {ach_id: tpm}   (used for every gene)
      by_gene[ensembl_id][table] -> {ach_id: tpm}   (overrides `expression` for one gene)
      protein                    -> {ach_id: intensity}
      geo_samples                -> list of (ach_id, gsm_id, gse_id, tpm)
      calls                      -> list of tables passed to get_expression_for_gene
    """
    from app.services import data_service

    data = {
        "expression": {
            "fact_expression_depmap": {"A": 100.0, "B": 50.0, "C": 10.0, "D": 1.0},
            "fact_expression_hpa":    {"A": 90.0,  "B": 60.0, "C": 5.0,  "D": 2.0},
            # GEO strongly prefers D: must have no effect on the ranking.
            "fact_expression_geo":    {"A": 1.0,   "B": 2.0,  "C": 3.0,  "D": 9999.0},
        },
        "by_gene": {},
        "protein": {},
        "geo_samples": [],
        "calls": [],
    }

    def get_expression_for_gene(ensembl_id, source_table):
        data["calls"].append(source_table)
        tables = data["by_gene"].get(ensembl_id, data["expression"])
        values = tables.get(source_table, {})
        return pd.DataFrame({"ach_id": list(values), "tpm": list(values.values())})

    def get_protein_for_gene(ensembl_id):
        values = data["protein"]
        return pd.DataFrame({"ach_id": list(values), "protein_intensity": list(values.values())})

    def geo_samples_for_gene(ensembl_id):
        rows = data["geo_samples"]
        return pd.DataFrame(
            [(a, gsm, tpm) for a, gsm, _gse, tpm in rows], columns=["ach_id", "gsm_id", "tpm"]
        )

    def geo_study_map():
        rows = data["geo_samples"]
        return pd.DataFrame(
            [(gsm, gse) for _a, gsm, gse, _tpm in rows], columns=["gsm_id", "gse_id"]
        ).drop_duplicates("gsm_id")

    monkeypatch.setattr(data_service, "get_expression_for_gene", get_expression_for_gene)
    monkeypatch.setattr(data_service, "get_protein_for_gene", get_protein_for_gene)
    monkeypatch.setattr(data_service, "_geo_samples_for_gene", geo_samples_for_gene)
    monkeypatch.setattr(data_service, "_geo_study_map", geo_study_map)
    monkeypatch.setattr(data_service, "_GEO_Z_CACHE", {})
    monkeypatch.setattr(data_service, "get_cell_lines_metadata", lambda ids: {i: {"cell_line_name": i} for i in ids})
    return data
