"""
CellLineFinder — build the GEO sample -> study map
===================================================
fact_expression_geo has one GSM_ID per sample but no study (GSE). The study of
each sample lives in the raw 10_GEOInfo.txt. This script turns it into a small
parquet (gsm_id, gse_id) that the backend uses to compute GEO support
*within* each study (GEO studies are never pooled).

Usage (from the repo root):
    python data_prepare/build_geo_study_map.py <path to 10_GEOInfo.txt>

Example (Git Bash):
    python data_prepare/build_geo_study_map.py ~/OneDrive/Desktop/data/nomenclature/10_GEOInfo.txt

Output: data_prepare/data/dim_geo_samples.parquet  (gsm_id, gse_id)
Samples without a GSE_ID are left out (they cannot be assigned to a study).
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "data" / "dim_geo_samples.parquet"
FACT_GEO = HERE / "data" / "fact_expression_geo.parquet"

GSM_RE = re.compile(r"^GSM\d+$")
GSE_RE = re.compile(r"^GSE\d+$")


def read_table(path: Path) -> pd.DataFrame:
    """Read a delimited text file, trying tab first, then comma / semicolon / whitespace."""
    last_err = None
    for enc in ("utf-8-sig", "latin-1"):
        for sep in ("\t", ",", ";", r"\s+"):
            try:
                df = pd.read_csv(path, sep=sep, dtype=str, encoding=enc, engine="python")
            except Exception as err:  # wrong delimiter / encoding: try the next one
                last_err = err
                continue
            if df.shape[1] >= 2:
                return df
    raise SystemExit(f"Could not read {path} as a table: {last_err}")


def find_column(df: pd.DataFrame, pattern: re.Pattern, name_hint: str) -> str | None:
    """Pick the column whose values look like the pattern (prefer one named like name_hint)."""
    best, best_share = None, 0.0
    for col in df.columns:
        vals = df[col].dropna().astype(str).str.strip()
        if len(vals) == 0:
            continue
        share = vals.map(lambda v: bool(pattern.match(v))).mean()
        named = name_hint in str(col).lower()
        if share > best_share or (share == best_share and named):
            best, best_share = col, share
    return best if best_share >= 0.5 else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("geo_info", type=Path, help="path to 10_GEOInfo.txt")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    src = args.geo_info.expanduser()
    if not src.exists():
        raise SystemExit(f"File not found: {src}")

    df = read_table(src)
    print(f"Read {len(df):,} rows, columns: {list(df.columns)}")

    gsm_col = find_column(df, GSM_RE, "gsm")
    gse_col = find_column(df, GSE_RE, "gse")
    if gsm_col is None or gse_col is None:
        print("\nCould not detect the GSM / GSE columns. First rows of the file:")
        print(df.head(5).to_string())
        raise SystemExit("Send me this output and I will adjust the script.")
    print(f"Using GSM column = {gsm_col!r}, GSE column = {gse_col!r}")

    out = pd.DataFrame({
        "gsm_id": df[gsm_col].fillna("").astype(str).str.strip(),
        "gse_id": df[gse_col].fillna("").astype(str).str.strip(),   # blank GSE_ID stays blank and is dropped below
    })
    n_all = len(out)
    out = out[out["gsm_id"].map(lambda v: bool(GSM_RE.match(v))) & out["gse_id"].map(lambda v: bool(GSE_RE.match(v)))]
    out = out.drop_duplicates("gsm_id")
    print(f"{len(out):,} of {n_all:,} samples have a GSE_ID  ->  {out['gse_id'].nunique()} studies")
    sizes = out.groupby("gse_id").size().sort_values(ascending=False)
    print("Samples per study (top 20):")
    print(sizes.head(20).to_string())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.output, index=False)
    print(f"\nWrote {args.output}")

    # Optional sanity check against the expression table
    if FACT_GEO.exists():
        try:
            import duckdb

            gsm_in_fact = duckdb.execute(
                f"SELECT DISTINCT CAST(GSM_ID AS VARCHAR) AS g FROM read_parquet('{FACT_GEO.as_posix()}')"
            ).df()["g"]
            covered = gsm_in_fact.isin(out["gsm_id"]).sum()
            print(f"{covered:,} of {len(gsm_in_fact):,} GSM_IDs in fact_expression_geo have a study")
        except Exception as err:  # the check is optional
            print(f"(skipped coverage check: {err})")


if __name__ == "__main__":
    sys.exit(main())
