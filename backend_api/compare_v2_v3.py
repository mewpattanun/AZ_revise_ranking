"""
Benchmark the current ranking (v3: robust z RNA consensus + geometric mean)
against the frozen v2.1 ranking (rank-based RRF with pooled GEO in the score).

Three checks per benchmark query:
  1. Change      — top-K overlap, Spearman rank correlation, biggest movers.
                   Says how much changed, not whether it improved.
  2. Controls    — where well-established positive-control cell lines rank
                   (e.g. HER2-amplified lines for ERBB2 high). Higher = better.
  3. GEO agree   — share of the top K that per-study GEO supports, against the
                   share for all ranked cell lines with GEO data (baseline).
                   A cell line counts as supported when its GEO status is
                   "supportive" for every gene GEO covers ("mixed" does not
                   count). v2.1 used pooled GEO in its score, so its figure is
                   partly circular.

Run from inside the backend_api folder (CLF_DATA_DIR defaults to
data_prepare/data):
    python compare_v2_v3.py
    python compare_v2_v3.py --top-k 20 --markdown ../docs/benchmarks/report.md
    python compare_v2_v3.py --query ERBB2:high,ESR1:high   (ad-hoc, no controls)
"""

import argparse
import sys
from dataclasses import dataclass, field
from datetime import date

from scipy.stats import spearmanr

from app.config import SCORING_VERSION
from app.services import data_service
from app.services import ranking_service as v3
from app.services import ranking_service_legacy as v2


@dataclass
class Benchmark:
    label: str
    genes: list[tuple[str, str]]          # [(hugo, direction)]
    controls: list[str] = field(default_factory=list)   # cell_line_name values
    rationale: str = ""


# Positive controls are textbook cell-line models for each target. Edit or
# extend this list as the team agrees on more controls.
BENCHMARKS = [
    Benchmark("ERBB2 high", [("ERBB2", "high")],
              ["SKBR3", "BT474", "HCC1954", "NCIN87"], "HER2-amplified"),
    Benchmark("ESR1 high", [("ESR1", "high")],
              ["MCF7", "T47D", "ZR751"], "ER-positive luminal breast"),
    Benchmark("EGFR high", [("EGFR", "high")],
              ["A431", "HCC827", "MDAMB468"], "EGFR-amplified"),
    Benchmark("VIM low", [("VIM", "low")],
              ["MCF7", "T47D", "ZR751"], "luminal epithelial, vimentin-negative"),
    Benchmark("MUC1 high", [("MUC1", "high")]),
    Benchmark("ERBB2 + ESR1 high", [("ERBB2", "high"), ("ESR1", "high")],
              ["BT474"], "ER-positive and HER2-amplified"),
]


def build_genes(pairs: list[tuple[str, str]]) -> list[dict]:
    """Resolve (hugo, direction) pairs to the gene dicts the pipeline expects."""
    genes = []
    for hugo, direction in pairs:
        resolved = data_service.resolve_gene(hugo)
        if resolved is None:
            print(f"Gene '{hugo}' not found in dim_genes.")
            sys.exit(1)
        genes.append({
            "hugo": resolved["hugo_symbol"],
            "ensembl_id": resolved["ensembl_id"],
            "direction": direction,
            "weight": 1.0,
        })
    return genes


def full_ranking(module, genes: list[dict]) -> list[str]:
    """Full ranked list of ach_ids (no top-N cut, no filters, module defaults)."""
    return [r["ach_id"] for r in module.run_ranking(genes, top_n=100000)]


def change_metrics(old: list[str], new: list[str], top_k: int) -> dict:
    old_rank = {a: i + 1 for i, a in enumerate(old)}
    new_rank = {a: i + 1 for i, a in enumerate(new)}
    union = [a for a in dict.fromkeys(old[:50] + new[:50]) if a in old_rank and a in new_rank]
    rho = spearmanr([old_rank[a] for a in union], [new_rank[a] for a in union]).statistic if len(union) > 2 else float("nan")

    movers_pool = [a for a in dict.fromkeys(old[:20] + new[:20]) if a in old_rank and a in new_rank]
    movers = sorted(movers_pool, key=lambda a: abs(old_rank[a] - new_rank[a]), reverse=True)[:3]
    return {
        "overlap": len(set(old[:top_k]) & set(new[:top_k])),
        "rho": rho,
        "movers": [(a, old_rank[a], new_rank[a]) for a in movers if old_rank[a] != new_rank[a]],
    }


def control_ranks(ranking: list[str], control_ids: list[str]) -> list[int | None]:
    pos = {a: i + 1 for i, a in enumerate(ranking)}
    return [pos.get(a) for a in control_ids]


def geo_agreement(ranking: list[str], genes: list[dict], top_k: int) -> tuple[float | None, float | None, int]:
    """
    (top-K supportive share, baseline supportive share, n top-K with GEO),
    using the same per-study GEO classification as the app's detail panel.
    """
    statuses = {}
    for g in genes:
        geo = data_service.get_geo_support(g["ensembl_id"], g["direction"])
        statuses[g["hugo"]] = geo[0]["geo_status"].to_dict() if geo is not None else {}

    def supported(ach):
        covered = [statuses[g["hugo"]][ach] for g in genes if ach in statuses[g["hugo"]]]
        return None if not covered else all(s == "supportive" for s in covered)

    flags = [supported(a) for a in ranking]
    top = [f for f in flags[:top_k] if f is not None]
    allf = [f for f in flags if f is not None]
    share = lambda xs: sum(xs) / len(xs) if xs else None
    return share(top), share(allf), len(top)


def fmt_pct(x):
    return "—" if x is None else f"{x:.0%}"


def fmt_rank(r, total):
    return "absent" if r is None else f"#{r}/{total}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top-k", type=int, default=10, help="Top-K used for overlap, control hits and GEO agreement")
    parser.add_argument("--query", help="Ad-hoc query instead of the benchmark suite, e.g. ERBB2:high,ESR1:high")
    parser.add_argument("--markdown", help="Also write the report as Markdown to this path")
    args = parser.parse_args()
    k = args.top_k
    # Windows consoles default to cp1252, which cannot print → or ρ.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    if args.query:
        pairs = [tuple(p.split(":")) if ":" in p else (p, "high") for p in args.query.split(",")]
        benchmarks = [Benchmark(args.query, pairs)]
    else:
        benchmarks = BENCHMARKS

    dim = data_service.query_df(
        "SELECT CAST(ach_id AS VARCHAR) AS ach_id, CAST(cell_line_name AS VARCHAR) AS name FROM dim_cell_lines"
    )
    names = dict(zip(dim["ach_id"], dim["name"].fillna(dim["ach_id"])))
    by_name = {v: a for a, v in names.items()}

    change_rows, control_rows, geo_rows = [], [], []
    for b in benchmarks:
        print(f"Running {b.label} ...", file=sys.stderr)
        genes = build_genes(b.genes)
        old, new = full_ranking(v2, genes), full_ranking(v3, genes)

        c = change_metrics(old, new, k)
        movers = "; ".join(f"{names.get(a, a)} #{o}→#{n}" for a, o, n in c["movers"]) or "—"
        change_rows.append([b.label, f"{c['overlap']}/{k}", f"{c['rho']:.3f}", movers])

        if b.controls:
            ids = [by_name.get(n) for n in b.controls]
            old_r, new_r = control_ranks(old, ids), control_ranks(new, ids)
            hits = lambda rs: sum(1 for r in rs if r is not None and r <= k)
            detail = ", ".join(
                f"{n} {fmt_rank(o, len(old))}→{fmt_rank(nw, len(new))}"
                for n, o, nw in zip(b.controls, old_r, new_r)
            )
            control_rows.append([b.label, b.rationale, f"{hits(old_r)}/{len(ids)}", f"{hits(new_r)}/{len(ids)}", detail])

        old_top, _, _ = geo_agreement(old, genes, k)
        new_top, base, n_geo = geo_agreement(new, genes, k)
        geo_rows.append([b.label, fmt_pct(new_top), fmt_pct(base), f"{n_geo}/{k}", fmt_pct(old_top)])

    sections = [
        (f"1. Change (v2.1 → v3, top {k})",
         ["Query", f"Top-{k} overlap", "Spearman ρ (top-50 union)", "Biggest movers (top-20 union)"], change_rows),
        (f"2. Positive controls in top {k} (rank shown as #rank/total)",
         ["Query", "Controls", "v2.1 hits", "v3 hits", "Control ranks v2.1→v3"], control_rows),
        (f"3. GEO agreement (share of top {k} that per-study GEO supports)",
         ["Query", "v3 top-K", "Baseline (all ranked)", "Top-K with GEO", "v2.1 top-K (partly circular)"], geo_rows),
    ]

    header = f"CellLineFinder benchmark — v2.1 (frozen) vs v3 ({SCORING_VERSION}) — {date.today().isoformat()}"
    print("\n" + header)
    for title, cols, rows in sections:
        if not rows:
            continue
        print(f"\n{title}")
        widths = [max(len(str(x)) for x in col) for col in zip(cols, *rows)]
        for row in [cols] + rows:
            print("  ".join(str(x).ljust(w) for x, w in zip(row, widths)))

    if args.markdown:
        lines = [f"# {header}", ""]
        for title, cols, rows in sections:
            if not rows:
                continue
            lines += [f"## {title}", "", "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
            lines += ["| " + " | ".join(str(x) for x in row) + " |" for row in rows]
            lines.append("")
        lines.append("GEO agreement for v2.1 is partly circular: v2.1 used pooled GEO in its score.")
        with open(args.markdown, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print(f"\nMarkdown report written to {args.markdown}")


if __name__ == "__main__":
    main()
