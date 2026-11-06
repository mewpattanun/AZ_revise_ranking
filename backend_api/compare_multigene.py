"""
Multi-gene validation matrix (priority 4, storyboard slide 13)
==============================================================
Which multi-gene rule retrieves candidates that satisfy the intended joint profile?

Builds 1,000 synthetic cell lines with 3 genes (per-gene desirability 0-1, higher =
better), injects the scenario profiles from the slide, and prints the rank each
profile gets (1 = best) under five methods:

  summed_rrf  legacy: sum of 1/(60 + rank) over the genes
  arithmetic  mean of the per-gene desirability
  geometric   weighted geometric mean (the new default)
  minimum     the weakest gene
  threshold   mean, but a cell line is only eligible when EVERY gene is >= 0.5
              (not eligible = ranked last)

A missing gene (NA) is skipped by every method.

Usage (from backend_api, no database needed):
    python compare_multigene.py
"""

import numpy as np
import pandas as pd

from app.config import RRF_K
from app.services import ranking_service as rk

N = 1000
GENES = ["G1", "G2", "G3"]
THRESHOLD = 0.5

rng = np.random.default_rng(0)
population = pd.DataFrame(rng.uniform(0, 1, (N, 3)), columns=GENES)
population.index = [f"pop{i:04d}" for i in range(N)]


def at_rank(gene: str, rank: int) -> float:
    """Desirability a cell line needs to sit at `rank` for this gene in the population."""
    return float(np.sort(population[gene].to_numpy())[::-1][rank - 1])


NA = np.nan
# For a LOW gene the slide's raw value is inverted first: desirability = 1 - raw
profiles = {
    "Balanced high         0.85 | 0.82 | 0.80": (0.85, 0.82, 0.80),
    "One exceptional       1.00 | 0.40 | 0.38": (1.00, 0.40, 0.38),
    "One failed            0.95 | 0.90 | 0.02": (0.95, 0.90, 0.02),
    "Mixed dir. satisfied  HIGH 0.90 | HIGH 0.85 | LOW 0.10": (0.90, 0.85, 1 - 0.10),
    "Mixed dir. violated   HIGH 0.90 | HIGH 0.85 | LOW 0.90": (0.90, 0.85, 1 - 0.90),
    "Missing gene          0.85 | NA | 0.82": (0.85, NA, 0.82),
    "Rank trade-off X      #1 + #300": (at_rank("G1", 1), at_rank("G2", 300), NA),
    "Rank trade-off Y      #15 + #18": (at_rank("G1", 15), at_rank("G2", 18), NA),
}
injected = pd.DataFrame(list(profiles.values()), columns=GENES, index=list(profiles))
D = pd.concat([population, injected])


def method_scores(D: pd.DataFrame) -> pd.DataFrame:
    ranks = {g: D[g].dropna().rank(ascending=False, method="first") for g in GENES}
    summed = rk._rrf(ranks, average=False) / (len(GENES) / (RRF_K + 1))
    eligible = (D.fillna(1.0) >= THRESHOLD).all(axis=1)
    return pd.DataFrame({
        "summed_rrf": summed.reindex(D.index),
        "arithmetic": D.mean(axis=1, skipna=True),
        "geometric": rk.weighted_geometric_mean(D),
        "minimum": D.min(axis=1, skipna=True),
        "threshold": D.mean(axis=1, skipna=True).where(eligible, -1.0),
    })


scores = method_scores(D)
ranks = scores.rank(ascending=False, method="min").astype(int)

pd.set_option("display.width", 200)
print(f"Rank of each scenario among {len(D)} cell lines (1 = best), and the score in brackets\n")
table = pd.DataFrame({
    m: [f"{ranks.loc[p, m]:>4d} ({scores.loc[p, m]:.3f})" if scores.loc[p, m] >= 0 else f"{ranks.loc[p, m]:>4d} (excl.)"
        for p in profiles]
    for m in scores.columns
}, index=list(profiles))
print(table.to_string())

print("""
How to read it
  One failed (0.95 | 0.90 | 0.02): the slide's problem case. Summed RRF keeps it near the top
  because two genes are excellent; geometric mean pushes it far down.
  Mixed direction: the LOW gene is inverted first, so only the 'satisfied' profile scores well.
  Missing gene: skipped, the other genes decide.
  Rank trade-off: X is #1 on one gene and #300 on the other, Y is #15 and #18.
""")
