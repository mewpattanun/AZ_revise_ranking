# CellLineFinder benchmark — v2.1 (frozen) vs v3 (3.0.0.dev3) — 2026-11-22

## 1. Change (v2.1 → v3, top 10)

| Query | Top-10 overlap | Spearman ρ (top-50 union) | Biggest movers (top-20 union) |
|---|---|---|---|
| ERBB2 high | 5/10 | 0.655 | COLO678 #14→#39; SNU119 #17→#42; MDAMB453 #13→#34 |
| ESR1 high | 8/10 | 0.876 | SKOV3 #44→#18; PEO1 #35→#16; UACC812 #33→#17 |
| EGFR high | 4/10 | 0.652 | OSRC2 #20→#46; CAL27 #8→#29; HSC4 #11→#31 |
| VIM low | 2/10 | -0.159 | EVSAT #574→#1; SUM185PE #513→#11; NCIH1770 #390→#9 |
| MUC1 high | 3/10 | 0.295 | NCIH2228 #16→#75; NCIH2126 #15→#57; DU145 #13→#52 |
| ERBB2 + ESR1 high | 8/10 | 0.541 | COLO205 #15→#89; NCIH1781 #37→#19; SKOV3 #16→#4 |

## 2. Positive controls in top 10 (rank shown as #rank/total)

| Query | Controls | v2.1 hits | v3 hits | Control ranks v2.1→v3 |
|---|---|---|---|---|
| ERBB2 high | HER2-amplified | 2/4 | 2/4 | SKBR3 #30/1581→#16/1546, BT474 #32/1581→#20/1546, HCC1954 #2/1581→#9/1546, NCIN87 #1/1581→#2/1546 |
| ESR1 high | ER-positive luminal breast | 1/3 | 3/3 | MCF7 #4/1581→#6/1544, T47D #12/1581→#7/1544, ZR751 #13/1581→#10/1544 |
| EGFR high | EGFR-amplified | 2/3 | 3/3 | A431 #19/1581→#5/1546, HCC827 #2/1581→#8/1546, MDAMB468 #5/1581→#9/1546 |
| VIM low | luminal epithelial, vimentin-negative | 0/3 | 0/3 | MCF7 #622/1581→#148/1546, T47D #627/1581→#266/1546, ZR751 #652/1581→#424/1546 |
| ERBB2 + ESR1 high | ER-positive and HER2-amplified | 1/1 | 1/1 | BT474 #9/1581→#6/1546 |

## 3. GEO agreement (share of top 10 that per-study GEO supports)

| Query | v3 top-K | Baseline (all ranked) | Top-K with GEO | v2.1 top-K (partly circular) |
|---|---|---|---|---|
| ERBB2 high | 100% | 17% | 7/10 | 100% |
| ESR1 high | 100% | 17% | 8/10 | 100% |
| EGFR high | 100% | 13% | 3/10 | 100% |
| VIM low | 100% | 25% | 2/10 | — |
| MUC1 high | 100% | 25% | 5/10 | 100% |
| ERBB2 + ESR1 high | 88% | 5% | 8/10 | 88% |

GEO agreement for v2.1 is partly circular: v2.1 used pooled GEO in its score.
