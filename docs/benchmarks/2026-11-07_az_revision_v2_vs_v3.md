# CellLineFinder benchmark — v2.1 (frozen) vs v3 (3.0.0.dev2) — 2026-11-07

## 1. Change (v2.1 → v3, top 10)

| Query | Top-10 overlap | Spearman ρ (top-50 union) | Biggest movers (top-20 union) |
|---|---|---|---|
| ERBB2 high | 9/10 | 0.584 | HCC1569 #19→#44; HCC38 #40→#20; ZR7530 #15→#34 |
| ESR1 high | 8/10 | 0.876 | SKOV3 #44→#18; PEO1 #35→#16; UACC812 #33→#17 |
| EGFR high | 8/10 | 0.411 | YSCCC #17→#64; A431 #19→#61; UMUC5 #16→#57 |
| VIM low | 10/10 | 0.846 | UPCISCC154 #5→#1; G415 #1→#4; LU135 #2→#3 |
| MUC1 high | 7/10 | 0.286 | UACC812 #10→#79; CAPAN2 #11→#80; CCC5 #7→#75 |
| ERBB2 + ESR1 high | 5/10 | 0.530 | COLO205 #15→#650; KYSE30 #62→#18; T47D #17→#5 |

## 2. Positive controls in top 10 (rank shown as #rank/total)

| Query | Controls | v2.1 hits | v3 hits | Control ranks v2.1→v3 |
|---|---|---|---|---|
| ERBB2 high | HER2-amplified | 2/4 | 2/4 | SKBR3 #30/1581→#50/1546, BT474 #32/1581→#58/1546, HCC1954 #2/1581→#4/1546, NCIN87 #1/1581→#1/1546 |
| ESR1 high | ER-positive luminal breast | 1/3 | 3/3 | MCF7 #4/1581→#6/1544, T47D #12/1581→#7/1544, ZR751 #13/1581→#10/1544 |
| EGFR high | EGFR-amplified | 2/3 | 2/3 | A431 #19/1581→#61/1546, HCC827 #2/1581→#2/1546, MDAMB468 #5/1581→#3/1546 |
| VIM low | luminal epithelial, vimentin-negative | 0/3 | 0/3 | MCF7 #622/1581→#276/1546, T47D #627/1581→#325/1546, ZR751 #652/1581→#490/1546 |
| ERBB2 + ESR1 high | ER-positive and HER2-amplified | 1/1 | 0/1 | BT474 #9/1581→#14/1546 |

## 3. GEO agreement (share of top 10 that per-study GEO supports)

| Query | v3 top-K | Baseline (all ranked) | Top-K with GEO | v2.1 top-K (partly circular) |
|---|---|---|---|---|
| ERBB2 high | 100% | 17% | 7/10 | 100% |
| ESR1 high | 100% | 17% | 8/10 | 100% |
| EGFR high | 100% | 13% | 8/10 | 100% |
| VIM low | — | 25% | 0/10 | — |
| MUC1 high | 100% | 25% | 6/10 | 100% |
| ERBB2 + ESR1 high | 60% | 5% | 10/10 | 88% |

GEO agreement for v2.1 is partly circular: v2.1 used pooled GEO in its score.

## Notes

v2.1 = frozen dissertation ranking (`ranking_service_legacy.py`: rank-based RRF, pooled GEO in the score). v3 = AZ revision points 1-4 (`ranking_service.py`: DepMap+HPA robust z consensus, asinh magnitude score, weighted geometric mean, GEO as per-study evidence only), plus the deterministic tie order and core-only fix from this PR. Both use default weights (RNA 0.7, protein 0.3), no filters.

- **Cell lines without protein are capped below every cell line with protein.** SKBR3 (#30→#50) and BT474 (#32→#58) for ERBB2, and A431 (#19→#61) for EGFR, have no proteomics data. Their score is "RNA only" = 0.7 × RNA score, so it can never exceed 0.7, while every top-20 line has RNA+Protein and scores above 0.9. Their RNA consensus alone ranks them #15, #19 and #6. v2.1 had the same cap, but there the scores were rank-based and closer together; with magnitude scores the gap is larger. This is the "separate coverage from score" issue (AZ priority 6) and is **not** changed here.
- **MUC1 (ρ 0.29) and EGFR (ρ 0.41) changed most.** For MUC1 the top three by score are only #20-#39 by RNA z, so protein is driving the order. Worth checking against MUC1 biology before relying on it.
- **COLO205 #15 → #650 (ERBB2 + ESR1):** COLO205 has no DepMap or HPA RNA for either gene; its v2.1 rank came only from pooled GEO. Intended effect of removing GEO from the score.
- **Ranked list shrinks (~1,581 → ~1,546):** about 35 cell lines had GEO as their only RNA evidence and no longer rank.
- **ESR1 controls improved** (1/3 → 3/3 in the top 10). **VIM low controls moved up** (MCF7 #622 → #276) but are still far from the top: many lines have near-zero VIM and tie (AZ priority 5, detection-aware low expression).
- **GEO agreement:** per-study GEO (at least 10 samples per study, |z| ≥ 1, see `data_service.classify_geo`) is stricter than the old pooled percentile, so the baseline is low (5-25 %). v3 top-10 lines with GEO data are all supported for single-gene queries; 60 % for ERBB2 + ESR1, where both genes must be supported.

Generated with `backend_api/compare_v2_v3.py --markdown`.
