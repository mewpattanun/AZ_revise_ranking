# Low-expression mode: thresholds and validation (3.0.0.dev3), 2026-11-08

This covers AZ priority 5: "Low expression may actually be assay noise; classify first, then rank."

## 1. Rule used

Each cell line gets an expression category from its raw TPM. The category is the highest band that any source puts the cell line in, so a line only counts as off when DepMap and HPA both agree.

| Category | Rule | Set in `config.py` |
|---|---|---|
| Missing | no RNA value in any selected source | – |
| Below detection | every source below its detection limit (DepMap < 0.1, HPA < 0.1; HPA reports in 0.1 steps, so this means 0) | `EXPR_DETECTION_TPM` |
| Negligible | every source < 1 TPM | `EXPR_OFF_TPM = 1.0` |
| Low measurable | highest source between 1 and 3 TPM | `EXPR_MODERATE_TPM = 3.0` |
| Moderate | highest source > 3 TPM | |

"Off" means Below detection or Negligible. A line is flagged ⚠ when one source says off (< 1 TPM) and another says expressed (≥ 1 TPM).

**Modes** (`low_expression_mode`, used only for LOW genes):

- `off` (default, button "TPM below 1"):
  - Off lines come first and tie on RNA. Inside that group, protein (where measured) orders them, then ACH ID.
  - Next come Low measurable lines, then Moderate lines.
  - Missing lines come last.
- `measurable` (button "TPM greater than 1"):
  - Low measurable lines where both sources agree come first, then Low measurable lines where the sources disagree.
  - Next come Moderate lines, then off lines.
  - Missing lines come last.
- `score`: the old behaviour (ascending score only). Categories are still shown.

In a multi-gene search the mode sets the RNA desirability of off lines for a LOW gene:
- `off` gives them the best value (1).
- `measurable` gives them 0, so the geometric mean treats them as a failed gene.

## 2. Why 1 TPM: the evidence

Scripts used: `low_dist.py`, `thresholds.py`, `negatives.py` and `low_bench.py`. They query the Parquet files read-only and are not committed.

### 2a. Do DepMap and HPA agree that a gene is off?

All cell line × gene pairs present in both sources:

| Cutoff | Off in DepMap | Off in HPA | Sources disagree, among pairs called off |
|---|---|---|---|
| 0.1 TPM | 28.9% | 26.9% | **21.9%** |
| 0.5 TPM | 37.4% | 37.5% | 11.2% |
| **1 TPM** | **41.1%** | **41.8%** | **9.1%** |
| 2 TPM | 45.2% | 46.5% | 8.4% |

Below 0.1 TPM, more than 1 in 5 off calls is disagreement or rounding. HPA has no values between 0 and 0.1. From 1 TPM upwards, agreement is stable.

### 2b. Is the protein detected? (Gygi mass-spec, 375 cell lines, about 12k genes)

| DepMap RNA band | 0 | 0–0.1 | 0.1–0.25 | 0.25–0.5 | 0.5–1 | 1–2 | 2–3 | 3–5 | 5–10 | > 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Protein detected | 26% | 35% | 43% | 46% | 48% | 51% | 56% | 62% | 71% | 88% |

Detection rises smoothly and there is no clear jump, so the data does not pick a cutoff for us. 1 TPM is the common "expressed" convention (HPA also uses nTPM ≥ 1).

### 2c. Known negative cell lines

| Gene | Cell line (why negative) | DepMap TPM | HPA TPM | Category | Off? |
|---|---|---|---|---|---|
| ESR1 | MDA-MB-231 (ER−) | 0.09 | 0.0 | Below detection | ✓ |
| ESR1 | BT549 (ER−) | 0.01 | 0.0 | Below detection | ✓ |
| ESR1 | SKBR3 (ER−) | 0.25 | 0.3 | Negligible | ✓ |
| ESR1 | HCC1954 (ER−) | 3.24 | 2.9 | Moderate | ✗ |
| PTEN | PC3 (deleted) | 0.88 | 0.2 | Negligible | ✓ |
| CDKN2A | A549 (deleted) | 0.11 | 0.1 | Negligible | ✓ |
| CDKN2A | U87MG (deleted) | 0.10 | 0.3 | Negligible | ✓ |
| CDKN2A | MCF7 (deleted) | 1.45 | 0.8 | Low measurable ⚠ | ✗ |
| CDH1 | MDA-MB-231 (lost) | 0.53 | 0.3 | Negligible | ✓ |
| CDH1 | HS578T (lost) | 0.83 | 0.4 | Negligible | ✓ |
| CDH1 | BT549 (lost) | 1.17 | 11.5 | Moderate | ✗ |
| VIM | MCF7 (epithelial) | 2.19 | 2.5 | Low measurable | ✗ |
| VIM | T47D (epithelial) | 2.34 | 4.4 | Moderate | ✗ |

Known positive lines (for example MCF7 and T47D for ESR1, HeLa for CDKN2A) are all above 20 TPM, so no cutoff calls them off.

Negatives caught when both sources must agree:

| Cutoff | Negatives caught |
|---|---|
| 0.1 TPM | 2/13 |
| 0.5 TPM | 5/13 |
| **1 TPM** | **8/13** |
| 2 TPM | 9/13 |

### 2d. Size of the off group for a LOW search (1 TPM)

| Query | Off cell lines | Of all ranked |
|---|---|---|
| ESR1 low | 1370 | 1544 |
| CDH1 low | 681 | 1546 |
| CDKN2A low | 423 | 1546 |
| VIM low | 58 | 1546 |
| PTEN low | 24 | 1546 |

Inside the off group, ranks are ties. The rank number inside the group only reflects protein (where measured) and the ACH-ID tie rule, so treat the whole group as equally "off".

### 2e. Multi-gene searches: firm push-down (added 2026-11-24)

After the RNA normalisation switch, z-scores use only expressed cell lines (>= 1 TPM) as the reference. For ESR1, the median of the expressed lines is 2.7 TPM. As a result, a line with 3 TPM looked "average" and lost very little desirability as a LOW gene (HCC1954: 0.78, compared with 1.00 for off lines). In ERBB2 HIGH + ESR1 LOW, its very high ERBB2 kept it at #17.

The fix: multi-gene searches now use the same "classify first" tiers as single-gene searches. Each line's tier is its worst tier over the LOW genes, and the geometric mean orders lines within a tier.

| Search | Line (LOW-gene TPM) | Before | After |
|---|---|---|---|
| ERBB2 HIGH + ESR1 LOW | HCC1954 (3.1) | #17 | #1465 |
| | EFM192A (6.0) | #29 | #1466 |
| | UACC893 (2.3) | #15 | #1371 |
| TAL1 HIGH + TAL2 LOW | F36P (6.0) | n/a | #1508 |
| | K562 (1.9) | n/a | #1426 |

1,370 lines are off for ESR1 and 1,424 for TAL2, so expressing lines now come straight after the off group. ERBB2 HIGH + CD3E LOW is unchanged, because no top ERBB2 line expresses CD3E. Mode `score` keeps the old, milder behaviour for comparison.

### 2f. Off cutoff methods: gene-relative off (added 2026-11-25)

A fixed 1 TPM cutoff cannot call a gene off when its "off" level is above 1 TPM. Examples: VIM is 2–4 TPM in epithelial MCF7/T47D against ~1,000 TPM in mesenchymal lines, and CDKN2A is 1.45 TPM in MCF7 (deleted). The API option `off_method` (sidebar: **Off cutoff**) chooses how "off" is decided. Below 1 TPM is always off. A line at 1 TPM or more that a method calls off gets the category `low_relative` and the flag "Low for GENE". With "TPM below 1", the lines below 1 TPM come first and tie, as before. The `low_relative` lines come next, ordered by their score: their RNA is not tied, because 2 vs 20 TPM is a real difference, not noise, and protein counts where it was measured. With "TPM greater than 1", they come just before the lines below 1 TPM at the end.

| Method | Rule (settings in `config.py`) |
|---|---|
| `fixed` (default) | below 1 TPM only |
| `percent_of_reference` | below 10% of the 90th percentile of lines ≥ 1 TPM |
| `mixture` | low group of a 2-group Gaussian split of log2(TPM + 1) (groups ≥ 4-fold apart, each ≥ 5% of lines) |
| `distance_from_on` | more than 3 SD below the high (expressing) group of that split |
| `strict_mixture` | mixture, only for a clear split (Ashman's D > 2, BIC gain > 10), and ≥ 90% likely in the low group |
| `custom` | below `off_cutoff_tpm` (≥ 1) |

Known cell lines: 13 known negatives and 16 lines that must not be called off. The 16 include normal, non-amplified ERBB2, to catch over-calling. Results are from the implemented code.

| Method | Negatives caught | Wrongly called off |
|---|---|---|
| fixed | 8/13 | 0/16 |
| percent_of_reference | **12/13** | 2/16 (EGFR in MCF7, T47D) |
| mixture | 10/13 | 0/16 |
| distance_from_on | 10/13 | 0/16 |
| strict_mixture | 10/13 | **0/16** |

Extra lines called off at ≥ 1 TPM (out of ~1,545):

| Gene | percent_of_reference | mixture | distance_from_on | strict_mixture |
|---|---|---|---|---|
| VIM | 674 | 479 | 403 | 433 |
| CDH1 | 271 | 56 | 0 | 20 |
| CDKN2A | 354 | 12 | 0 | 0 |
| ERBB2 | 312 | 209 | 0 | 0 |
| EGFR | 208 | 0 | 0 | 0 |
| ESR1 | 62 | 0 | 0 | 0 |
| PTEN | 20 | 0 | 0 | 0 |

`strict_mixture` catches the VIM cases (MCF7, T47D) without calling normal ERBB2 lines off; `mixture` and `percent_of_reference` over-call (ERBB2 +209 / +312). No data-driven method catches HCC1954 (ESR1 3.2 TPM): ESR1-high lines are few and many lines sit between 1 and 10 TPM, so ESR1 has no clear split. A `custom` cutoff (e.g. 10 TPM) handles that case. The default stays `fixed`, as the AZ team asked; the other methods are opt-in.

## 3. Known limitations

1. **Off only relative to the gene's normal level.** Examples are HCC1954 ESR1 (~3 TPM against ~40 TPM in ER+ lines) and VIM in epithelial lines (2–4 TPM against 2,000+ TPM in mesenchymal lines). A fixed TPM cutoff cannot call these off. The `off_method` option (section 2f) now offers gene-relative rules; `strict_mixture` catches the VIM cases, and ESR1 needs a `custom` cutoff.
2. **Protein lost but RNA normal.** MDA-MB-468 has no PTEN protein (PTEN is mutated) but 19 TPM of RNA. RNA thresholds cannot catch this; use the mutation filter or protein data.
3. **Missing protein still counts as zero for low searches.** Inside the off group, RNA-only lines score 0 and come before lines with measured protein. This bias was already there before this change (the old `score` mode does the same) and belongs to AZ priority 6 (separate coverage from the score).
4. **Same cutoffs for both sources.** These are TPM values from two processing pipelines of the same data. A per-source calibration would need an independent assay.
