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

## 3. Known limitations

1. **Off only relative to the gene's normal level.** Examples are HCC1954 ESR1 (~3 TPM against ~40 TPM in ER+ lines) and VIM in epithelial lines (2–4 TPM against 2,000+ TPM in mesenchymal lines). A fixed TPM cutoff cannot call these off. A gene-relative rule ("below X% of the level in expressing lines") could be added later as a flag.
2. **Protein lost but RNA normal.** MDA-MB-468 has no PTEN protein (PTEN is mutated) but 19 TPM of RNA. RNA thresholds cannot catch this; use the mutation filter or protein data.
3. **Missing protein still counts as zero for low searches.** Inside the off group, RNA-only lines score 0 and come before lines with measured protein. This bias was already there before this change (the old `score` mode does the same) and belongs to AZ priority 6 (separate coverage from the score).
4. **Same cutoffs for both sources.** These are TPM values from two processing pipelines of the same data. A per-source calibration would need an independent assay.
