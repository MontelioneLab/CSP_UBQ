# HN CSP significance threshold columns

This document describes the boolean significance columns written to each target’s
`csp_table.csv` and merged into `master_alignment.csv`. All masks apply to the
combined HN CSP (`csp_A`) after referencing (when enabled). A residue is marked
significant when `csp_A >= threshold` for that method.

Defaults live in [`scripts/config.py`](../scripts/config.py) (`Thresholds`).
Masks are assigned in [`scripts/csp.py`](../scripts/csp.py)
(`_apply_hn_significance_masks`).

## Value set used for statistics

Unless noted, mean / SD / percentile statistics use the **final HN CSP value
list** for the target: referenced CSPs when mean-based or grid referencing is
applied, otherwise the unre-referenced recorded CSPs. Residues excluded for
large raw `|ΔN|` are not in this list.

**Population SD** (`statistics.pstdev`) is used wherever SD appears below.

## Cleaned-mean path (iterative outlier removal)

Iterative outlier removal (`compute_threshold_with_outlier_removal`):

1. Start with all recorded HN CSP values.
2. Remove values above `mean + outlier_z_score · SD` (default `outlier_z_score = 3.0`).
3. Repeat until stable or limits hit (`max_outlier_iterations`, `max_outlier_fraction`).
4. Compute **cleaned mean** and **cleaned SD** on the remaining set.

The cleaned-mean path cutoff before flooring is:

`cleaned_mean + significance_z_score · cleaned_SD`

With default `significance_z_score = 0.0`, that equals the cleaned mean. The
**primary** pipeline cutoff then applies a floor of `cutoff_05_ppm` (default
0.05 ppm):

`max(cleaned_mean + significance_z_score · cleaned_SD, cutoff_05_ppm)`

With default `significance_z_score = 0.0`, primary equals
`max(cleaned_mean, 0.05 ppm)`.

| Column | Threshold | Outlier removal |
|--------|-----------|-----------------|
| `significant` | Primary cutoff (`max(cleaned_mean + significance_z · cleaned_SD, cutoff_05_ppm)`) | Yes |
| `significant_sigma_0` | Cleaned mean (no floor) | Yes |
| `significant_1sd` / `significant_sigma_1` | Cleaned mean + 1·cleaned SD | Yes |
| `significant_2sd` / `significant_sigma_2` | Cleaned mean + 2·cleaned SD | Yes |

When `significance_z_score = 0` and cleaned mean ≥ 0.05 ppm, `significant`
matches `significant_sigma_0`. When cleaned mean is below 0.05 ppm, primary
uses the 0.05 ppm floor (same as `significant_max_05_cleaned_mean`).

`csp_z` (when present) is the z-score relative to cleaned mean/SD.

Residues with `csp_z` above `Thresholds.max_classification_csp_z` (default **100**)
are left unclassified: empty `classification` in the master CSV, omitted from
CSP classification bars, and gray apo→holo connectors on aligned HSQC overlays.

## Raw-mean path (no outlier removal)

Statistics are computed on the full final HN CSP value list with **no** iterative
outlier removal.

| Column | Threshold | Outlier removal |
|--------|-----------|-----------------|
| `significant_raw_mean` | Raw mean | No |
| `significant_raw_1sd` | Raw mean + 1·raw SD | No |
| `significant_raw_2sd` | Raw mean + 2·raw SD | No |

## Floor on cleaned mean

| Column | Threshold | Outlier removal |
|--------|-----------|-----------------|
| `significant_max_05_cleaned_mean` | `max(cutoff_05_ppm, cleaned_mean)` | Cleaned mean uses outlier removal; floor is fixed |

Default `cutoff_05_ppm = 0.05`. This is stricter than cleaned mean alone when the
cleaned mean is below 0.05 ppm, and identical to cleaned mean when cleaned mean
is already ≥ 0.05 ppm.

## Rank-based (percentile) masks

Thresholds are linear-interpolation percentiles of the final HN CSP value list
(same list as raw-mean stats). A residue is significant if
`csp_A >= percentile`.

| Column | Threshold | Config |
|--------|-----------|--------|
| `significant_top_10_percentile` | 90th percentile (top 10%) | `top_10_percentile_fraction = 0.10` |
| `significant_top_5_percentile` | 95th percentile (top 5%) | `top_5_percentile_fraction = 0.05` |

## Fixed absolute cutoffs (ppm)

| Column | Threshold | Config |
|--------|-----------|--------|
| `significant_03_ppm` | 0.03 ppm | `cutoff_03_ppm` |
| `significant_05_ppm` | 0.05 ppm | `cutoff_05_ppm` |
| `significant_10_ppm` | 0.10 ppm | `cutoff_10_ppm` |

These do not use mean/SD or outlier removal.

## Encoding in CSV outputs

In `csp_table.csv` and `master_alignment.csv`, boolean significance columns are
written as `1` / `0` when defined, or left blank when the residue has no
recorded `csp_A` (or when stats could not be computed).

## Related figure

[`scripts/create_fig_3_thresholds.py`](../scripts/create_fig_3_thresholds.py)
compares these masks in a multi-panel TP/FP and confusion-matrix distance
histogram (`figures/figure_3_thresholds.png`).

## Continuous ROC and intensity analysis (SI Fig. S17C and D)

`scripts/analyze_csp_thresholds.py` reads existing master tables without rerunning
alignment, referencing, or structural calculations. The canonical binding-site
union is `scripts.merge_csv.is_binding_site_row`, including direct atom contacts.
Positive evidence establishes a positive label; a negative label requires all
structural components to be observed. Missing structural evidence is unknown,
never an automatic negative. Eligibility also requires a finite nonnegative HN
CSP and a saved significance mask, excludes prolines, and preserves the existing
`csp_z > 100` exclusion. Eligibility is fixed throughout every sweep.

Each target's cutoff is reconstructed from **all** recorded CSPs using the pinned
cleaned-mean/floor procedure, before eligibility exclusions. Outlier removal
estimates the cutoff; it does not remove ordinary high CSPs from the ROC.
The two score scales are raw `csp_A` (ppm) and `csp_A / primary_cutoff`.
Sweeps use inclusive `score >= threshold`, group ties, and include the all-negative
and all-positive endpoints. Per-target AUC is unchanged by positive scaling;
pooled AUC can differ. Equal-target mean AUC is the arithmetic mean of valid
target AUCs; the mean ROC interpolates target TPRs on a common 1,001-point FPR grid.
Single-class targets have undefined AUC and are listed explicitly.

The relative sweep predicts significant CSP when `csp_A >= m * primary_cutoff`.
`m=1` reconstructs the primary rule, within saved CSV rounding. Saved `significant`
flags define the exact baseline and TP/FP violin membership. They are never
silently replaced by reconstructed predictions. `residue_eligibility.csv` records
boundary disagreements. ROC FPR is `FP/(FP+TN)`; distal-CSP fraction (FDR) is
`FP/(FP+TP)`. Undefined precision, FDR, and TP/FP ratios are blank, not zero.

Uncertainty uses 2,000 seeded whole-target bootstrap replicates, retaining all
residues within each resampled target. Intervals are percentile 95%; bands are
pointwise. They describe variation across target pairs and do not account for
dependence between different targets from the same receptor family.

Intensity violins use a Scott-bandwidth KDE in log10 CSP space, bounded by the
observed range and renormalized to unit area, sharing one width scale within each
TP/FP pair. Medians and IQRs are in ppm on a log axis. Groups with fewer than five
observations or zero variance show points; empty groups are labeled. Pooled
residue distributions and per-target median distributions answer different
weighting questions; the latter shows paired target observations.

The original 2026-11-04 snapshot contains 135 targets and 12,267 eligible residues:
522 rows lack recorded CSPs, 103 have unknown structural labels, and one has
`csp_z > 100`. Baseline counts are TP=1,331, FP=3,023, TN=7,031, FN=882.
`2JMF_6262` has no TP at baseline; every target still has both structural classes.
One rounded boundary observation in `6QXZ_27250` differs from the reconstructed
cutoff. The run manifest records input/source hashes, source snapshots, versions,
settings, seed, creation time with timezone, and predecessor when supplied.
Changing inputs/settings requires a new immutable run linked with `--predecessor`.

SI Fig. S17C shows the ROC analysis and S17D shows pooled significant-CSP
magnitudes, alongside the two distance panels in S17A and B. The original
four-panel threshold overview and 135-target atlas are available in the frozen run.

S17D's symmetric KL test uses 20 shared log-spaced bins and a pseudocount of
0.5 per bin. It compares the observed divergence against 99,999 random
permutations that preserve each target's TP and FP counts. The result and bin
sensitivity checks are in `outputs/2026-11-11_csp-intensity-kl_001/result.json`.
Reproduce them with `scripts/analyze_csp_intensity_kl.py`, passing the frozen
publication run as `--analysis-dir` and a new JSON path as `--output`.
