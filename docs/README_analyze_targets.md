# analyze_targets.py

Reference for `scripts/analyze_targets.py`.

This script computes per-target classification summaries from `master_alignment.csv` files in `outputs/`, renders NH-based figures, and orchestrates related analyses (CA-inclusive, same-author subset, and 1D single-atom summaries).

## What it analyzes

- **Primary NH analysis** from each target's `master_alignment.csv`
  - per-target F1 and MCC
  - histogram/stacked histogram by CA distance
  - confusion-matrix stacked histogram
- **Orchestrated analyses**:
  - `analyze_targets_ca` (CA-inclusive CSP summary; requires target selection CSV)
  - `analyze_targets_same_author`
  - `analyze_targets_single_atom_shifts`

## Usage

```bash
python scripts/analyze_targets.py [options]
```

## CLI options

- `--outputs-dir` (default: `outputs`)
  - Root directory containing target folders.
- `--output-image` (default: `outputs/f1_heatmap.png`)
  - NH heatmap output path.
- `--summary-csv` (optional)
  - Write NH per-target summary table.
- `--histogram-image` (default: `outputs/significant_ca_distance_hist.png`)
- `--stacked-histogram-image` (default: `outputs/significant_ca_distance_stacked_hist.png`)
- `--confusion-matrix-stacked-histogram-image` (default: `outputs/confusion_matrix_stacked_histogram.png`)
- `--summary-dir` (default: `outputs/summary_statistics`)
  - Canonical location for orchestrated summary artifacts.

### Target selection controls

- `--targets-csv <path>`
  - CSV with a `holo_pdb` column; filters analyzed targets.
- `--targets <id1,id2,...>`
  - Comma-separated target IDs.
  - If the same set already exists as `outputs/targets_<NAME>.csv`, that selection name is reused.
  - Otherwise prompts for a selection name unless `--selection-name` is provided.
- `--targets-from <NAME>`
  - Reuse an existing selection from `outputs/targets_<NAME>.csv`.
- `--selection-name <NAME>`
  - Name to persist a new `--targets` selection.

### Plot formatting controls

- `--no-plot-titles`
- `--axis-fontsize <float>`
- `--tick-fontsize <float>`
- `--legend-fontsize <float>`
- `--legend-title-fontsize <float>`
- `--compact-legend-labels`

## Output behavior

- Writes the requested top-level NH outputs (for example `--output-image`, `--summary-csv`).
- Always writes canonical NH summary artifacts into `--summary-dir`:
  - `f1_heatmap_nh.png`
  - `significant_ca_distance_hist_nh.png`
  - `significant_ca_distance_stacked_hist_nh.png`
  - `confusion_matrix_stacked_histogram_nh.png`
  - `f1_summary_nh.csv`
- Then runs orchestrated analyses and writes their outputs into the same summary directory.
- With named target selections (`--targets`/`--targets-from`), it also writes selection-specific files in `outputs/` (for example `f1_heatmap_<NAME>.png`, `f1_summary_<NAME>.csv`, `targets_<NAME>.csv`).

## Examples

```bash
# Analyze all targets with defaults
python scripts/analyze_targets.py

# Analyze a named subset from an existing selection file
python scripts/analyze_targets.py --targets-from TFIIH

# Analyze explicit targets and save under a new selection name
python scripts/analyze_targets.py --targets 1cf4,1d5g,1l8c --selection-name small_set

# Use a curated target CSV and custom summary directory
python scripts/analyze_targets.py --targets-csv data/targets_example.csv --summary-dir outputs/summary_statistics_example
```
