# Offset Analysis

This document describes `scripts/analyze_offsets.py`, which summarizes holo referencing offsets across targets.

The script pulls two data sources from `outputs/`:

- Per-target `csp_table.csv` files (`H_offset`, `N_offset`)
- Per-target grid-search summary CSV files named:
  `offset_grid_H_-0.12_0.12_0.01__N_-1.2_1.2_0.05__C_0.05.csv`

For duplicate systems, suffixed target directories (for example `2K7A_1`, `2K7A_2`) are treated as separate targets.

## Usage

```bash
python scripts/analyze_offsets.py [options]
```

## Options

- `--outputs-dir` (default: `outputs`)
  - Root directory scanned recursively for `csp_table.csv` and fixed-name grid CSV files.
- `--output-dir` (default: `offset_analysis`)
  - Destination directory for generated figures and summary CSV.

## Examples

```bash
# Defaults
python scripts/analyze_offsets.py

# Custom input/output locations
python scripts/analyze_offsets.py --outputs-dir outputs --output-dir outputs/offset_analysis
```

## Generated outputs

The script writes the following files to `--output-dir`:

- `offset_grid_heatmap.png` (when grid CSVs are found)
- `offset_distributions_all_targets.png`
- `offset_values_by_target.png`
- `offset_correlation.png`
- `offset_statistics.csv`

## Notes

- The script currently expects the fixed grid-search filename shown above.
- If no valid offsets are found in `csp_table.csv`, plotting is skipped with a console message.
