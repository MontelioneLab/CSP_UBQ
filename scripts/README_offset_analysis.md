# Offset Analysis Script

This script analyzes H_offset and N_offset values from CSP table CSV files across all targets and creates histograms showing the distribution of these values.

## Usage

```bash
python scripts/analyze_offsets.py [options]
```

### Options

- `--outputs-dir`: Path to outputs directory containing CSP table files (default: `outputs`)
- `--output-dir`: Directory to save analysis results (default: `offset_analysis`)

### Examples

```bash
# Use default settings
python scripts/analyze_offsets.py

# Specify custom directories
python scripts/analyze_offsets.py --outputs-dir /path/to/outputs --output-dir /path/to/results
```

## Output Files

The script generates the following files in the output directory:

1. **`offset_distributions_all_targets.png`**: Histograms showing the overall distribution of H_offset and N_offset values across all targets
2. **`offset_distributions_by_target.png`**: Box plots showing the distribution of offset values for each individual target
3. **`offset_correlation.png`**: Scatter plot showing the correlation between H_offset and N_offset values
4. **`offset_statistics.csv`**: Detailed statistics including mean, standard deviation, min, max, median, and quartiles for each target and offset type

## Key Findings

From the analysis of 6,736 offset values across 64 targets:

- **H_offset**: Mean = 0.0036 ppm, Std = 0.0445 ppm, Range = [-0.1062, 0.1192] ppm
- **N_offset**: Mean = -0.0119 ppm, Std = 0.3121 ppm, Range = [-1.0872, 1.2836] ppm

The H_offset values are generally small and centered around zero, while N_offset values show more variation and tend to be slightly negative on average.

## Requirements

- Python 3.6+
- numpy
- matplotlib
- Standard library modules (os, csv, glob, pathlib, argparse, typing)
