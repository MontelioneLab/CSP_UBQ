# Supplemental Figures and Tables Wrapper

This page documents `scripts/create_all_supplement_figures_and_tables.py`, the wrapper that runs the supplementary figure/table generators in sequence.

## What it does

The wrapper orchestrates the SI generation scripts (figures, tables, and equation renderings) using a single command. It passes consistent paths to sub-scripts and reports failures at the end (or stops early with `--stop-on-error`).

## Usage

```bash
python scripts/create_all_supplement_figures_and_tables.py \
  --outputs-dir outputs \
  --figures-dir figures \
  --csv data/CSP_UBQ.csv
```

## CLI options

- `--outputs-dir` (default: `outputs`)
  - Root outputs directory used by downstream scripts.
- `--figures-dir` (default: `figures`)
  - Destination for generated figure/table assets.
- `--csv` (default: `data/CSP_UBQ.csv`)
  - Main CSP input CSV path.
- `--include-placeholders`
  - Also run placeholder/optional scripts that may not yet emit final outputs.
- `--stop-on-error`
  - Exit immediately on first failing sub-script (default behavior is continue and report all failures).

## Prerequisites

- Pipeline outputs available under `--outputs-dir` (including per-target artifacts like `master_alignment.csv`).
- `confusion_matrix_per_system.csv` under `--outputs-dir` (typically produced by pipeline/follow-on analysis).
- `data/CSP_UBQ_ph0.5_temp5C.csv` for SI Table S1 generation.
- Additional domain-selection CSVs required by specific table scripts (for example `targets_CBP.csv`, etc., when expected by those scripts).

## Notes

- If a sub-script file is missing, the wrapper logs a skip for that script.
- PyMOL is not required by the wrapper itself, but some upstream artifacts it depends on may require PyMOL during their own generation.
