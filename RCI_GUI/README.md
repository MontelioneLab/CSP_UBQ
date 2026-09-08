# RCI GUI and CLI

This repository contains the active Random Coil Index workflow, a Qt desktop GUI for running it, regression tests, and archived legacy scripts and docs.

The main entry points are:

- GUI: [scripts/rci_gui.py](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/scripts/rci_gui.py)
- Active analysis script: [scripts/rci_v_1c_PyNMR-STAR.py](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/scripts/rci_v_1c_PyNMR-STAR.py)
- Legacy reference script: [archive/rci_v_1c.py](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/archive/rci_v_1c.py)
- Original project README: [archive/README](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/archive/README)

## Project layout

- `scripts/` contains the active runnable code.
- `inputs/` contains example and working input files for local runs.
- `outputs/` contains generated artifacts, including GUI run folders.
- `tests/` contains regression tests.
- `archive/` contains preserved legacy scripts and the original README.
- `RCI-Wishart-calculation/` contains papers, examples, and older supporting material from the historical project layout.

## Supported input formats

The active script accepts:

- NMR-STAR 3 files, parsed with `pynmrstar`
- SHIFTY-style tables with a header like `#NUM AA HA CA CB CO N HN`

The current PyNMR-STAR path does not support older NMR-STAR 2.1 files directly. If you need to compare behavior against the old parser, use the archived legacy script in `archive/`.

For multi-entity NMR-STAR 3 inputs, both the converter and the main analysis workflow can now target a specific `Entity_ID` instead of collapsing overlapping residue numbering across chains or constructs.

## Environment setup

The easiest setup path is Conda, using [environment.yaml](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/environment.yaml).

```bash
conda env create -f environment.yaml
conda activate rci
```

If you already created the environment and want to refresh it:

```bash
conda env update -f environment.yaml --prune
conda activate rci
```

This environment installs:

- `python`
- `numpy`
- `matplotlib` for modern PNG and SVG plots
- `pyside6` for the desktop GUI
- `gnuplot` for legacy PS, GIF, and JPG plot generation
- `pytest` for tests
- `pynmrstar` for NMR-STAR 3 input

## Plot generation

The active script now supports two plotting paths:

- Legacy `gnuplot` output, preserved for compatibility
- New `matplotlib` output, used for modern `png` and `svg` files

The `gnuplot` path still generates legacy image formats through flags like `-gif`, `-jpg`, and `-ps`. The new `matplotlib` path is enabled with `-mpl` and writes cleaner `png` and `svg` outputs in parallel with the existing numeric text outputs.

The Qt GUI now prefers the `matplotlib` SVG files for display instead of the older `gnuplot` GIF/JPG outputs.

## Default RCI behavior

The active local workflow now defaults to:

- `incomplete_data_use = 0`
- `end_corr = 2`
- `points_to_delete = 0`
- `gap_limit = 0`
- `sw_neighbor_flag = 1`

In practice, this means local runs exclude residues that do not have a complete set of supported chemical shifts, use the additive terminal correction implemented by `end_effect2`, keep terminal points instead of deleting them up front, disable gap-filling of nearby residues, and enable the combined Schwarz-Wang neighboring residue correction.

This default is aligned across:

- direct runs of [scripts/rci_v_1c_PyNMR-STAR.py](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/scripts/rci_v_1c_PyNMR-STAR.py)
- GUI-launched runs from [scripts/rci_gui.py](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/scripts/rci_gui.py)

If you want to remove artificially weak terminal residues instead of correcting them additively, use `-end_corr5`. That mode is more aggressive and can improve agreement with some server traces, but it may drop terminal residues from the output instead of preserving them.

## Running the GUI

From the project root:

```bash
conda activate rci
python scripts/rci_gui.py
```

## Using the GUI

The GUI is a wrapper around the active RCI script.

1. Click `Choose Input File`.
2. Pick an input file from `inputs/` or another location on disk.
3. Click `Run Script`.
4. Wait for the status to change from `running` to `completed`, `failed`, or `cancelled`.
5. Review the generated SVG plots, logs, and text outputs in the right-hand panel.

If the selected input is a multi-entity NMR-STAR 3 file, the GUI also exposes an entity selector on the analysis tab and the format-conversion tab so you can choose which entity to analyze or export.

GUI runs currently launch the script with the equivalent of:

```bash
python scripts/rci_v_1c_PyNMR-STAR.py -b <input-file> -mpl -r 4 -d 1 -no_i -end_corr2 -gapfill2 -dynamr
```

## What the GUI stores

Each GUI run stores a self-contained workspace under:

`outputs/gui_runs/<run-id>/`

A run directory includes:

- `run.json` with run metadata and selected options
- `command.log` with the executed command
- `stdout.log`
- `stderr.log`
- `workspace/` containing the copied input and generated analysis files

Inside the GUI:

- `Previous runs` lists earlier executions
- the image area previews generated SVG files when present, otherwise other discovered image artifacts
- the `Logs` tab shows run logs
- the `Text Outputs` tab shows generated `.RCI.txt`, `.S2.txt`, `.MD_RMSD.txt`, and `.NMR_RMSD.txt` files when present

## Running from the command line

You can run the active script directly:

```bash
conda activate rci
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/bmr15476_3.str.txt
```

By default, direct CLI runs now behave as if incomplete residues are excluded and terminal correction mode `2` is selected. If you want a stronger terminal filter that removes artificially weak terminal residues, use:

```bash
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/bmr15476_3.str.txt -end_corr5
```

For a SHIFTY input:

```bash
conda activate rci
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/bmr15476_shifty.txt
```

For a specific entity in a multi-entity NMR-STAR 3 file:

```bash
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/30786_3.str.txt -entity 1
```

The active SHIFTY conversion path now uses the simpler interoperable schema:

- SHIFTY files contain a single `HA` column
- for glycine, converting from SHIFTY back to NMR-STAR duplicates `HA` into both `HA2` and `HA3`
- this loses the distinction between `HA2` and `HA3`, but keeps the files compatible with downstream tools that expect the simpler SHIFTY layout

For modern `matplotlib` plots:

```bash
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/bmr15476_3.str.txt -mpl
```

For legacy `gnuplot` GIF output:

```bash
python scripts/rci_v_1c_PyNMR-STAR.py -b inputs/bmr15476_3.str.txt -gif
```

## Example inputs

Example files already included in [inputs/](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs):

- [inputs/bmr15476_3.str.txt](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs/bmr15476_3.str.txt)
- [inputs/bmr15476_shifty.txt](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs/bmr15476_shifty.txt)
- [inputs/bmr15476_21.str.txt](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs/bmr15476_21.str.txt)
- [inputs/30786_3.str.txt](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs/30786_3.str.txt)
- [inputs/PyJCScorr](/Users/tiburon/Desktop/GTM_Research/01_Active_Projects/RCI/inputs/PyJCScorr)

For the PyNMR-STAR-backed script, prefer the NMR-STAR 3 or SHIFTY examples.

## Tests

Run the regression tests from the project root:

```bash
conda activate rci
pytest -q
```

The test suite currently covers:

- GUI artifact discovery helpers
- GUI entity selection and format conversion helpers
- SVG preference in the GUI image selector
- NMR-STAR 3 compatibility against legacy parser outputs
- SHIFTY compatibility against the `nmrstar3toSHIFTY.py` conversion workflow
- multi-entity conversion and entity-specific output equivalence for `30786`
- clear failure behavior for unsupported NMR-STAR 2.1 input on the PyNMR-STAR path

## Troubleshooting

- `PySide6 is required for the Qt GUI`
  Activate the Conda environment from `environment.yaml` before launching the GUI.

- `Unable to parse NMR-STAR file ... with pynmrstar`
  The file is likely an older NMR-STAR 2.1 file or malformed for the PyNMR-STAR parser. Use an NMR-STAR 3 file or a SHIFTY export.

- No outputs appear in the GUI
  Open the selected run and inspect `stdout.log`, `stderr.log`, and `command.log`.

- No modern SVG or PNG files appear
  Check whether `matplotlib` is installed in the environment. The script logs a warning and skips modern plot generation if `matplotlib` cannot be imported.

- Legacy plot files fail to generate
  Check whether `gnuplot` is installed in the active environment.

- The selected file disappears or cannot be opened
  The GUI copies the chosen input into the run workspace before execution. Re-run with a stable input file path if needed.
