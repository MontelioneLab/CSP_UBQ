# CSP pipeline — reference

This document supplements the [project README](../README.md) with **CLI details, I/O layout, configuration, and troubleshooting** for the CSP_UBQ analysis pipeline.

## Overview

The pipeline processes pairs of apo and holo protein structures to:

- Extract chemical shift data from BMRB NMR-STAR files  
- Align sequences and compute CSPs using the HN-weighted equation  
- Analyze backbone NH occlusion and interaction-based binding-site criteria  
- Map CSPs to PDB structures for visualization  
- Generate CSV tables, plots, and PyMOL scripts  

### Features

- Multi-saveframe support for BMRB entries  
- Needleman–Wunsch global alignment with best-score selection  
- SASA occlusion and additional interaction filters  
- Concurrent processing with configurable workers  
- Verbose logging via `CSP_VERBOSE=1`  

## Installation

Install dependencies from the repository root:

```bash
pip install -r requirements.txt
```

Optional: **PyMOL** for running generated `.pml` scripts or `--case-study` mode (GUI).

## Usage

### Basic usage

```bash
python scripts/pipeline.py --input data/CSP_UBQ.csv --out outputs
```

Optional metadata columns (`ec_classes`, `scope_fold_type`):

```bash
python scripts/pipeline.py --input data/CSP_UBQ.csv --out outputs --run-csv-metadata-annotation
```

Subset of rows (BMRB IDs as in your CSV usage, if applicable):

```bash
python scripts/pipeline.py --input data/CSP_UBQ.csv --ids 18251,4700 --out outputs
```

Single holo PDB (interactive PyMOL for case-study figures unless disabled):

```bash
python scripts/pipeline.py --input data/CSP_UBQ.csv --holo-pdb 1H8B --out outputs
```

By default the pipeline attempts case-study figure generation per target. Use `--no-case-study` to skip. With PyMOL available, it writes `outputs/<target>/<pdb_id>_case_study.png`: open `color_csp_mask.pml`, set the view, press **F5** to save the camera, then the pipeline reuses it for ray-traced panels.

### Verbose and parallel

```bash
CSP_VERBOSE=1 python scripts/pipeline.py --input data/CSP_UBQ.csv --out outputs
python scripts/pipeline.py --input data/CSP_UBQ.csv --out outputs --workers 4
```

### Apo BMRB 17769 sequence alignment (utility)

```bash
python scripts/align_apo_bmrb_sequences.py
```

Reads `data/CSP_UBQ.csv`, aligns apo sequences to the `apo_bmrb == 17769` reference (`match_seq`), writes `apo_bmrb_17769_alignments.csv` with scores, overlap, identity, approximate E-value.

### Example

```bash
CSP_VERBOSE=1 python scripts/pipeline.py --input data/CSP_UBQ.csv --ids 18251,4700 --out outputs --workers 2
```

## Referencing method (mean vs grid)

The pipeline can determine optimal holo \(^1\text{H}\) / \(^{15}\text{N}\) offsets by **grid search**, maximizing the count of aligned residues with CSP below a cutoff after referencing.

Configure in `scripts/config.py` via the `Referencing` dataclass:

- `method`: `"mean"` (legacy) or `"grid"` (default `"grid"`)  
- `grid_h_min`, `grid_h_max`, `grid_h_step` (defaults include step **0.01** ppm)  
- `grid_n_min`, `grid_n_max`, `grid_n_step` (defaults include step **0.05** ppm)  
- `grid_cutoff` (default **0.05** ppm)  
- `cache_results`, `save_heatmap`  

Programmatic entry points in `scripts/csp.py`:

- `compute_csp_A(..., enable_referencing=True, referencing_method=None, grid_params=None, target_id=None)`  
- `compute_csp_from_aligned_sequences(..., enable_referencing=True, referencing_method=None, grid_params=None, target_id=None)`  

With `referencing_method="grid"` and `target_id` set, results cache to `outputs/<target_id>/offset_grid_<params>.csv` and heatmaps to `offset_grid_<params>.png`.

Override example:

```python
grid_params = {
    "h_min": -0.12, "h_max": 0.12, "h_step": 0.01,
    "n_min": -1.2,  "n_max": 1.2,  "n_step": 0.05,
    "cutoff": 0.05,
}
```

## Input format

The pipeline expects a CSV (e.g. `data/CSP_UBQ.csv`) with at least:

- `apo_bmrb` — BMRB ID for apo  
- `holo_bmrb` — BMRB ID for holo  
- `holo_pdb` — PDB ID for holo  
- `match_seq` — sequence matching info (alignment uses pipeline logic; column may be ignored depending on path)  

See your project CSV for any additional columns.

## Output files

Per processed entry, under `outputs/{holo_pdb}/` (or suffixed directory if duplicates):

### Data

- `csp_table.csv` — CSPs, masks, interaction/occlusion fields  
- `occlusion_analysis.csv` — SASA occlusion details where applicable  
- `alignment.txt` — alignment notes  

### Visualizations

- `color_csp_mask.pml`, `color_occlusion.pml`, and other `.pml` / plots as produced by the run  

### Parsed data

- `parsed/` — `{bmrb_id}_combined.csv`, per-saveframe CSVs  

A full run also aggregates **`outputs/confusion_matrix_per_system.csv`** for downstream figure scripts.

## CSP calculation (implementation)

The HN CSP uses the `Compute` settings in `scripts/config.py`:

$$\text{CSP} = \sqrt{\tfrac{1}{2}\left(\Delta\delta_H^2 + (w_N\,\Delta\delta_N)^2\right)}$$

with default **\(w_N = 0.14\)** (`csp_delta_n_scale`). CA-inclusive CSPs use a separate weight (`csp_delta_ca_scale`).

### Significance threshold (defaults)

Iterative outlier removal in `Thresholds`:

1. Mean and SD of CSPs; remove values with Z-score above `outlier_z_score` (default **3.0**).  
2. Repeat until stable or limits hit (`max_outlier_iterations`, `max_outlier_fraction`).  
3. Final cutoff: mean + `significance_z_score` × SD of the cleaned set (default **0.0** × SD → **mean** of cleaned CSPs).  

CLI overrides include `--outlier-z`, `--significance-z`, `--max-outlier-iterations`, `--max-outlier-fraction`, `--absolute-cutoff`.

## SASA occlusion analysis

Backbone N/H SASA is compared for holo complex vs protein-only (mdtraj Shrake–Rupley). Thresholds live in `scripts/config.py` (`SASAAnalysis.sasa_threshold`; default **0.0** in code — confirm before interpreting “occluded” flags).

## Pipeline architecture

| Module | Role |
|--------|------|
| `pipeline.py` | CLI orchestrator |
| `bmrb_io.py` | BMRB download / parse |
| `align.py` | Sequence alignment |
| `csp.py` | CSP and significance |
| `rcsb_io.py` | PDB download / parse |
| `sasa_analysis.py` | SASA occlusion |
| `interaction_analysis.py` | H-bond, distance, and related filters |
| `visualize.py` | Plots and PyMOL scripts |
| `config.py` | Central configuration |

### Data flow

1. Read/annotate input CSV  
2. Fetch and parse apo/holo BMRB  
3. Global sequence alignment  
4. CSP computation with referencing  
5. Holo PDB mapping  
6. SASA and interaction analyses  
7. Write tables, figures, PyMOL scripts; refresh confusion-matrix summary  

## Configuration

Edit `scripts/config.py` for paths, network timeouts, alignment scoring, thresholds, referencing grid, SASA and distance cutoffs.

## PyMOL

Example:

```bash
pymol outputs/<holo_pdb>/color_occlusion.pml
```

Scripts typically color the receptor chain, binding-site-related residues, and ligand (e.g. cyan) according to the generator in `visualize.py`.

## Troubleshooting

1. **Dependencies** — `pip install -r requirements.txt`  
2. **Network** — BMRB / RCSB availability; retries in `config.py`  
3. **Parsing** — `CSP_VERBOSE=1`  
4. **Memory / time** — reduce `--workers`  

```bash
CSP_VERBOSE=1 python scripts/pipeline.py --input data/CSP_UBQ.csv --ids 18251 --out outputs
```

## Tests

```bash
pytest scripts/test_pipeline.py
pytest scripts/test_pipeline.py -m integration
pytest scripts/test_pipeline.py -m "not integration"
```

## Repository layout (abbreviated)

```
CSP_UBQ/
├── scripts/           # Python package and CLIs
├── docs/              # This reference
├── CS_Lists/          # BMRB cache
├── PDB_FILES/         # PDB cache
├── outputs/           # Pipeline outputs
├── figures/           # Generated or committed figures
├── data/              # Curated input CSV tables (CSP_UBQ.csv, targets_*.csv, etc.)
├── requirements.txt
├── environment.yml
└── README.md
```
