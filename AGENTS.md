# AGENTS.md — CSP_UBQ

This file orients AI assistants and human contributors working in this repository.

## What this project is

A Python pipeline for analyzing **Chemical Shift Perturbations (CSPs)** from BMRB NMR-STAR data and PDB structures: sequence alignment, HN-weighted CSPs, SASA-based backbone occlusion, interaction filters, plots, and PyMOL visualization scripts.

For the public overview, methods summary, and how to run the main scripts, see [README.md](README.md). For CLI flags, outputs, grid referencing, and troubleshooting, see [docs/pipeline_reference.md](docs/pipeline_reference.md).

## Repository map

| Path | Role |
|------|------|
| [scripts/](scripts/) | All Python code. Main entry: [scripts/pipeline.py](scripts/pipeline.py). Core modules include [scripts/csp.py](scripts/csp.py), [scripts/bmrb_io.py](scripts/bmrb_io.py), [scripts/config.py](scripts/config.py) (e.g. `Referencing`, grid referencing). |
| [outputs/](outputs/) | Generated per-target data under `outputs/<holo_pdb>/`. Treat as build artifacts unless the task is to regenerate outputs or fix pipeline behavior. |
| [CS_Lists/](CS_Lists/) | NMR-STAR / chemical shift list inputs. |
| [figures/](figures/) | Publication and supplementary figures (often built with `scripts/create_*` helpers). |
| Root CSVs (e.g. `CSP_UBQ.csv`) | Pipeline inputs; column definitions are in [docs/pipeline_reference.md](docs/pipeline_reference.md). |

## Environment and how to run

Install dependencies from [requirements.txt](requirements.txt):

```bash
pip install -r requirements.txt
```

Alternatively use [environment.yml](environment.yml) with conda/mamba (`conda env create -f environment.yml`).

Typical pipeline run:

```bash
python scripts/pipeline.py --input CSP_UBQ.csv --out outputs
```

Additional flags (`--ids`, `--workers`, `--no-case-study`, metadata annotation, etc.) and verbose logging (`CSP_VERBOSE=1`) are documented in [docs/pipeline_reference.md](docs/pipeline_reference.md).

Grid referencing options and defaults live in [scripts/config.py](scripts/config.py); [docs/pipeline_reference.md](docs/pipeline_reference.md) describes the grid workflow in detail.

## Testing

Tests live under [scripts/test_pipeline.py](scripts/test_pipeline.py) with shared config in [scripts/conftest.py](scripts/conftest.py). From the repository root:

```bash
pytest scripts/test_pipeline.py
```

Integration-marked tests:

```bash
pytest scripts/test_pipeline.py -m integration
```

Omit integration tests:

```bash
pytest scripts/test_pipeline.py -m "not integration"
```

## Notes for agents

- **Case-study figures** (on by default; disable with `--no-case-study`) expect a **PyMOL GUI** and manual steps (orientation, F5 to save the view) when views are not cached. Do not assume a headless environment will succeed.
- **Avoid large unrelated diffs** in `outputs/`, `figures/`, and other generated or binary assets unless the user explicitly asked to regenerate or change them.
- Prefer focused code changes in [scripts/](scripts/) that match existing style and patterns; extend existing functions rather than duplicating logic.

For full project documentation, start from [README.md](README.md); use [docs/pipeline_reference.md](docs/pipeline_reference.md) for pipeline implementation details.
