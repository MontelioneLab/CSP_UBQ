# Supplemental Figures and Tables Wrapper

This page documents `scripts/create_all_supplement_figures_and_tables.py`, the wrapper that runs the supplementary figure/table generators in sequence and rebuilds `SI_documents/SI merged.pdf`.

## Prerequisites

1. **Pair sync (required before regeneration after apo/holo edits)**  
   Treat `data/CSP_UBQ_ph0.5_temp5C.csv` as the source of truth for apo/holo pair edits, then propagate into other data CSVs:

   ```bash
   python scripts/sync_pairs_from_ph05.py
   python scripts/apo_holo_exp_conditions.py --fetch
   ```

   Do **not** overwrite `CSP_UBQ_ph0.5_temp5C.csv` via `filter_csp_ubq_by_buffer.py` as part of this sync.

2. **Pipeline outputs** under `--outputs-dir` (including per-target `master_alignment.csv` / `1d_analysis.csv` and `confusion_matrix_per_system.csv`):

   ```bash
   python scripts/pipeline.py \
     --input data/CSP_UBQ.csv \
     --out outputs \
     --no-case-study
   ```

3. Domain-selection CSVs for ST7–ST9 (`data/targets_BET_ET.csv`, `targets_TFIIH.csv`, `targets_ubiquitin.csv`).

4. **LaTeX toolchain** (`latexmk` or `pdflatex`) and **pypdf** for the final SI PDF merge.

## Usage

```bash
python scripts/create_all_supplement_figures_and_tables.py \
  --outputs-dir outputs \
  --figures-dir figures \
  --csv data/CSP_UBQ.csv
```

Optional: `--stop-on-error` exits on the first failing sub-script (default: continue and report all failures).

## SI figure scripts (S1–S27)

| SI Fig. | Script |
|--------|--------|
| S1–S8 | `scripts/create_si_figs_s1_s9.py` |
| S9 | `scripts/create_si_st10_fig_s9_dissimilar_conditions.py` |
| S10 | `scripts/create_si_fig_s10.py` (HA/CA confusion histograms) |
| S11 | `scripts/create_si_fig_s11.py` (CA-inclusive confusion histograms) |
| S12 | `scripts/create_si_fig_s12.py` |
| S13 | `scripts/create_si_fig_s13.py` |
| S14 | `scripts/create_si_fig_s14.py` |
| S15 | `scripts/create_si_fig_s15.py` |
| S16 | `scripts/create_si_fig_s16.py` (static PDB search screenshot) |
| S17 | `scripts/create_si_fig_s17.py` |
| S18 | `scripts/create_si_fig_s18.py` |
| S19 | `scripts/create_si_fig_s19.py` |
| S20 | `scripts/create_si_fig_s20.py` |
| S21 | `scripts/create_si_fig_s21.py` |
| S22 | `scripts/create_si_fig_s22.py` |
| S23 | `scripts/create_si_fig_s23.py` |
| S24 | `scripts/create_si_fig_s24.py` |
| S25 | `scripts/create_si_fig_s25.py` |
| S26 | `scripts/create_si_fig_s26.py` |
| S27 | `scripts/create_si_fig_s27.py` (2FIN case-study z panel) |

The wrapper invokes each of these in order before compiling/merging the SI PDF.

## What it does

1. Regenerates SI figure/table assets into `figures/` (`SF*.png`, `ST*.tex`, `SE*.png`).
2. Compiles `SI_documents/tex/main.tex` (which `\input`s `figures/ST1`–`ST10`) plus the C-term figures master via `scripts/build_si_merged_pdf.py`.
3. Merges PDFs in this order into `SI_documents/SI merged.pdf` (and copies `SI_merged.pdf` at the repo root):

```text
SI N term.pdf (title)
→ SI_TOC.pdf (from tex/si_toc.tex)
→ Supplementary Text.pdf
→ CSP_UBQ_SUPPL_TABLES.pdf   (from tex/main.tex → ST1–ST10)
→ SI C term.pdf
```

Compile the tables portion alone:

```bash
cd SI_documents/tex && latexmk -pdf -interaction=nonstopmode main.tex
```

`SI N term.pdf` title page is kept; its stale TOC page is replaced by a compiled
`SI_documents/tex/si_toc.tex` listing Tables S1–S10, Figs S1–S27, and Equations S1–S3.
`SI_documents/Supplementary Text.pdf` is inserted immediately after that TOC.
The C term is rebuilt from regenerated SF10–SF27 (and SE1–SE3) figures; SF16 uses
`figures/SF16_pdb_search.png` when present. References are appended from
`SI_documents/SI_C_term_references.pdf` only (not the last N pages of a prior C-term,
which previously re-introduced a duplicate SF19).

SI Fig. S10 uses the pH/temp-matched subset that also has apo and holo HA and CA
shifts (`csp_table_HA_CA.csv`; $n=101$ on the current outputs tree).
Analysis figures S11, S12, S15, and S18–S21 default to the $n=145$ buffer-similar
cohort (`data/CSP_UBQ_ph0.5_temp5C.csv`). S13/S14 use the rematch-excluded
pH/temp-matched list (`ph05_n137_targets.csv`, $n=136`).
SI Fig. S24 uses the five selected apo–apo control panels
`figures/SF24{A–E}_apo_apo_*.png`
(sources: `figures/selected_apo_apo_controls/{query}_{match}.png`,
falling back to
`outputs/apo_apo_matches/{query}_{match}/aggregate_csp_hsqc_grid.png`).
Panels: A 52080/52079, B 28070/28071, C 34000/34001, D 34394/6354, E 17769/51725.
SI Fig. S25 uses the precomputed terminal-anchor vs global-offset panel
`figures/SF25_terminal_anchor_vs_global.png`
(source: `outputs/hsqc_overlay_anchor_vs_global_all.png`).
SI Fig. S26 is composed by `scripts/create_si_fig_s26_csp_vs_distance.py` from
`p_significant_vs_ca_distance_ph05.png` (panel a) and the atom–atom
`max(0.05 ppm, cleaned mean)` CSP-$z$ scatter (panel b; bottom-left of the
CA vs any-atom $2\times 2$).
SI Fig. S27 copies the case-study z panel from
`outputs/2FIN_6809/2FIN_case_study_z.png`.

## CLI options

- `--outputs-dir` (default: `outputs`) — root outputs directory used by downstream scripts.
- `--figures-dir` (default: `figures`) — destination for generated figure/table assets.
- `--csv` (default: `data/CSP_UBQ.csv`) — main CSP input CSV path.
- `--stop-on-error` — exit immediately on first failing sub-script.

## Notes

- If a sub-script file is missing, the wrapper logs a skip for that script.
- PyMOL is not required by the wrapper itself, but some upstream artifacts it depends on may require PyMOL during their own generation.
- CSV-path behavior is mixed by design in the current wrapper:
  - `create_si_table_s1.py` always uses `data/CSP_UBQ_ph0.5_temp5C.csv` (the wrapper does not forward `--csv` to this script).
  - `create_si_fig_s17.py` regenerates `outputs/buffer_threshold_sweep/sweep_metrics.csv` plus the `heatmap_n.png` and `heatmap_pct_allosteric.png` panel assets before composing `figures/SF17_buffer_sweep.png`.
  - `create_si_fig_s20.py` receives `--outputs-dir`, `--input` (`<outputs-dir>/confusion_matrix_per_system.csv`), and `--output-image`.
  - `create_custom_selection_latex_tables.py` reads selection CSVs from `outputs/si_figs_s1_s9_aux/` (written by `create_si_figs_s1_s9.py`) plus domain CSVs under `data/`.
- Rebuild only the merged PDF (after figures/tables already exist):

  ```bash
  python scripts/build_si_merged_pdf.py
  ```
