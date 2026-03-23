#!/usr/bin/env python3
"""
SI Table S14 + SI Fig. S13 — targets in CSP_UBQ.csv but not in the buffer-similar
subset (CSP_UBQ_ph0.5_temp5C.csv).

Writes:
  - A filtered CSP CSV (full study rows with dissimilar apo/holo conditions)
  - targets CSV (holo_pdb) for analyze_targets.py
  - figures/ST14_dissimilar_apo_holo_conditions.tex (via create_csp_latex_table.py)
  - figures/SF13_dissimilar_apo_holo_conditions.png (two-panel confusion histograms)

Row identity matches CSP_UBQ vs the filtered file on (holo_pdb, apo_bmrb) after the
same normalization used for LaTeX tables.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

try:
    from .create_csp_latex_table import _bmrb_to_str, _normalize_pdb_id
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.create_csp_latex_table import _bmrb_to_str, _normalize_pdb_id  # type: ignore

SELECTION_NAME = "dissimilar_apo_holo_conditions"
FILE_TOKEN = "dissimilar_apo_holo_conditions"
SF_ID = "SF13"


def _row_key_series(df: pd.DataFrame) -> pd.Series:
    holo = df["holo_pdb"].map(_normalize_pdb_id)
    apo = df["apo_bmrb"].map(_bmrb_to_str)
    return holo + "|" + apo


def derive_dissimilar_df(full: pd.DataFrame, similar: pd.DataFrame) -> pd.DataFrame:
    """Rows in full whose (holo_pdb, apo_bmrb) key is not in similar."""
    sim_keys = set(_row_key_series(similar))
    keys = _row_key_series(full)
    out = full.loc[~keys.isin(sim_keys)].copy()
    return out


def load_targets_holo_pdbs(dissimilar: pd.DataFrame) -> list[str]:
    """PDB strings as in CSP_UBQ (strip only) so they match outputs/ directory names."""
    seen: set[str] = set()
    ordered: list[str] = []
    for v in dissimilar["holo_pdb"].tolist():
        h = str(v).strip()
        if not h or h.lower() == "nan" or h in seen:
            continue
        seen.add(h)
        ordered.append(h)
    return ordered


def compose_two_panel_figure(
    panel_a_path: Path,
    panel_b_path: Path,
    output_image: Path,
) -> None:
    image_a = mpimg.imread(panel_a_path)
    image_b = mpimg.imread(panel_b_path)

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for ax in axes:
        ax.set_axis_off()

    axes[0].imshow(image_a)
    axes[1].imshow(image_b)

    axes[0].text(
        0.01,
        0.99,
        "A.",
        transform=axes[0].transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    axes[1].text(
        0.01,
        0.99,
        "B.",
        transform=axes[1].transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
    )

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="SI Table S14 + SI Fig. S13 (dissimilar apo/holo experimental conditions)."
    )
    parser.add_argument(
        "--full-csp-csv",
        type=Path,
        default=Path("CSP_UBQ.csv"),
        help="Full study table (default: CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--similar-csp-csv",
        type=Path,
        default=Path("CSP_UBQ_ph0.5_temp5C.csv"),
        help="Similar-conditions subset (default: CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument(
        "--confusion-csv",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Per-system confusion metrics.",
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--aux-dir",
        type=Path,
        default=Path("outputs") / "si_dissimilar_aux",
        help="Intermediate CSVs and analyze_targets aux outputs.",
    )
    parser.add_argument(
        "--output-tex",
        type=Path,
        default=Path("figures") / "ST14_dissimilar_apo_holo_conditions.tex",
    )
    parser.add_argument(
        "--subset-csv",
        type=Path,
        default=None,
        help="Optional path to write the dissimilar CSP rows (default: aux-dir).",
    )
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    full_path = args.full_csp_csv if args.full_csp_csv.is_absolute() else repo_root / args.full_csp_csv
    sim_path = args.similar_csp_csv if args.similar_csp_csv.is_absolute() else repo_root / args.similar_csp_csv
    confusion_path = (
        args.confusion_csv if args.confusion_csv.is_absolute() else repo_root / args.confusion_csv
    )
    outputs_dir = repo_root / args.outputs_dir
    figures_dir = repo_root / args.figures_dir
    aux_dir = repo_root / args.aux_dir
    out_tex = args.output_tex if args.output_tex.is_absolute() else repo_root / args.output_tex

    if not full_path.is_file():
        print(f"Error: full CSP CSV not found: {full_path}", file=sys.stderr)
        return 1
    if not sim_path.is_file():
        print(f"Error: similar-conditions CSV not found: {sim_path}", file=sys.stderr)
        return 1
    if not confusion_path.is_file():
        print(f"Error: confusion CSV not found: {confusion_path}", file=sys.stderr)
        return 1

    full_df = pd.read_csv(full_path, dtype=str)
    sim_df = pd.read_csv(sim_path, dtype=str)
    for col in ("holo_pdb", "apo_bmrb"):
        if col not in full_df.columns or col not in sim_df.columns:
            print(f"Error: CSP CSVs must contain columns {col!r}.", file=sys.stderr)
            return 1

    dissimilar = derive_dissimilar_df(full_df, sim_df)
    if dissimilar.empty:
        print("No dissimilar-conditions rows (full minus similar subset is empty).", file=sys.stderr)
        return 1

    aux_dir.mkdir(parents=True, exist_ok=True)
    subset_path = args.subset_csv
    if subset_path is None:
        subset_path = aux_dir / "CSP_UBQ_dissimilar_experimental_conditions.csv"
    elif not subset_path.is_absolute():
        subset_path = repo_root / subset_path
    subset_path.parent.mkdir(parents=True, exist_ok=True)
    dissimilar.to_csv(subset_path, index=False)
    print(f"Wrote {len(dissimilar)} dissimilar rows to {subset_path}")

    targets_csv = aux_dir / "targets_dissimilar_apo_holo_conditions.csv"
    holo_list = load_targets_holo_pdbs(dissimilar)
    with targets_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["holo_pdb"])
        w.writeheader()
        for h in holo_list:
            w.writerow({"holo_pdb": h})
    print(f"Wrote {len(holo_list)} holo_pdb targets to {targets_csv}")

    latex_script = repo_root / "scripts" / "create_csp_latex_table.py"
    tex_cmd = [
        args.python,
        str(latex_script),
        "--csp-csv",
        str(subset_path),
        "--confusion-csv",
        str(confusion_path),
        "--output",
        str(out_tex),
        "--all",
        "--selection-name",
        SELECTION_NAME,
    ]
    r = subprocess.run(tex_cmd, cwd=str(repo_root))
    if r.returncode != 0:
        print("create_csp_latex_table.py failed.", file=sys.stderr)
        return r.returncode

    analyze_script = repo_root / "scripts" / "analyze_targets.py"
    stacked = figures_dir / f"{SF_ID}_{FILE_TOKEN}_stacked_hist.png"
    confusion_hist = figures_dir / f"{SF_ID}_{FILE_TOKEN}_confusion_stacked_hist.png"
    summary_dir = aux_dir / f"si_figs_{FILE_TOKEN}_summary_statistics"
    an_cmd = [
        args.python,
        str(analyze_script),
        "--outputs-dir",
        str(outputs_dir),
        "--targets-csv",
        str(targets_csv),
        "--output-image",
        str(aux_dir / f"si_fig_{FILE_TOKEN}.png"),
        "--summary-csv",
        str(aux_dir / f"si_fig_{FILE_TOKEN}_summary.csv"),
        "--histogram-image",
        str(aux_dir / f"si_fig_{FILE_TOKEN}_hist.png"),
        "--stacked-histogram-image",
        str(stacked),
        "--confusion-matrix-stacked-histogram-image",
        str(confusion_hist),
        "--summary-dir",
        str(summary_dir),
        "--no-plot-titles",
        "--axis-fontsize",
        "16",
        "--tick-fontsize",
        "14",
        "--legend-fontsize",
        "14",
        "--legend-title-fontsize",
        "14",
        "--compact-legend-labels",
    ]
    r2 = subprocess.run(an_cmd, cwd=str(repo_root))
    if r2.returncode != 0:
        print("analyze_targets.py failed.", file=sys.stderr)
        return r2.returncode

    collated = figures_dir / f"{SF_ID}_{FILE_TOKEN}.png"
    if stacked.is_file() and confusion_hist.is_file():
        compose_two_panel_figure(stacked, confusion_hist, collated)
        print(f"SI Fig. S13 saved to {collated.resolve()}")
    else:
        print("Warning: stacked or confusion histogram PNGs missing.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
