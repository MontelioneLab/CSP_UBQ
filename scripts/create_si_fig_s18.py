#!/usr/bin/env python3
"""
Legacy helper: BMRB-scrape target list (optional).

For the main supplement, SI Table S14 / SI Fig. S13 use
create_si_st14_fig_s13_dissimilar_conditions.py (dissimilar apo/holo conditions).

This script still runs analyze_targets.py for targets_BMRB.csv and writes
figures/SF13_bmrb.png if you maintain that CSV.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def load_targets_csv(targets_csv: Path) -> list[str]:
    targets: set[str] = set()
    with open(targets_csv, "r", newline="", encoding="utf-8") as f:
        import csv
        reader = csv.DictReader(f)
        for row in reader:
            holo_pdb = (row.get("holo_pdb") or "").strip().lower()
            if holo_pdb:
                targets.add(holo_pdb)
    return sorted(targets)


def compose_two_panel_figure(
    panel_a_path: Path,
    panel_b_path: Path,
    output_image: Path,
) -> None:
    """Create a two-panel figure: A (left) = stacked_hist, B (right) = confusion_stacked_hist."""
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
        description="Legacy: confusion histograms for targets_BMRB.csv (not the default SI Fig. S13)."
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--targets-csv", type=Path, default=Path("data/targets_BMRB.csv"))
    parser.add_argument("--aux-dir", type=Path, default=Path("outputs") / "si_figs_s1_s12_aux")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else repo_root / args.targets_csv
    outputs_dir = repo_root / args.outputs_dir
    figures_dir = repo_root / args.figures_dir
    aux_dir = repo_root / args.aux_dir

    if not targets_csv.exists():
        print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
        return 1

    targets = load_targets_csv(targets_csv)
    if not targets:
        print("No targets found in targets CSV.", file=sys.stderr)
        return 1

    file_token = "bmrb"
    sf_id = "SF13"
    analyze_script = repo_root / "scripts" / "analyze_targets.py"
    stacked_histogram_image = figures_dir / f"{sf_id}_{file_token}_stacked_hist.png"
    confusion_histogram_image = figures_dir / f"{sf_id}_{file_token}_confusion_stacked_hist.png"
    summary_dir = aux_dir / f"si_figs_s1_s12_{file_token}_summary_statistics"

    cmd = [
        args.python,
        str(analyze_script),
        "--outputs-dir",
        str(outputs_dir),
        "--targets-csv",
        str(targets_csv),
        "--output-image",
        str(aux_dir / f"si_figs_s1_s12_{file_token}.png"),
        "--summary-csv",
        str(aux_dir / f"si_figs_s1_s12_{file_token}_summary.csv"),
        "--histogram-image",
        str(aux_dir / f"si_figs_s1_s12_{file_token}_hist.png"),
        "--stacked-histogram-image",
        str(stacked_histogram_image),
        "--confusion-matrix-stacked-histogram-image",
        str(confusion_histogram_image),
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

    subprocess.run(cmd, check=True, cwd=str(repo_root))

    panel_a = stacked_histogram_image
    panel_b = confusion_histogram_image
    collated_path = figures_dir / f"{sf_id}_{file_token}.png"
    if panel_a.exists() and panel_b.exists():
        compose_two_panel_figure(panel_a, panel_b, collated_path)
        print(f"BMRB collated figure saved to {collated_path.resolve()}")
    else:
        print("Warning: stacked or confusion histogram not created.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
