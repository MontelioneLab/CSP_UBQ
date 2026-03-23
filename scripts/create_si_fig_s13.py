#!/usr/bin/env python3
"""
Generator script filename (s13) ≠ SI index: this produces SI Fig. S14.

SI Fig. S14 — CA-inclusive CSP confusion matrix histograms (two-panel image).

- Panel A: stacked histogram of significant residues (TP/FP only)
- Panel B: stacked confusion-matrix histogram (TN/FP/FN/TP)

Output:
- figures/SF14_ca_inclusive.png

Default targets list: CSP_UBQ_ph0.5_temp5C.csv (buffer-filtered subset). Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

try:
    from .config import classification_colors
except Exception:
    _root = Path(__file__).resolve().parent.parent
    if str(_root) not in sys.path:
        sys.path.insert(0, str(_root))
    from scripts.config import classification_colors  # type: ignore

try:
    from .analyze_targets_ca import (
        collect_results,
        load_targets_from_csv,
        render_confusion_matrix_stacked_histogram,
        render_stacked_histogram,
    )
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_targets_ca import (  # type: ignore
        collect_results,
        load_targets_from_csv,
        render_confusion_matrix_stacked_histogram,
        render_stacked_histogram,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S14 (CA-inclusive CSP confusion matrix histograms)."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root outputs directory containing per-target folders.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with holo_pdb target IDs (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument("--output-image", type=Path, default=Path("figures") / "SF14_ca_inclusive.png")
    return parser.parse_args()


def compose_two_panel_figure(panel_a_path: Path, panel_b_path: Path, output_image: Path) -> None:
    image_a = mpimg.imread(panel_a_path)
    image_b = mpimg.imread(panel_b_path)

    fig, axes = plt.subplots(2, 1, figsize=(10, 12))
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

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.02)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else project_root / args.targets_csv
    output_image = (
        args.output_image if args.output_image.is_absolute() else project_root / args.output_image
    )

    if not outputs_dir.exists():
        print(f"Error: outputs directory does not exist: {outputs_dir}", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    allowed_targets = load_targets_from_csv(targets_csv)
    _, distances, confusion_records = collect_results(outputs_dir, allowed_targets)

    positive_count = sum(1 for record in distances if record.is_predicted_positive)
    negative_count = sum(1 for record in distances if not record.is_predicted_positive)

    output_image.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="si_fig_s13_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        panel_a_path = tmp_dir_path / "panel_a.png"
        panel_b_path = tmp_dir_path / "panel_b.png"

        render_stacked_histogram(
            distances,
            panel_a_path,
            positive_count,
            negative_count,
            tp_color=classification_colors.TP,
            show_title=False,
            ylabel="Number of Residues",
            bold_axes=False,
            axis_label_fontsize=20,
            legend_fontsize=14,
        )
        render_confusion_matrix_stacked_histogram(
            confusion_records,
            panel_b_path,
            tp_color=classification_colors.TP,
            show_title=False,
            bold_axes=False,
            axis_label_fontsize=20,
            legend_fontsize=14,
        )
        if not panel_a_path.is_file() or not panel_b_path.is_file():
            print(
                "Cannot compose SI Fig. S14: one or both panel PNGs were not written "
                "(no CA-distance data for selected targets / outputs).",
                file=sys.stderr,
            )
            return 1
        compose_two_panel_figure(panel_a_path, panel_b_path, output_image)

    print(f"SI Fig. S14 saved to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
