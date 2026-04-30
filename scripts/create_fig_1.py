#!/usr/bin/env python3
"""
Create publication-ready Figure 1 panels as two separate files:

- figure_1_a.png: the chemical-shift offset grid-search heatmap for PDB 7JQ8,
  copied from `outputs/7jq8/offset_grid_H_-0.12_0.12_0.01__N_-1.2_1.2_0.05__C_0.05.png`.
- figure_1_b.png: a 3x3 confusion-matrix table summarizing the TP/FP/FN/TN
  classification scheme used throughout the study, with colors pulled from
  `scripts/config.py` so they stay consistent with the rest of the project.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import classification_colors
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors


DEFAULT_PANEL_A = (
    _REPO_ROOT
    / "outputs/7jq8/offset_grid_H_-0.12_0.12_0.01__N_-1.2_1.2_0.05__C_0.05.png"
)


def _draw_confusion_matrix(ax: plt.Axes) -> None:
    """Render the 3x3 confusion-matrix table onto `ax`, making it wider for better text fit."""

    header_bg = "#ffffff"
    header_fg = "#000000"
    cell_fg = "#ffffff"

    cells = [
        [("", header_bg, header_fg, "normal"),
         ("Orthosteric Residue", header_bg, header_fg, "normal"),
         ("Allosteric Residue", header_bg, header_fg, "normal")],
        [("Significant CSP", header_bg, header_fg, "normal"),
         ("True Positive (TP)", classification_colors.TP, cell_fg, "normal"),
         ("False Positive (FP)", classification_colors.FP, cell_fg, "normal")],
        [("No Significant CSP", header_bg, header_fg, "normal"),
         ("False Negative (FN)", classification_colors.FN, cell_fg, "normal"),
         ("True Negative (TN)", classification_colors.TN, cell_fg, "normal")],
    ]

    n_rows = len(cells)
    n_cols = len(cells[0])

    col_width = 3.0
    row_height = 0.5

    ax.set_xlim(0, n_cols * col_width)
    ax.set_ylim(0, n_rows * row_height)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    for r, row in enumerate(cells):
        y = n_rows - 1 - r
        for c, (text, bg, fg, weight) in enumerate(row):
            is_header = (r == 0) or (c == 0)
            edgecolor = "#000000"
            linewidth = 1.2
            ax.add_patch(
                Rectangle(
                    (c * col_width, y * row_height),
                    col_width,
                    row_height,
                    facecolor=bg,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )
            )
            ax.text(
                c * col_width + col_width / 2,
                y * row_height + row_height / 2,
                text,
                ha="center",
                va="center",
                color=fg,
                fontsize=13 if is_header else 14,
                fontweight=weight,
                wrap=True,
            )


def save_panel_a(panel_a_path: Path, output_path: Path) -> None:
    """Copy the source Panel A PNG to `output_path` (preserves original resolution)."""
    if not panel_a_path.exists():
        raise FileNotFoundError(f"Panel A image not found: {panel_a_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(panel_a_path, output_path)


def save_panel_b(
    output_path: Path,
    *,
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    """Render the 3x3 confusion-matrix table as a standalone figure."""
    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 14,
            "axes.labelsize": 14,
            "figure.dpi": dpi,
        }
    )

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    _draw_confusion_matrix(ax)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Figure 1A (offset-grid heatmap) and Figure 1B (confusion-matrix table) as separate PNGs."
    )
    parser.add_argument(
        "--panel-a",
        type=Path,
        default=DEFAULT_PANEL_A,
        help="Path to the Panel A source PNG (default: outputs/7jq8/offset_grid_*.png).",
    )
    parser.add_argument(
        "--output-a",
        type=Path,
        default=_REPO_ROOT / "figures/figure_1_a.png",
        help="Output path for Panel A (default: figures/figure_1_a.png).",
    )
    parser.add_argument(
        "--output-b",
        type=Path,
        default=_REPO_ROOT / "figures/figure_1_b.png",
        help="Output path for Panel B (default: figures/figure_1_b.png).",
    )
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument(
        "--fig-width",
        type=float,
        default=8.0,
        help="Figure width (inches) for Panel B.",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=3.0,
        help="Figure height (inches) for Panel B.",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    panel_a = args.panel_a if args.panel_a.is_absolute() else _REPO_ROOT / args.panel_a
    output_a = args.output_a if args.output_a.is_absolute() else _REPO_ROOT / args.output_a
    output_b = args.output_b if args.output_b.is_absolute() else _REPO_ROOT / args.output_b

    save_panel_a(panel_a, output_a)
    save_panel_b(
        output_b,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )

    print(f"Figure 1A written to {output_a.resolve()}")
    print(f"Figure 1B written to {output_b.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
