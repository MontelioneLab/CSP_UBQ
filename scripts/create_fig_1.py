#!/usr/bin/env python3
"""
Create publication-ready Figure 1 panels and a stacked combined figure:

- figure_1_a.png: the chemical-shift offset grid-search heatmap for PDB 7JQ8,
  copied from the pipeline ``outputs/<HOLO>_<apo_bmrb>/offset_grid_*.png`` (resolved
  from ``data/CSP_UBQ.csv`` when the legacy ``outputs/7jq8/`` path is absent).
- figure_1_b.png: a 3x3 confusion-matrix table summarizing the TP/FP/FN/TN
  classification scheme used throughout the study, with colors pulled from
  `scripts/config.py` so they stay consistent with the rest of the project.
- figure_1.png: panel A above panel B, labeled ``A.`` and ``B.``.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import Referencing, classification_colors
    from .csp import _build_param_slug
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import Referencing, classification_colors
    from scripts.csp import _build_param_slug
    from scripts.target_resolution import load_target_rows, resolve_target_rows


_FIG1_HOLO_PDB = "7jq8"


def _default_hn_offset_grid_basename(*, ext: str = "png") -> str:
    cfg = Referencing()
    slug = _build_param_slug(
        h_min=cfg.grid_h_min,
        h_max=cfg.grid_h_max,
        h_step=cfg.grid_h_step,
        n_min=cfg.grid_n_min,
        n_max=cfg.grid_n_max,
        n_step=cfg.grid_n_step,
        cutoff=float(cfg.grid_cutoff),
    )
    return f"offset_grid_{slug}.{ext}"


_FIG1_OFFSET_GRID_BASENAME = _default_hn_offset_grid_basename(ext="png")

DEFAULT_PANEL_A = (
    _REPO_ROOT
    / "outputs"
    / "7jq8"
    / _FIG1_OFFSET_GRID_BASENAME
)


def resolve_fig1_panel_a_from_pipeline(outputs_dir: Path, targets_csv: Path) -> Path | None:
    """Locate Panel A PNG under canonical ``outputs/<HOLO>_<apo_bmrb>/`` for PDB 7JQ8."""
    if not outputs_dir.is_dir() or not targets_csv.is_file():
        return None
    rows = [
        r
        for r in load_target_rows(targets_csv)
        if r.holo_pdb.strip().lower() == _FIG1_HOLO_PDB
    ]
    if not rows:
        return None
    for out_path in resolve_target_rows(rows, outputs_dir, log_warnings=False):
        cand = out_path / _FIG1_OFFSET_GRID_BASENAME
        if cand.is_file():
            return cand
    return None


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
    if panel_a_path.resolve() != output_path.resolve():
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


def compose_two_panel_figure(
    panel_a_path: Path,
    panel_b_path: Path,
    output_image: Path,
    *,
    dpi: int,
) -> None:
    """Stack existing Panel A/B PNGs with ``A.`` / ``B.`` labels."""
    image_a = mpimg.imread(panel_a_path)
    image_b = mpimg.imread(panel_b_path)
    aspect_a = image_a.shape[0] / max(image_a.shape[1], 1)
    aspect_b = image_b.shape[0] / max(image_b.shape[1], 1)
    fig_width = 8.0
    fig_height = fig_width * (aspect_a + aspect_b) + 0.6

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": [aspect_a, aspect_b]},
    )
    for ax, image, label in ((axes[0], image_a, "A."), (axes[1], image_b, "B.")):
        ax.imshow(image)
        ax.set_axis_off()
        ax.text(
            0.01,
            1.02,
            label,
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=20,
            fontweight="bold",
            clip_on=False,
        )

    plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.02, hspace=0.08)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=dpi)
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Figure 1A, Figure 1B, and stacked figure_1.png (A. above B.)."
    )
    parser.add_argument(
        "--panel-a",
        type=Path,
        default=DEFAULT_PANEL_A,
        help=(
            "Path to the Panel A source PNG (default: legacy outputs/7jq8/…; "
            "if missing, resolved from --targets-csv under --outputs-dir)."
        ),
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=_REPO_ROOT / "outputs",
        help="Pipeline outputs root for resolving Panel A (default: outputs).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ.csv",
        help="Study CSV with holo_pdb / apo_bmrb for resolving Panel A (default: data/CSP_UBQ.csv).",
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
    parser.add_argument(
        "--output-combined",
        type=Path,
        default=_REPO_ROOT / "figures/figure_1.png",
        help="Stacked combined figure path (default: figures/figure_1.png).",
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
    outputs_dir = (
        args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    )
    targets_csv = (
        args.targets_csv if args.targets_csv.is_absolute() else _REPO_ROOT / args.targets_csv
    )
    output_a = args.output_a if args.output_a.is_absolute() else _REPO_ROOT / args.output_a
    output_b = args.output_b if args.output_b.is_absolute() else _REPO_ROOT / args.output_b
    output_c = (
        args.output_combined
        if args.output_combined.is_absolute()
        else _REPO_ROOT / args.output_combined
    )

    if not panel_a.exists() and panel_a.resolve() == DEFAULT_PANEL_A.resolve():
        alt = resolve_fig1_panel_a_from_pipeline(outputs_dir, targets_csv)
        if alt is not None:
            panel_a = alt

    save_panel_a(panel_a, output_a)
    save_panel_b(
        output_b,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )

    compose_two_panel_figure(output_a, output_b, output_c, dpi=args.dpi)

    print(f"Figure 1A written to {output_a.resolve()}")
    print(f"Figure 1B written to {output_b.resolve()}")
    print(f"Figure 1 written to {output_c.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
