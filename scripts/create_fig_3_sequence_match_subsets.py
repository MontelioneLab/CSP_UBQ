#!/usr/bin/env python3
"""
Create a 2×3 Figure 3-style grid comparing sequence-match cohorts:

  columns: Ideal | Non-ideal | Combined (union of the two)
  rows:    Fig 3a (TP/FP) | Fig 3b (TN/FP/FN/TP)

Default CSVs:
  data/CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv
  data/CSP_UBQ_ph0.5_temp5C_non_ideal_sequence_match.csv

Default output: figures/figure_3_sequence_match_subsets.png
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Patch

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import classification_colors
    from .create_fig_3 import (
        FIGURE_3_DISTANCE_XLABEL,
        _get_bins,
        _resolve_selected_dirs,
        _set_plot_style,
        collect_distance_categories,
    )
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors
    from scripts.create_fig_3 import (
        FIGURE_3_DISTANCE_XLABEL,
        _get_bins,
        _resolve_selected_dirs,
        _set_plot_style,
        collect_distance_categories,
    )


DistanceCategories = Tuple[List[float], List[float], List[float], List[float]]

DEFAULT_IDEAL_CSV = _REPO_ROOT / "data" / "CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv"
DEFAULT_NON_IDEAL_CSV = (
    _REPO_ROOT / "data" / "CSP_UBQ_ph0.5_temp5C_non_ideal_sequence_match.csv"
)


def _read_csv_rows(path: Path) -> Tuple[List[str], List[dict]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return fieldnames, rows


def _write_union_csv(ideal_csv: Path, non_ideal_csv: Path, out_path: Path) -> int:
    """Write unique (apo_bmrb, holo_bmrb, holo_pdb) rows from both CSVs; return n rows."""
    fields_a, rows_a = _read_csv_rows(ideal_csv)
    fields_b, rows_b = _read_csv_rows(non_ideal_csv)
    if fields_a != fields_b:
        raise ValueError(
            f"CSV fieldnames differ:\n  {ideal_csv}: {fields_a}\n  {non_ideal_csv}: {fields_b}"
        )
    seen: Set[Tuple[str, str, str]] = set()
    merged: List[dict] = []
    for row in rows_a + rows_b:
        key = (
            row.get("apo_bmrb", "").strip(),
            row.get("holo_bmrb", "").strip(),
            row.get("holo_pdb", "").strip().upper(),
        )
        if key in seen:
            continue
        seen.add(key)
        merged.append(row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields_a)
        writer.writeheader()
        writer.writerows(merged)
    return len(merged)


def _plot_3a_on_ax(
    ax: Axes,
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    bins: Sequence[float],
    x_max: float,
) -> None:
    ax.hist(
        [list(tp_distances), list(fp_distances)],
        bins=bins,
        stacked=True,
        color=[classification_colors.TP, classification_colors.FP],
        edgecolor="black",
        linewidth=0.6,
    )
    ax.set_xlim(0, x_max)

    n_tp = len(tp_distances)
    n_fp = len(fp_distances)
    denom = n_tp + n_fp
    if denom > 0:
        fp_pct = 100.0 * n_fp / denom
        ax.text(
            0.50,
            0.55,
            f"{fp_pct:.1f}% FP",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color="red",
            zorder=10,
        )


def _plot_3b_on_ax(
    ax: Axes,
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    fn_distances: Sequence[float],
    tn_distances: Sequence[float],
    bins: Sequence[float],
    x_max: float,
) -> None:
    ax.hist(
        [
            list(tn_distances),
            list(fp_distances),
            list(fn_distances),
            list(tp_distances),
        ],
        bins=bins,
        stacked=True,
        color=[
            classification_colors.TN,
            classification_colors.FP,
            classification_colors.FN,
            classification_colors.TP,
        ],
        edgecolor="black",
        linewidth=0.6,
    )
    ax.set_xlim(0, x_max)


def plot_subset_grid(
    cohorts: Sequence[Tuple[str, DistanceCategories, int]],
    output_path: Path,
    *,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    """Draw 2 rows (3a, 3b) × N columns (one per cohort)."""
    n_cols = len(cohorts)
    _set_plot_style(dpi)
    fig, axes = plt.subplots(
        2,
        n_cols,
        figsize=(fig_width, fig_height),
        sharex="col",
    )
    if n_cols == 1:
        axes = axes.reshape(2, 1)

    all_distances: List[float] = []
    for _, (tp, fp, fn, tn), _ in cohorts:
        all_distances.extend(tp)
        all_distances.extend(fp)
        all_distances.extend(fn)
        all_distances.extend(tn)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)

    for col_idx, (title, (tp, fp, fn, tn), n_targets) in enumerate(cohorts):
        ax_a = axes[0, col_idx]
        ax_b = axes[1, col_idx]
        _plot_3a_on_ax(ax_a, tp, fp, bins, x_max)
        _plot_3b_on_ax(ax_b, tp, fp, fn, tn, bins, x_max)

        ax_a.set_title(f"{title} (n={n_targets})", fontsize=13)
        if col_idx == 0:
            ax_a.set_ylabel("Number of Residues")
            ax_b.set_ylabel("Number of Residues")
        ax_b.set_xlabel(FIGURE_3_DISTANCE_XLABEL)

        ax_a.text(
            -0.08,
            1.02,
            "a",
            transform=ax_a.transAxes,
            ha="left",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            clip_on=False,
        )
        ax_b.text(
            -0.08,
            1.02,
            "b",
            transform=ax_b.transAxes,
            ha="left",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            clip_on=False,
        )

    legend_handles = [
        Patch(facecolor=classification_colors.TP, edgecolor="black", label="TP"),
        Patch(facecolor=classification_colors.FP, edgecolor="black", label="FP"),
        Patch(facecolor=classification_colors.FN, edgecolor="black", label="FN"),
        Patch(facecolor=classification_colors.TN, edgecolor="black", label="TN"),
    ]
    fig.subplots_adjust(left=0.06, right=0.88, top=0.92, bottom=0.08, wspace=0.22, hspace=0.28)
    fig.legend(
        handles=legend_handles,
        title="Classification",
        loc="upper left",
        bbox_to_anchor=(0.90, 0.92),
        frameon=True,
        fontsize=11,
        title_fontsize=12,
        borderaxespad=0.0,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create 2×3 Figure 3 grid for ideal / non-ideal / combined sequence-match subsets."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory of per-target outputs (default: outputs).",
    )
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--ideal-csv",
        type=Path,
        default=DEFAULT_IDEAL_CSV,
        help="Ideal sequence-match CSV (default: data/CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv).",
    )
    parser.add_argument(
        "--non-ideal-csv",
        type=Path,
        default=DEFAULT_NON_IDEAL_CSV,
        help="Non-ideal sequence-match CSV (default: data/CSP_UBQ_ph0.5_temp5C_non_ideal_sequence_match.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: figures/figure_3_sequence_match_subsets.png).",
    )
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--max-distance", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--fig-width", type=float, default=18.0)
    parser.add_argument("--fig-height", type=float, default=10.0)
    return parser.parse_args(list(argv))


def _abs(path: Path) -> Path:
    return path if path.is_absolute() else _REPO_ROOT / path


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    outputs_dir = _abs(args.outputs_dir)
    figures_dir = _abs(args.figures_dir)
    ideal_csv = _abs(args.ideal_csv)
    non_ideal_csv = _abs(args.non_ideal_csv)
    output_path = _abs(args.output or (figures_dir / "figure_3_sequence_match_subsets.png"))

    ideal_dirs = _resolve_selected_dirs(outputs_dir, ideal_csv, None) or set()
    non_ideal_dirs = _resolve_selected_dirs(outputs_dir, non_ideal_csv, None) or set()

    with tempfile.TemporaryDirectory(prefix="fig3_seqmatch_") as tmp:
        union_csv = Path(tmp) / "combined_sequence_match.csv"
        n_combined_rows = _write_union_csv(ideal_csv, non_ideal_csv, union_csv)
        combined_dirs = _resolve_selected_dirs(outputs_dir, union_csv, None) or set()

        print(f"[fig3-seqmatch] Ideal dirs: {len(ideal_dirs)} (csv rows from {ideal_csv.name})")
        print(
            f"[fig3-seqmatch] Non-ideal dirs: {len(non_ideal_dirs)} "
            f"(csv rows from {non_ideal_csv.name})"
        )
        print(
            f"[fig3-seqmatch] Combined dirs: {len(combined_dirs)} "
            f"(union csv rows={n_combined_rows})"
        )

        ideal_cats = collect_distance_categories(outputs_dir, ideal_dirs)
        non_ideal_cats = collect_distance_categories(outputs_dir, non_ideal_dirs)
        combined_cats = collect_distance_categories(outputs_dir, combined_dirs)

    cohorts: List[Tuple[str, DistanceCategories, int]] = [
        ("Ideal", ideal_cats, len(ideal_dirs)),
        ("Non-ideal", non_ideal_cats, len(non_ideal_dirs)),
        ("Combined", combined_cats, len(combined_dirs)),
    ]

    for title, (tp, fp, fn, tn), n in cohorts:
        print(
            f"[fig3-seqmatch] {title}: targets={n} "
            f"TP={len(tp)} FP={len(fp)} FN={len(fn)} TN={len(tn)}"
        )

    plot_subset_grid(
        cohorts,
        output_path,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )
    print(f"Figure written to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
