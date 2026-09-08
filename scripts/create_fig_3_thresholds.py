#!/usr/bin/env python3
"""
Create a multi-panel Figure 3 comparing alternate CSP significance thresholds.

Layout is a 6×4 grid: each threshold occupies an adjacent 3a|3b pair in
columns 1–2 or 3–4, with similar thresholds placed side-by-side in a row.

Thresholds (columns in master_alignment.csv):
  - max(0.05 ppm, cleaned mean)       → significant (primary)
  - cleaned mean + 1 SD               → significant_1sd
  - cleaned mean + 2 SD               → significant_2sd
  - top 10%                           → significant_top_10_percentile
  - top 5%                            → significant_top_5_percentile
  - 0.03 / 0.05 / 0.10 ppm            → significant_03/05/10_ppm
  - raw mean (no outlier removal)     → significant_raw_mean
  - raw mean + 1 SD                   → significant_raw_1sd
  - raw mean + 2 SD                   → significant_raw_2sd
  - cleaned mean (no floor)           → significant_sigma_0

By default, targets are restricted to holo_pdb IDs in CSP_UBQ_ph0.5_temp5C.csv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import classification_colors
    from .merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors
    from scripts.merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from scripts.target_resolution import load_target_rows, resolve_target_rows


CA_DISTANCE_COLUMN = "min_ca_distance_distance"
CSP_COLUMN = "csp_A"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)

FIGURE_3_DISTANCE_XLABEL = r"Minimum $C_\alpha$ distance (Å)"

# Each entry is one threshold (3a+3b pair occupies two adjacent columns).
# Grid layout is 6 rows × 2 threshold-pairs (4 columns total):
#   cols 1–2 = left pair, cols 3–4 = right pair.
# Similar thresholds are placed side-by-side within a row.
THRESHOLD_GRID: Sequence[Sequence[Tuple[str, str]]] = (
    # Row 0: primary floored cutoff + sigma ladder
    (
        ("significant", "max(0.05 ppm, cleaned mean)"),
        ("significant_1sd", "Cleaned mean + 1 SD"),
    ),
    # Row 1: remaining sigma + rank
    (
        ("significant_2sd", "Cleaned mean + 2 SD"),
        ("significant_top_10_percentile", "Top 10%"),
    ),
    # Row 2: rank + fixed ppm
    (
        ("significant_top_5_percentile", "Top 5%"),
        ("significant_03_ppm", "0.03 ppm"),
    ),
    # Row 3: fixed ppm
    (
        ("significant_05_ppm", "0.05 ppm"),
        ("significant_10_ppm", "0.10 ppm"),
    ),
    # Row 4: raw mean ladder (no outlier removal)
    (
        ("significant_raw_mean", "Raw mean"),
        ("significant_raw_1sd", "Raw mean + 1 SD"),
    ),
    # Row 5: remaining raw + unfloored cleaned mean
    (
        ("significant_raw_2sd", "Raw mean + 2 SD"),
        ("significant_sigma_0", "Cleaned mean"),
    ),
)

THRESHOLD_SPECS: Sequence[Tuple[str, str]] = tuple(
    spec for row in THRESHOLD_GRID for spec in row
)

DistanceCategories = Tuple[List[float], List[float], List[float], List[float]]


def _as_bool(value: object) -> bool:
    """Parse CSV / pandas bool-like cells; missing/invalid → False (for predictor flags)."""
    parsed = parse_optional_bool(value)
    return bool(parsed)


def _resolve_selected_dirs(
    outputs_dir: Path,
    targets_csv: Optional[Path],
    targets_str: Optional[str],
) -> Optional[Set[str]]:
    if targets_csv is None and not targets_str:
        return None
    extras: List[str] = []
    if targets_str:
        extras = [t for t in targets_str.split(",") if t.strip()]
    rows = load_target_rows(targets_csv, extra_holo_pdbs=extras)
    if not rows:
        return set()
    paths = resolve_target_rows(rows, outputs_dir)
    return {p.name for p in paths}


def _get_bins(
    data: Sequence[float],
    bin_width: float,
    max_distance: Optional[float],
) -> Tuple[List[float], float]:
    if not data:
        raise ValueError("No CA-distance data found for selected targets.")
    data_max = max(data) if max_distance is None else max_distance
    n_bins = max(1, int((data_max + bin_width) / bin_width))
    bins = [i * bin_width for i in range(n_bins + 1)]
    return bins, data_max


def collect_distance_categories(
    outputs_dir: Path,
    significant_column: str,
    selected_dir_names: Optional[Set[str]] = None,
) -> DistanceCategories:
    """Collect (TP, FP, FN, TN) min-CA distances for one significance mask."""
    tp_distances: List[float] = []
    fp_distances: List[float] = []
    fn_distances: List[float] = []
    tn_distances: List[float] = []

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if selected_dir_names is not None and target_name not in selected_dir_names:
            continue

        df = pd.read_csv(alignment_path)
        required_cols = {significant_column, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"[WARN] Skipping {alignment_path}: missing columns {missing}")
            continue

        # Exclude residues without a recorded CSP from TN/FN/TP/FP
        df = filter_recorded_csp_dataframe(
            df, csp_column=CSP_COLUMN, significant_column=significant_column
        )
        if df.empty:
            continue

        is_significant = df[significant_column].map(_as_bool)
        is_binding = pd.DataFrame({c: df[c].map(_as_bool) for c in PREDICTOR_COLUMNS}).any(axis=1)
        distances = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")
        valid = distances.notna()
        if not valid.any():
            continue

        sig = is_significant[valid]
        bind = is_binding[valid]
        dist = distances[valid]

        tp_distances.extend(dist[sig & bind].astype(float).tolist())
        fp_distances.extend(dist[sig & ~bind].astype(float).tolist())
        fn_distances.extend(dist[~sig & bind].astype(float).tolist())
        tn_distances.extend(dist[~sig & ~bind].astype(float).tolist())

    return tp_distances, fp_distances, fn_distances, tn_distances


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
            0.50,
            f"{fp_pct:.1f}%",
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


def _annotation_box(ax: Axes, lines: Sequence[str]) -> None:
    ax.text(
        0.98,
        0.95,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "alpha": 0.85,
            "edgecolor": "0.7",
            "linewidth": 0.5,
        },
    )


def plot_threshold_grid(
    categories_by_threshold: Dict[str, DistanceCategories],
    output_path: Path,
    *,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    panel_width: float,
    panel_height: float,
) -> None:
    n_rows = len(THRESHOLD_GRID)
    n_cols = 4  # two thresholds × (3a, 3b) per row
    # Leave right margin for the figure-level legend.
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(panel_width * n_cols + 2.2, panel_height * n_rows),
        sharex="col",
    )
    if n_rows == 1:
        axes = axes.reshape(1, n_cols)

    # Shared bins / x-limits across all panels (fair visual comparison).
    all_distances: List[float] = []
    for tp, fp, fn, tn in categories_by_threshold.values():
        all_distances.extend(tp)
        all_distances.extend(fp)
        all_distances.extend(fn)
        all_distances.extend(tn)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)

    for row_idx, row_specs in enumerate(THRESHOLD_GRID):
        for pair_idx, (column, title) in enumerate(row_specs):
            col_a = pair_idx * 2
            col_b = col_a + 1
            tp, fp, fn, tn = categories_by_threshold[column]
            ax_a = axes[row_idx, col_a]
            ax_b = axes[row_idx, col_b]

            _plot_3a_on_ax(ax_a, tp, fp, bins, x_max)
            _plot_3b_on_ax(ax_b, tp, fp, fn, tn, bins, x_max)

            if pair_idx == 0:
                ax_a.set_ylabel("Number of Residues")

            ax_a.set_title(title, fontsize=11)
            ax_b.set_title(title, fontsize=11)

            _annotation_box(ax_a, [f"TP = {len(tp)}", f"FP = {len(fp)}"])
            _annotation_box(
                ax_b,
                [
                    f"TP = {len(tp)}",
                    f"FP = {len(fp)}",
                    f"FN = {len(fn)}",
                    f"TN = {len(tn)}",
                ],
            )

    for col_idx in range(n_cols):
        axes[-1, col_idx].set_xlabel(FIGURE_3_DISTANCE_XLABEL)

    legend_handles = [
        Patch(facecolor=classification_colors.TP, edgecolor="black", label="TP"),
        Patch(facecolor=classification_colors.FP, edgecolor="black", label="FP"),
        Patch(facecolor=classification_colors.FN, edgecolor="black", label="FN"),
        Patch(facecolor=classification_colors.TN, edgecolor="black", label="TN"),
    ]
    # Reserve right-side space so the legend sits clear of all axes.
    fig.subplots_adjust(left=0.05, right=0.86, top=0.97, bottom=0.05, wspace=0.28, hspace=0.35)
    fig.legend(
        handles=legend_handles,
        title="Classification",
        loc="upper left",
        bbox_to_anchor=(0.88, 0.97),
        frameon=True,
        fontsize=10,
        title_fontsize=11,
        borderaxespad=0.0,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create multi-panel Figure 3 for alternate CSP significance thresholds."
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: figures/figure_3_thresholds.png)",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ_ph0.5_temp5C.csv",
        help="CSV with 'holo_pdb' column (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument("--targets", type=str, help="Optional comma-separated holo_pdb list.")
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--max-distance", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--panel-width", type=float, default=4.0)
    parser.add_argument("--panel-height", type=float, default=3.0)
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    targets_csv = args.targets_csv
    if targets_csv is not None and not targets_csv.is_absolute():
        targets_csv = _REPO_ROOT / targets_csv
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    figures_dir = args.figures_dir if args.figures_dir.is_absolute() else _REPO_ROOT / args.figures_dir
    selected_dirs = _resolve_selected_dirs(outputs_dir, targets_csv, args.targets)

    categories_by_threshold: Dict[str, DistanceCategories] = {}
    for column, title in THRESHOLD_SPECS:
        print(f"[fig3-thresholds] Collecting distances for {title} ({column})")
        categories_by_threshold[column] = collect_distance_categories(
            outputs_dir,
            column,
            selected_dirs,
        )
        tp, fp, fn, tn = categories_by_threshold[column]
        print(f"  TP={len(tp)} FP={len(fp)} FN={len(fn)} TN={len(tn)}")

    output_path = args.output or (figures_dir / "figure_3_thresholds.png")
    if not output_path.is_absolute():
        output_path = _REPO_ROOT / output_path

    plot_threshold_grid(
        categories_by_threshold,
        output_path,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
        dpi=args.dpi,
        panel_width=args.panel_width,
        panel_height=args.panel_height,
    )
    print(f"Figure written to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
