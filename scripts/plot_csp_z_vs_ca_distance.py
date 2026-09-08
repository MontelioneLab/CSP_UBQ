"""
CSP Z-Score vs Nearest Interchain Distance Scatterplots

Creates confusion-matrix-colored scatterplots of CSP z-score (Y) vs nearest
interchain distance (X) — CA–CA or atom–atom — for a single target or across
all targets under an outputs directory.

Y uses the stored pipeline ``csp_z`` (relative to cleaned mean / cleaned SD after
iterative outlier removal). Point colors use TP/FP/TN/FN derived from a chosen
significance column plus the standard binding-site ground truth.

Usage:
    python scripts/plot_csp_z_vs_ca_distance.py --outputs outputs --all
    python scripts/plot_csp_z_vs_ca_distance.py --outputs outputs --target outputs/1CF4_18251
    python scripts/plot_csp_z_vs_ca_distance.py --outputs outputs --all \\
        --output-combined outputs/csp_z_vs_ca_distance_scatter_cleaned_mean_vs_max05.png \\
        --output-quad outputs/csp_z_vs_distance_scatter_ca_vs_any_atom_2x2.png
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

try:
    from .config import classification_colors, paths
    from .merge_csv import compute_classification
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors, paths
    from scripts.merge_csv import compute_classification

Point = Tuple[float, float, str]  # (distance, csp_z, classification)

_VALID_CLASSES = ("TP", "FP", "TN", "FN")

_CLASS_COLORS = {
    "TP": classification_colors.TP,
    "FP": classification_colors.FP,
    "TN": classification_colors.TN,
    "FN": classification_colors.FN,
}

_LEGEND_LABELS = {
    "TP": "(TP) Sig. CSP in Binding Site",
    "FP": "(FP) Sig. CSP -- Allosteric",
    "TN": "(TN) low CSP -- Allosteric",
    "FN": "(FN) low CSP in Binding Site",
}

DEFAULT_SIGNIFICANT_COLUMN = "significant"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
ANY_ATOM_DISTANCE_COLUMN = "min_any_atom_distance"
CA_DISTANCE_XLABEL = "Nearest interchain CA–CA distance (Å)"
ANY_ATOM_DISTANCE_XLABEL = "Nearest interchain atom–atom distance (Å)"


def collect_points_from_master_csv(
    path: str,
    *,
    significant_column: str = DEFAULT_SIGNIFICANT_COLUMN,
    distance_column: str = ANY_ATOM_DISTANCE_COLUMN,
) -> List[Point]:
    """Extract (distance, csp_z, classification) points from one master_alignment.csv.

    Classification is recomputed from ``significant_column`` and binding-site
    fields via :func:`merge_csv.compute_classification`. ``csp_z`` is taken as
    stored (cleaned-mean z-score).
    """
    points: List[Point] = []
    try:
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            if significant_column not in fieldnames:
                print(
                    f"Warning: {path} missing significance column "
                    f"{significant_column!r}; skipping"
                )
                return points
            if distance_column not in fieldnames:
                print(
                    f"Warning: {path} missing distance column "
                    f"{distance_column!r}; skipping"
                )
                return points
            for row in reader:
                csp_z_str = (row.get("csp_z") or "").strip()
                distance_str = (row.get(distance_column) or "").strip()
                if not csp_z_str or not distance_str:
                    continue
                try:
                    csp_z = float(csp_z_str)
                    distance = float(distance_str)
                except (ValueError, TypeError):
                    continue
                if not (np.isfinite(csp_z) and np.isfinite(distance)):
                    continue

                row_for_cls = dict(row)
                row_for_cls["significant"] = row.get(significant_column, "")
                classification = compute_classification(row_for_cls).strip().upper()
                if classification not in _VALID_CLASSES:
                    continue
                points.append((distance, csp_z, classification))
    except Exception as e:
        print(f"Warning: Error reading {path}: {e}")
    return points


def find_master_alignment_files(outputs_dir: str) -> List[str]:
    """Find all master_alignment.csv files under outputs_dir."""
    pattern = os.path.join(outputs_dir, "**", "master_alignment.csv")
    return sorted(glob.glob(pattern, recursive=True))


def collect_points_from_outputs(
    outputs_dir: str,
    *,
    significant_column: str = DEFAULT_SIGNIFICANT_COLUMN,
    distance_column: str = ANY_ATOM_DISTANCE_COLUMN,
    selected_dir_names: Optional[Set[str]] = None,
) -> List[Point]:
    """Collect points from master_alignment.csv files under outputs_dir.

    If ``selected_dir_names`` is set, only those per-target subdirectories
    (basenames) are included.
    """
    points: List[Point] = []
    for csv_path in find_master_alignment_files(outputs_dir):
        if selected_dir_names is not None:
            parent = os.path.basename(os.path.dirname(os.path.abspath(csv_path)))
            if parent not in selected_dir_names:
                continue
        points.extend(
            collect_points_from_master_csv(
                csv_path,
                significant_column=significant_column,
                distance_column=distance_column,
            )
        )
    return points


def _selected_dir_names(
    outputs_dir: str, targets_csv: Optional[str]
) -> Optional[Set[str]]:
    if not targets_csv:
        return None
    try:
        from .target_resolution import load_target_rows, resolve_target_rows
    except Exception:
        from scripts.target_resolution import load_target_rows, resolve_target_rows

    csv_path = Path(targets_csv)
    rows = load_target_rows(csv_path)
    paths = resolve_target_rows(rows, Path(outputs_dir))
    names = {p.name for p in paths}
    print(f"Restricting scatter to {len(names)} target dirs from {csv_path.name}")
    return names


def _fp_quadrant_count(points: Sequence[Point]) -> int:
    return sum(
        1 for dist, csp_z, cls in points if cls == "FP" and csp_z > 3.0 and dist > 10.0
    )


def _panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.02,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        color="black",
        clip_on=False,
    )


def _shared_ylim(axes, point_groups: Sequence[Sequence[Point]]) -> None:
    y_vals = [p[1] for group in point_groups for p in group]
    if not y_vals:
        return
    y_min, y_max = min(y_vals), max(y_vals)
    pad = 0.05 * (y_max - y_min or 1.0)
    for ax in axes:
        ax.set_ylim(y_min - pad, y_max + pad)


def _draw_csp_z_vs_distance_scatter(
    ax,
    points: Sequence[Point],
    *,
    title: Optional[str] = None,
    xlabel: str = CA_DISTANCE_XLABEL,
    legend_fontsize: float = 9,
    label_fontsize: float = 12,
    title_fontsize: float = 14,
    marker_size: float = 28,
    annotate_fp_quadrant: bool = True,
    fp_text_fontsize: float = 14,
) -> None:
    """Draw CSP z vs distance scatter onto an existing axes."""
    for cls in _VALID_CLASSES:
        xs = [p[0] for p in points if p[2] == cls]
        ys = [p[1] for p in points if p[2] == cls]
        if not xs:
            continue
        ax.scatter(
            xs,
            ys,
            c=_CLASS_COLORS[cls],
            alpha=0.55,
            edgecolors="none",
            s=marker_size,
            label=f"{_LEGEND_LABELS[cls]} ({len(xs)})",
            zorder=3 if cls in ("TP", "FP") else 2,
        )

    ax.axhline(3.0, color="red", linestyle="--", linewidth=1.5, zorder=4)
    ax.axvline(10.0, color="red", linestyle="--", linewidth=1.5, zorder=4)

    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel("CSP z-score", fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
    ax.legend(loc="best", framealpha=0.9, fontsize=legend_fontsize)

    if annotate_fp_quadrant:
        fp_n = _fp_quadrant_count(points)
        ax.text(
            0.97,
            0.97,
            f"{fp_n} FP\n(d > 10 Å, z > 3)",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=fp_text_fontsize,
            fontweight="bold",
            color="red",
            zorder=5,
        )


# Backward-compatible alias used by older imports/tests.
_draw_csp_z_vs_ca_distance_scatter = _draw_csp_z_vs_distance_scatter


def plot_csp_z_vs_ca_distance_scatter(
    points: Sequence[Point],
    out_path: str,
    *,
    title: Optional[str] = None,
    xlabel: str = CA_DISTANCE_XLABEL,
) -> None:
    """Draw and save a single CSP z vs distance scatterplot."""
    if not points:
        print(f"Warning: No valid points to plot; skipping {out_path}")
        return

    fig, ax = plt.subplots(figsize=(10, 8))
    _draw_csp_z_vs_distance_scatter(
        ax,
        points,
        title=title or "CSP Z-Score vs Nearest Interchain Distance",
        xlabel=xlabel,
    )
    plt.tight_layout()

    output_dir = os.path.dirname(out_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved scatterplot to {out_path} ({len(points)} points)")
    print(
        f"FP with CSP z-score > 3 and distance > 10 Å: {_fp_quadrant_count(points)}"
    )


def plot_dual_significance_scatters(
    points_cleaned_mean: Sequence[Point],
    points_max05: Sequence[Point],
    out_path: str,
    *,
    xlabel: str = CA_DISTANCE_XLABEL,
    fig_width: float = 18.0,
    fig_height: float = 8.0,
) -> None:
    """Side-by-side scatters: (a) max(0.05, cleaned mean), (b) cleaned mean."""
    if not points_max05 and not points_cleaned_mean:
        print(f"Warning: No valid points to plot; skipping {out_path}")
        return

    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(fig_width, fig_height),
        constrained_layout=True,
    )

    _draw_csp_z_vs_distance_scatter(
        ax_a,
        points_max05,
        title="max(0.05 ppm, cleaned mean)",
        xlabel=xlabel,
        legend_fontsize=8,
        marker_size=22,
    )
    _draw_csp_z_vs_distance_scatter(
        ax_b,
        points_cleaned_mean,
        title="Cleaned mean",
        xlabel=xlabel,
        legend_fontsize=8,
        marker_size=22,
    )

    _panel_label(ax_a, "a")
    _panel_label(ax_b, "b")
    _shared_ylim((ax_a, ax_b), (points_max05, points_cleaned_mean))

    output_dir = os.path.dirname(out_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    print(
        f"Saved combined scatterplot to {out_path} "
        f"(a={len(points_max05)} points, b={len(points_cleaned_mean)} points)"
    )
    print(
        f"Panel a FP (z>3, d>10 Å): {_fp_quadrant_count(points_max05)}; "
        f"panel b: {_fp_quadrant_count(points_cleaned_mean)}"
    )


def plot_quad_distance_scatters(
    ca_max05: Sequence[Point],
    ca_cleaned: Sequence[Point],
    atom_max05: Sequence[Point],
    atom_cleaned: Sequence[Point],
    out_path: str,
    *,
    fig_width: float = 18.0,
    fig_height: float = 14.0,
) -> None:
    """2×2: a/b CA–CA (max05 / cleaned); c/d atom–atom (max05 / cleaned)."""
    groups = (ca_max05, ca_cleaned, atom_max05, atom_cleaned)
    if not any(groups):
        print(f"Warning: No valid points to plot; skipping {out_path}")
        return

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(fig_width, fig_height),
        constrained_layout=True,
    )
    ax_a, ax_b = axes[0, 0], axes[0, 1]
    ax_c, ax_d = axes[1, 0], axes[1, 1]

    panel_specs = (
        (ax_a, ca_max05, "a", "CA–CA · max(0.05 ppm, cleaned mean)", CA_DISTANCE_XLABEL),
        (ax_b, ca_cleaned, "b", "CA–CA · cleaned mean", CA_DISTANCE_XLABEL),
        (ax_c, atom_max05, "c", "Atom–atom · max(0.05 ppm, cleaned mean)", ANY_ATOM_DISTANCE_XLABEL),
        (ax_d, atom_cleaned, "d", "Atom–atom · cleaned mean", ANY_ATOM_DISTANCE_XLABEL),
    )
    for ax, pts, label, title, xlabel in panel_specs:
        _draw_csp_z_vs_distance_scatter(
            ax,
            pts,
            title=title,
            xlabel=xlabel,
            legend_fontsize=7,
            marker_size=18,
            fp_text_fontsize=12,
        )
        _panel_label(ax, label)

    _shared_ylim((ax_a, ax_b, ax_c, ax_d), groups)

    output_dir = os.path.dirname(out_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    print(
        f"Saved 2×2 scatterplot to {out_path} "
        f"(a={len(ca_max05)}, b={len(ca_cleaned)}, "
        f"c={len(atom_max05)}, d={len(atom_cleaned)} points)"
    )
    print(
        "FP (z>3, d>10 Å): "
        f"a={_fp_quadrant_count(ca_max05)}, "
        f"b={_fp_quadrant_count(ca_cleaned)}, "
        f"c={_fp_quadrant_count(atom_max05)}, "
        f"d={_fp_quadrant_count(atom_cleaned)}"
    )


def plot_per_target(
    tgt_dir: str,
    out_path: Optional[str] = None,
    *,
    significant_column: str = DEFAULT_SIGNIFICANT_COLUMN,
    distance_column: str = ANY_ATOM_DISTANCE_COLUMN,
) -> Optional[str]:
    """Create the scatterplot for a single target directory."""
    master_csv = os.path.join(tgt_dir, "master_alignment.csv")
    if not os.path.isfile(master_csv):
        print(f"Warning: master_alignment.csv not found in {tgt_dir}; skipping scatterplot")
        return None

    if out_path is None:
        out_path = os.path.join(tgt_dir, "csp_z_vs_ca_distance_scatter.png")

    points = collect_points_from_master_csv(
        master_csv,
        significant_column=significant_column,
        distance_column=distance_column,
    )
    target_name = os.path.basename(os.path.normpath(tgt_dir))
    xlabel = (
        ANY_ATOM_DISTANCE_XLABEL
        if distance_column == ANY_ATOM_DISTANCE_COLUMN
        else CA_DISTANCE_XLABEL
    )
    plot_csp_z_vs_ca_distance_scatter(
        points,
        out_path,
        title=f"{target_name} CSP Z-Score vs Distance",
        xlabel=xlabel,
    )
    return out_path if points else None


def plot_dataset(
    outputs_dir: str,
    out_path: Optional[str] = None,
    *,
    significant_column: str = DEFAULT_SIGNIFICANT_COLUMN,
    distance_column: str = ANY_ATOM_DISTANCE_COLUMN,
    selected_dir_names: Optional[Set[str]] = None,
    title: Optional[str] = None,
) -> Optional[str]:
    """Create the scatterplot for residues across targets under outputs_dir."""
    if out_path is None:
        out_path = os.path.join(outputs_dir, "csp_z_vs_ca_distance_scatter.png")

    points = collect_points_from_outputs(
        outputs_dir,
        significant_column=significant_column,
        distance_column=distance_column,
        selected_dir_names=selected_dir_names,
    )
    xlabel = (
        ANY_ATOM_DISTANCE_XLABEL
        if distance_column == ANY_ATOM_DISTANCE_COLUMN
        else CA_DISTANCE_XLABEL
    )
    n_tgt = len(selected_dir_names) if selected_dir_names is not None else None
    default_title = (
        f"pH/temp matched (n={n_tgt}): CSP Z-Score vs Nearest Interchain Distance"
        if n_tgt is not None
        else "All Targets: CSP Z-Score vs Nearest Interchain Distance"
    )
    plot_csp_z_vs_ca_distance_scatter(
        points,
        out_path,
        title=title or default_title,
        xlabel=xlabel,
    )
    return out_path if points else None


def main() -> None:
    default_outputs = paths.outputs_dir if hasattr(paths, "outputs_dir") else "outputs"
    parser = argparse.ArgumentParser(
        description="Scatterplot of CSP z-score vs nearest interchain distance, "
        "colored by confusion-matrix classification"
    )
    parser.add_argument(
        "--outputs",
        default=default_outputs,
        help="Path to outputs directory (default: outputs/)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--target",
        metavar="DIR",
        help="Single target directory containing master_alignment.csv",
    )
    mode.add_argument(
        "--all",
        action="store_true",
        help="Plot residues across master_alignment.csv files under --outputs",
    )
    parser.add_argument(
        "--targets-csv",
        default=None,
        help="With --all: restrict to target rows in this CSV (holo_pdb/apo_bmrb).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PNG path (default depends on --target / --all)",
    )
    parser.add_argument(
        "--significant-column",
        default=DEFAULT_SIGNIFICANT_COLUMN,
        help=(
            "Master-alignment boolean column used for significance when deriving "
            "TP/FP/TN/FN (default: significant = max(cleaned mean, 0.05); "
            "use significant_sigma_0 for unfloored cleaned mean)."
        ),
    )
    parser.add_argument(
        "--distance-column",
        default=ANY_ATOM_DISTANCE_COLUMN,
        choices=(CA_DISTANCE_COLUMN, ANY_ATOM_DISTANCE_COLUMN),
        help=(
            f"Distance column for X-axis (default: {ANY_ATOM_DISTANCE_COLUMN}; "
            f"use {CA_DISTANCE_COLUMN} for nearest interchain CA–CA)."
        ),
    )
    parser.add_argument(
        "--output-combined",
        default=None,
        help=(
            "With --all: write CA–CA dual figure (a=max05, b=cleaned mean) and the "
            "matching atom–atom dual companion under --outputs."
        ),
    )
    parser.add_argument(
        "--output-quad",
        default=None,
        help=(
            "With --all: write 2×2 figure a/b=CA–CA (max05/cleaned), "
            "c/d=atom–atom (max05/cleaned)."
        ),
    )

    args = parser.parse_args()
    sig_col = args.significant_column
    selected = _selected_dir_names(args.outputs, args.targets_csv)

    if args.target:
        result = plot_per_target(
            args.target,
            out_path=args.output,
            significant_column=sig_col,
            distance_column=args.distance_column,
        )
        if result is None:
            print("Error: No scatterplot produced (missing data or no valid points)")
            sys.exit(1)
        print("Done!")
        return

    if args.output_combined is not None or args.output_quad is not None:
        ca_max05 = collect_points_from_outputs(
            args.outputs,
            significant_column="significant",
            distance_column=CA_DISTANCE_COLUMN,
            selected_dir_names=selected,
        )
        ca_cleaned = collect_points_from_outputs(
            args.outputs,
            significant_column="significant_sigma_0",
            distance_column=CA_DISTANCE_COLUMN,
            selected_dir_names=selected,
        )
        atom_max05 = collect_points_from_outputs(
            args.outputs,
            significant_column="significant",
            distance_column=ANY_ATOM_DISTANCE_COLUMN,
            selected_dir_names=selected,
        )
        atom_cleaned = collect_points_from_outputs(
            args.outputs,
            significant_column="significant_sigma_0",
            distance_column=ANY_ATOM_DISTANCE_COLUMN,
            selected_dir_names=selected,
        )

        out_ca_max05 = os.path.join(
            args.outputs, "csp_z_vs_ca_distance_scatter_max05.png"
        )
        out_ca_cleaned = os.path.join(
            args.outputs, "csp_z_vs_ca_distance_scatter_cleaned_mean.png"
        )
        out_atom_max05 = os.path.join(
            args.outputs, "csp_z_vs_any_atom_distance_scatter_max05.png"
        )
        out_atom_cleaned = os.path.join(
            args.outputs, "csp_z_vs_any_atom_distance_scatter_cleaned_mean.png"
        )
        out_ca_dual = args.output_combined or os.path.join(
            args.outputs, "csp_z_vs_ca_distance_scatter_cleaned_mean_vs_max05.png"
        )
        out_atom_dual = os.path.join(
            args.outputs,
            "csp_z_vs_any_atom_distance_scatter_cleaned_mean_vs_max05.png",
        )
        out_quad = args.output_quad or os.path.join(
            args.outputs, "csp_z_vs_distance_scatter_ca_vs_any_atom_2x2.png"
        )

        if args.output:
            print(
                "Note: with --output-combined/--output-quad, writing standard "
                "paths under --outputs (ignoring --output)."
            )

        # Individual panels
        plot_csp_z_vs_ca_distance_scatter(
            ca_max05,
            out_ca_max05,
            title="All Targets: CSP Z-Score vs CA–CA Distance "
            "(max(0.05 ppm, cleaned mean))",
            xlabel=CA_DISTANCE_XLABEL,
        )
        plot_csp_z_vs_ca_distance_scatter(
            ca_cleaned,
            out_ca_cleaned,
            title="All Targets: CSP Z-Score vs CA–CA Distance (cleaned mean)",
            xlabel=CA_DISTANCE_XLABEL,
        )
        plot_csp_z_vs_ca_distance_scatter(
            atom_max05,
            out_atom_max05,
            title="All Targets: CSP Z-Score vs Atom–Atom Distance "
            "(max(0.05 ppm, cleaned mean))",
            xlabel=ANY_ATOM_DISTANCE_XLABEL,
        )
        plot_csp_z_vs_ca_distance_scatter(
            atom_cleaned,
            out_atom_cleaned,
            title="All Targets: CSP Z-Score vs Atom–Atom Distance (cleaned mean)",
            xlabel=ANY_ATOM_DISTANCE_XLABEL,
        )

        # Dual companions
        plot_dual_significance_scatters(
            ca_cleaned,
            ca_max05,
            out_ca_dual,
            xlabel=CA_DISTANCE_XLABEL,
        )
        plot_dual_significance_scatters(
            atom_cleaned,
            atom_max05,
            out_atom_dual,
            xlabel=ANY_ATOM_DISTANCE_XLABEL,
        )

        # Four-panel figure
        plot_quad_distance_scatters(
            ca_max05,
            ca_cleaned,
            atom_max05,
            atom_cleaned,
            out_quad,
        )
        print("Done!")
        return

    result = plot_dataset(
        args.outputs,
        out_path=args.output,
        significant_column=sig_col,
        distance_column=args.distance_column,
        selected_dir_names=selected,
    )
    if result is None:
        print("Error: No scatterplot produced (missing data or no valid points)")
        sys.exit(1)
    print("Done!")


if __name__ == "__main__":
    main()
