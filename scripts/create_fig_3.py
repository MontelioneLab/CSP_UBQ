#!/usr/bin/env python3
"""
Create publication-ready Figure 3 panels from `outputs/*/master_alignment.csv`:

- figure_3_a.png: stacked histogram of significant residues (TP vs FP only)
- figure_3_b.png: stacked histogram of full confusion matrix (TN, FP, FN, TP)
- figure_3.png: panel a above panel b, with labels and FP/(TP+FP) on a

By default, targets are restricted to holo_pdb IDs listed in CSP_UBQ_ph0.5_temp5C.csv
(same buffer similarity filter as that file). Override with --targets-csv or --targets.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set

import matplotlib.pyplot as plt
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import classification_colors
    from .merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors
    from scripts.merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from scripts.target_resolution import load_target_rows, resolve_target_rows


SIGNIFICANT_COLUMN = "significant"
CSP_COLUMN = "csp_A"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)

# Matplotlib mathtext: C with α subscript (not plaintext "CA").
FIGURE_3_DISTANCE_XLABEL = r"Minimum $C_\alpha$ distance (Å)"


def _as_bool(value: object) -> bool:
    """Parse CSV / pandas bool-like cells; missing/invalid → False (for predictor flags)."""
    parsed = parse_optional_bool(value)
    return bool(parsed)


def _resolve_selected_dirs(
    outputs_dir: Path,
    targets_csv: Optional[Path],
    targets_str: Optional[str],
) -> Optional[Set[str]]:
    """Resolve targets-CSV rows + comma-separated --targets to a set of outputs/<dir> basenames.

    Each row is matched to an outputs/ subdir by congruent ``apo_bmrb`` and
    ``holo_bmrb`` (delegated to :mod:`scripts.target_resolution`); the first
    matching dir wins when several share the same BMRB pair.
    """
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


def _get_bins(data: Sequence[float], bin_width: float, max_distance: Optional[float]) -> tuple[List[float], float]:
    if not data:
        raise ValueError("No CA-distance data found for selected targets.")
    data_max = max(data) if max_distance is None else max_distance
    n_bins = max(1, int((data_max + bin_width) / bin_width))
    bins = [i * bin_width for i in range(n_bins + 1)]
    return bins, data_max


def collect_distance_categories(
    outputs_dir: Path,
    selected_dir_names: Optional[Set[str]] = None,
) -> tuple[List[float], List[float], List[float], List[float]]:
    """Collect (TP, FP, FN, TN) min-CA distances across selected outputs/<dir> targets.

    ``selected_dir_names`` is the set of resolved outputs basenames (e.g.
    ``{"2KWI_2", "2L29_1", ...}``) returned by
    :func:`scripts.target_resolution.resolve_target_rows`. Pass ``None`` to
    include every target subdirectory.
    """
    tp_distances: List[float] = []
    fp_distances: List[float] = []
    fn_distances: List[float] = []
    tn_distances: List[float] = []

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if selected_dir_names is not None and target_name not in selected_dir_names:
            continue

        df = pd.read_csv(alignment_path)
        required_cols = {SIGNIFICANT_COLUMN, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"[WARN] Skipping {alignment_path}: missing columns {missing}")
            continue

        # Exclude residues without a recorded CSP from TN/FN/TP/FP
        df = filter_recorded_csp_dataframe(
            df, csp_column=CSP_COLUMN, significant_column=SIGNIFICANT_COLUMN
        )
        if df.empty:
            continue

        is_significant = df[SIGNIFICANT_COLUMN].map(_as_bool)
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


def _set_plot_style(dpi: int) -> None:
    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 13,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "figure.dpi": dpi,
        }
    )


def plot_figure_3a(
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    output_path: Path,
    *,
    title: str,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    all_distances = list(tp_distances) + list(fp_distances)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)
    _set_plot_style(dpi)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(
        [tp_distances, fp_distances],
        bins=bins,
        stacked=True,
        color=[classification_colors.TP, classification_colors.FP],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TP) CSP --  in binding site ({len(tp_distances)})",
            f"(FP) CSP --  not in binding site ({len(fp_distances)})",
        ],
    )
    ax = plt.gca()
    plt.xlabel(FIGURE_3_DISTANCE_XLABEL)
    plt.ylabel("Number of Residues")
    if title.strip():
        plt.title(title)
    ax.tick_params(axis="both", labelsize=14)
    plt.legend( frameon=True, fontsize=14, title_fontsize=14)
    plt.xlim(0, x_max)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_figure_3b(
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    fn_distances: Sequence[float],
    tn_distances: Sequence[float],
    output_path: Path,
    *,
    title: str,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    all_distances = list(tp_distances) + list(fp_distances) + list(fn_distances) + list(tn_distances)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)
    _set_plot_style(dpi)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(
        [tn_distances, fp_distances, fn_distances, tp_distances],
        bins=bins,
        stacked=True,
        color=[
            classification_colors.TN,
            classification_colors.FP,
            classification_colors.FN,
            classification_colors.TP,
        ],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TN) low CSP -- not in binding site ({len(tn_distances)})",
            f"(FP) CSP -- not in binding site ({len(fp_distances)})",
            f"(FN) low CSP --  in binding site ({len(fn_distances)})",
            f"(TP) CSP --  in binding site ({len(tp_distances)})",
        ],
    )
    ax = plt.gca()
    plt.xlabel(FIGURE_3_DISTANCE_XLABEL)
    plt.ylabel("Number of Residues")
    if title.strip():
        plt.title(title)
    ax.tick_params(axis="both", labelsize=14)
    plt.legend(frameon=True, fontsize=14, title_fontsize=14)
    plt.xlim(0, x_max)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_figure_3_combined(
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    fn_distances: Sequence[float],
    tn_distances: Sequence[float],
    output_path: Path,
    *,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    """Stacked Figure 3 (panel a above panel b) with labels and FP/(TP+FP) on a."""
    _set_plot_style(dpi)
    fig, (ax_a, ax_b) = plt.subplots(
        2,
        1,
        figsize=(fig_width, fig_height),
        constrained_layout=True,
    )

    # --- Panel a: TP / FP ---
    all_a = list(tp_distances) + list(fp_distances)
    bins_a, x_max_a = _get_bins(all_a, bin_width, max_distance)
    ax_a.hist(
        [tp_distances, fp_distances],
        bins=bins_a,
        stacked=True,
        color=[classification_colors.TP, classification_colors.FP],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TP) CSP --  in binding site ({len(tp_distances)})",
            f"(FP) CSP --  not in binding site ({len(fp_distances)})",
        ],
    )
    ax_a.set_xlabel(FIGURE_3_DISTANCE_XLABEL)
    ax_a.set_ylabel("Number of Residues")
    ax_a.tick_params(axis="both", labelsize=14)
    ax_a.legend(frameon=True, fontsize=12)
    ax_a.set_xlim(0, x_max_a)

    n_tp, n_fp = len(tp_distances), len(fp_distances)
    denom = n_tp + n_fp
    fp_of_sig = (100.0 * n_fp / denom) if denom else float("nan")
    ax_a.text(
        0.50,
        0.55,
        f"{fp_of_sig:.1f}% FP",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=20,
        fontweight="bold",
        color="red",
        zorder=5,
    )
    ax_a.text(
        -0.08,
        1.02,
        "A.",
        transform=ax_a.transAxes,
        ha="left",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        color="black",
        clip_on=False,
    )

    # --- Panel b: full confusion matrix ---
    all_b = all_a + list(fn_distances) + list(tn_distances)
    bins_b, x_max_b = _get_bins(all_b, bin_width, max_distance)
    ax_b.hist(
        [tn_distances, fp_distances, fn_distances, tp_distances],
        bins=bins_b,
        stacked=True,
        color=[
            classification_colors.TN,
            classification_colors.FP,
            classification_colors.FN,
            classification_colors.TP,
        ],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TN) low CSP -- not in binding site ({len(tn_distances)})",
            f"(FP) CSP -- not in binding site ({len(fp_distances)})",
            f"(FN) low CSP --  in binding site ({len(fn_distances)})",
            f"(TP) CSP --  in binding site ({len(tp_distances)})",
        ],
    )
    ax_b.set_xlabel(FIGURE_3_DISTANCE_XLABEL)
    ax_b.set_ylabel("Number of Residues")
    ax_b.tick_params(axis="both", labelsize=14)
    ax_b.legend(frameon=True, fontsize=11)
    ax_b.set_xlim(0, x_max_b)
    ax_b.text(
        -0.08,
        1.02,
        "B.",
        transform=ax_b.transAxes,
        ha="left",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        color="black",
        clip_on=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Figure 3A and Figure 3B histograms.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--output-a", type=Path, default=None, help="Output path for figure_3_a.png")
    parser.add_argument("--output-b", type=Path, default=None, help="Output path for figure_3_b.png")
    parser.add_argument(
        "--output-combined",
        type=Path,
        default=None,
        help="Path for stacked combined figure (a above b). Default: figures/figure_3.png.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ_ph0.5_temp5C.csv",
        help="CSV with 'holo_pdb' column (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument("--targets", type=str, help="Optional comma-separated holo_pdb list.")
    parser.add_argument(
        "--title-a",
        type=str,
        default="",
    )
    parser.add_argument(
        "--title-b",
        type=str,
        default="",
    )
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--max-distance", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--fig-width", type=float, default=8.0)
    parser.add_argument("--fig-height", type=float, default=5.0)
    parser.add_argument(
        "--combined-fig-width",
        type=float,
        default=8.0,
        help="Figure width for --output-combined (default: %(default)s).",
    )
    parser.add_argument(
        "--combined-fig-height",
        type=float,
        default=10.5,
        help="Figure height for --output-combined (default: %(default)s).",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    targets_csv = args.targets_csv
    if targets_csv is not None and not targets_csv.is_absolute():
        targets_csv = _REPO_ROOT / targets_csv
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    figures_dir = args.figures_dir if args.figures_dir.is_absolute() else _REPO_ROOT / args.figures_dir
    selected_dirs = _resolve_selected_dirs(outputs_dir, targets_csv, args.targets)
    tp, fp, fn, tn = collect_distance_categories(outputs_dir, selected_dirs)

    output_a = args.output_a or (figures_dir / "figure_3_a.png")
    output_b = args.output_b or (figures_dir / "figure_3_b.png")
    output_c = args.output_combined or (figures_dir / "figure_3.png")
    if not output_a.is_absolute():
        output_a = _REPO_ROOT / output_a
    if not output_b.is_absolute():
        output_b = _REPO_ROOT / output_b
    if not output_c.is_absolute():
        output_c = _REPO_ROOT / output_c

    plot_figure_3a(
        tp,
        fp,
        output_a,
        title=args.title_a,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )
    plot_figure_3b(
        tp,
        fp,
        fn,
        tn,
        output_b,
        title=args.title_b,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )

    print(f"Figure 3A written to {output_a.resolve()}")
    print(f"Figure 3B written to {output_b.resolve()}")

    plot_figure_3_combined(
        tp,
        fp,
        fn,
        tn,
        output_c,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
        dpi=args.dpi,
        fig_width=args.combined_fig_width,
        fig_height=args.combined_fig_height,
    )
    print(f"Combined Figure 3 written to {output_c.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
