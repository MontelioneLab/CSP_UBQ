#!/usr/bin/env python3
"""
Create publication-ready Figure 1 panels from `outputs/*/master_alignment.csv`:

- figure_1_a.png: stacked histogram of significant residues (TP vs FP only)
- figure_1_b.png: stacked histogram of full confusion matrix (TN, FP, FN, TP)

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
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors


SIGNIFICANT_COLUMN = "significant"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _parse_targets(targets_csv: Optional[Path], targets_str: Optional[str]) -> Optional[Set[str]]:
    selected: Set[str] = set()
    if targets_csv:
        df = pd.read_csv(targets_csv)
        if "holo_pdb" not in df.columns:
            raise ValueError(f"Missing 'holo_pdb' column in {targets_csv}")
        selected.update(df["holo_pdb"].astype(str).str.strip().str.lower().tolist())
    if targets_str:
        selected.update(t.strip().lower() for t in targets_str.split(",") if t.strip())
    return selected or None


def _target_is_selected(target_name: str, selected: Optional[Set[str]]) -> bool:
    if selected is None:
        return True
    target_lower = target_name.lower()
    # Support duplicated-output folders like 2mur_1, 2mur_2.
    base_target = target_lower.split("_")[0]
    return target_lower in selected or base_target in selected


def _get_bins(data: Sequence[float], bin_width: float, max_distance: Optional[float]) -> tuple[List[float], float]:
    if not data:
        raise ValueError("No CA-distance data found for selected targets.")
    data_max = max(data) if max_distance is None else max_distance
    n_bins = max(1, int((data_max + bin_width) / bin_width))
    bins = [i * bin_width for i in range(n_bins + 1)]
    return bins, data_max


def collect_distance_categories(
    outputs_dir: Path,
    selected_targets: Optional[Set[str]] = None,
) -> tuple[List[float], List[float], List[float], List[float]]:
    tp_distances: List[float] = []
    fp_distances: List[float] = []
    fn_distances: List[float] = []
    tn_distances: List[float] = []

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if not _target_is_selected(target_name, selected_targets):
            continue

        df = pd.read_csv(alignment_path)
        required_cols = {SIGNIFICANT_COLUMN, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"[WARN] Skipping {alignment_path}: missing columns {missing}")
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


def plot_figure_1a(
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
    plt.xlabel("Minimum CA Distance (Å)")
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


def plot_figure_1b(
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
            f"(TN) No CSP -- not in binding site ({len(tn_distances)})",
            f"(FP) CSP -- not in binding site ({len(fp_distances)})",
            f"(FN) No CSP --  in binding site ({len(fn_distances)})",
            f"(TP) CSP --  in binding site ({len(tp_distances)})",
        ],
    )
    ax = plt.gca()
    plt.xlabel("Minimum CA Distance (Å)")
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


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Figure 1A and Figure 1B histograms.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--output-a", type=Path, default=None, help="Output path for figure_1_a.png")
    parser.add_argument("--output-b", type=Path, default=None, help="Output path for figure_1_b.png")
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO_ROOT / "CSP_UBQ_ph0.5_temp5C.csv",
        help="CSV with 'holo_pdb' column (default: CSP_UBQ_ph0.5_temp5C.csv).",
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
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    targets_csv = args.targets_csv
    if not targets_csv.is_absolute():
        targets_csv = _REPO_ROOT / targets_csv
    selected_targets = _parse_targets(targets_csv, args.targets)
    tp, fp, fn, tn = collect_distance_categories(args.outputs_dir, selected_targets)

    output_a = args.output_a or (args.figures_dir / "figure_1_a.png")
    output_b = args.output_b or (args.figures_dir / "figure_1_b.png")

    plot_figure_1a(
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
    plot_figure_1b(
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

    print(f"Figure 1A written to {output_a.resolve()}")
    print(f"Figure 1B written to {output_b.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
