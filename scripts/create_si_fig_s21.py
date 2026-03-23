#!/usr/bin/env python3
"""
SI Fig. S20 — CSP significance threshold histogram.

For each target, recomputes the significance threshold from csp_A values using
the same logic as the pipeline (compute_threshold_with_outlier_removal with
default config). Plots a histogram of these per-target thresholds.

Output: figures/SF20_significance_threshold.png

By default only targets listed in CSP_UBQ_ph0.5_temp5C.csv (holo_pdb) are included.
Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path
from typing import List, Set

import matplotlib.pyplot as plt
import pandas as pd

try:
    from .config import thresholds
    from .csp import compute_threshold_with_outlier_removal
    from .plot_f1_vs_mcc import load_holo_pdb_ids_from_targets_csv
except Exception:
    import os as _os
    import sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import thresholds
    from scripts.csp import compute_threshold_with_outlier_removal
    from scripts.plot_f1_vs_mcc import load_holo_pdb_ids_from_targets_csv


def _output_dir_matches_targets(dir_name: str, holo_lower: Set[str]) -> bool:
    t = dir_name.lower()
    return t in holo_lower or t.split("_")[0] in holo_lower


def collect_thresholds(outputs_dir: Path, holo_lower: Set[str]) -> List[tuple[str, float]]:
    """Collect (target_name, threshold) for each target with valid csp_A data."""
    results: List[tuple[str, float]] = []

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if not _output_dir_matches_targets(target_name, holo_lower):
            continue
        df = pd.read_csv(alignment_path)

        if "csp_A" not in df.columns:
            continue

        csp_vals = pd.to_numeric(df["csp_A"], errors="coerce")
        valid = csp_vals.notna()
        values = csp_vals[valid].astype(float).tolist()

        if len(values) < 2:
            continue

        info = compute_threshold_with_outlier_removal(
            values,
            outlier_z=thresholds.outlier_z_score,
            significance_z=thresholds.significance_z_score,
            max_iterations=thresholds.max_outlier_iterations,
            max_outlier_fraction=thresholds.max_outlier_fraction,
        )
        results.append((target_name, info.threshold))

    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S20 (CSP significance threshold histogram)"
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with holo_pdb column (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument("--output", type=Path, default=Path("figures") / "SF20_significance_threshold.png")
    parser.add_argument("--bin-width", type=float, default=0.02)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--fig-width", type=float, default=8.0)
    parser.add_argument("--fig-height", type=float, default=5.0)
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    targets_csv = args.targets_csv
    if not targets_csv.is_absolute():
        targets_csv = project_root / targets_csv
    output_path = args.output if args.output.is_absolute() else project_root / args.output

    if not outputs_dir.exists():
        print(f"Error: outputs directory does not exist: {outputs_dir}", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    holo_set = load_holo_pdb_ids_from_targets_csv(targets_csv)
    data = collect_thresholds(outputs_dir, holo_set)
    if not data:
        print("No threshold data found.", file=sys.stderr)
        return 1

    sorted_data = sorted(data, key=lambda x: x[1], reverse=True)
    _, thresh_vals = zip(*sorted_data)
    thresholds_list = list(thresh_vals)

    print(f"Collected {len(thresholds_list)} thresholds")
    print("\nHolo PDBs by decreasing significance threshold:")
    for target, thresh in sorted_data:
        print(f"  {target}: {thresh:.4f}")
    mean_val = statistics.mean(thresholds_list)
    median_val = statistics.median(thresholds_list)
    print(f"  Min: {min(thresholds_list):.4f}")
    print(f"  Max: {max(thresholds_list):.4f}")
    print(f"  Mean: {mean_val:.4f}")
    print(f"  Median: {median_val:.4f}")

    bin_width = args.bin_width
    t_min, t_max = min(thresholds_list), max(thresholds_list)
    n_bins = max(1, int((t_max - t_min) / bin_width) + 1)
    bins = [t_min + i * bin_width for i in range(n_bins + 1)]

    plt.figure(figsize=(args.fig_width, args.fig_height))
    plt.hist(
        thresholds_list,
        bins=bins,
        edgecolor="black",
        linewidth=0.8,
    )
    plt.axvline(mean_val, color="red", linestyle="--", linewidth=2, label=f"Mean: {mean_val:.4f}")
    plt.axvline(median_val, color="blue", linestyle="--", linewidth=2, label=f"Median: {median_val:.4f}")
    plt.xlabel("CSP Significance Threshold (ppm)")
    plt.ylabel("Number of Targets")
    plt.title("Distribution of CSP Significance Thresholds Across Dataset")
    plt.legend()
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=args.dpi)
    plt.close()

    print(f"SI Fig. S20 written to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
