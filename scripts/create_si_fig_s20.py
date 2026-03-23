#!/usr/bin/env python3
"""
SI Fig. S19 — ideal N/H offsets (offset grid heatmap from grid search).

Reads grid search CSV files from outputs/, creates the heatmap, and saves to
./figures/SF19_ideal_offsets.png.

By default only targets listed in CSP_UBQ_ph0.5_temp5C.csv (holo_pdb) are included.
Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Set

try:
    from .analyze_offsets import collect_best_grid_offsets, create_grid_heatmap
    from .plot_f1_vs_mcc import load_holo_pdb_ids_from_targets_csv
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_offsets import collect_best_grid_offsets, create_grid_heatmap  # type: ignore
    from scripts.plot_f1_vs_mcc import load_holo_pdb_ids_from_targets_csv  # type: ignore


def _output_dir_matches_targets(dir_name: str, holo_lower: Set[str]) -> bool:
    t = dir_name.lower()
    return t in holo_lower or t.split("_")[0] in holo_lower


def main() -> int:
    parser = argparse.ArgumentParser(description="Create SI Fig. S19 (ideal N/H offsets heatmap).")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"), help="Path to outputs directory")
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"), help="Path to figures directory")
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with holo_pdb column (default: CSP_UBQ_ph0.5_temp5C.csv).",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    outputs_dir = project_root / args.outputs_dir if not args.outputs_dir.is_absolute() else args.outputs_dir
    figures_dir = project_root / args.figures_dir if not args.figures_dir.is_absolute() else args.figures_dir

    targets_csv = args.targets_csv
    if not targets_csv.is_absolute():
        targets_csv = project_root / targets_csv

    if not outputs_dir.exists():
        print(f"Error: outputs directory '{outputs_dir}' does not exist")
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    holo_set = load_holo_pdb_ids_from_targets_csv(targets_csv)

    (
        _all_h,
        _all_n,
        h_by_target,
        n_by_target,
        h_grid_values,
        n_grid_values,
    ) = collect_best_grid_offsets(str(outputs_dir))

    keys = sorted(
        k
        for k in h_by_target
        if k in n_by_target and _output_dir_matches_targets(k, holo_set)
    )
    best_h_offsets = [h_by_target[k] for k in keys]
    best_n_offsets = [n_by_target[k] for k in keys]

    if not best_h_offsets or not best_n_offsets:
        print(
            "No grid offset data found for selected targets. Run the pipeline and/or check --targets-csv.",
            file=sys.stderr,
        )
        return 1

    figures_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = str(figures_dir)
    heatmap_path = create_grid_heatmap(
        best_h_offsets,
        best_n_offsets,
        temp_dir,
        h_grid_values=h_grid_values,
        n_grid_values=n_grid_values,
    )

    if not heatmap_path:
        print("Failed to create heatmap.")
        return 1

    suppl_path = figures_dir / "SF19_ideal_offsets.png"
    os.rename(heatmap_path, suppl_path)
    print(f"SI Fig. S19 saved to {suppl_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
