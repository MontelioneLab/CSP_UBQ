#!/usr/bin/env python3
"""
SI Fig. S15 — ideal N/H offsets (offset grid heatmap from grid search).

Reads grid search CSV files from outputs/, creates the heatmap, and saves to
./figures/SF15_ideal_offsets.png.

By default only targets listed in CSP_UBQ_ph0.5_temp5C.csv (holo_pdb) are included.
Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
try:
    from .analyze_offsets import collect_best_grid_offsets, create_grid_heatmap
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_offsets import collect_best_grid_offsets, create_grid_heatmap  # type: ignore
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def main() -> int:
    parser = argparse.ArgumentParser(description="Create SI Fig. S15 (ideal N/H offsets heatmap).")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"), help="Path to outputs directory")
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"), help="Path to figures directory")
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with holo_pdb column (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
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

    rows = load_target_rows(targets_csv)
    selected_dir_names = {p.name for p in resolve_target_rows(rows, outputs_dir)}

    (
        _all_h,
        _all_n,
        h_by_target,
        n_by_target,
        h_grid_values,
        n_grid_values,
    ) = collect_best_grid_offsets(str(outputs_dir))

    keys = sorted(
        k for k in h_by_target if k in n_by_target and k in selected_dir_names
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

    suppl_path = figures_dir / "SF15_ideal_offsets.png"
    os.rename(heatmap_path, suppl_path)
    print(f"SI Fig. S15 saved to {suppl_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
