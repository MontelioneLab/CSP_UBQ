#!/usr/bin/env python3
"""
SI Fig. S17 — ideal N/H and HA/CA offsets (offset grid heatmaps from grid search).

Reads grid search CSV files from outputs/, creates the heatmaps, and saves to
./figures/SF17_ideal_offsets.png (N/H) and ./figures/SF17_ideal_offsets_ha_ca.png
(HA/CA). HA/CA uses grid CSVs only (no csp_table_HA_CA.csv fallback).

By default only targets listed in CSP_UBQ_ph0.5_temp5C.csv are included.
Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
try:
    from .analyze_offsets import (
        collect_best_grid_offsets,
        collect_best_grid_offsets_ha_ca,
        create_grid_heatmap,
    )
    from .config import Referencing
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_offsets import (  # type: ignore
        collect_best_grid_offsets,
        collect_best_grid_offsets_ha_ca,
        create_grid_heatmap,
    )
    from scripts.config import Referencing  # type: ignore
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S17 (ideal N/H and HA/CA offset heatmaps)."
    )
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
            "No N/H grid offset data found for selected targets. Run the pipeline and/or check --targets-csv.",
            file=sys.stderr,
        )
        return 1

    print(
        f"[SF17] Aggregating ideal N/H offsets for n={len(keys)} targets "
        f"from {targets_csv.name} (resolved {len(selected_dir_names)} dirs)"
    )

    figures_dir.mkdir(parents=True, exist_ok=True)
    suppl_path = figures_dir / "SF17_ideal_offsets.png"
    heatmap_path = create_grid_heatmap(
        best_h_offsets,
        best_n_offsets,
        str(figures_dir),
        h_grid_values=h_grid_values,
        n_grid_values=n_grid_values,
        output_name=suppl_path.name,
    )
    if not heatmap_path:
        print("Failed to create N/H heatmap.", file=sys.stderr)
        return 1
    print(f"SI Fig. S17 N/H panel saved to {suppl_path.resolve()} (n={len(keys)})")

    (
        _all_ha,
        _all_ca,
        ha_by_target,
        ca_by_target,
        ha_grid_values,
        ca_grid_values,
    ) = collect_best_grid_offsets_ha_ca(str(outputs_dir))

    ha_ca_keys = sorted(
        k for k in ha_by_target if k in ca_by_target and k in selected_dir_names
    )
    best_ha_offsets = [ha_by_target[k] for k in ha_ca_keys]
    best_ca_offsets = [ca_by_target[k] for k in ha_ca_keys]

    if not best_ha_offsets or not best_ca_offsets:
        print(
            f"No HA/CA grid offset CSVs found for selected targets "
            f"(resolved {len(selected_dir_names)} dirs; "
            f"{len(ha_by_target)} HA grids, {len(ca_by_target)} CA grids overall). "
            "HA/CA panel requires cached offset_grid_HA_*__CA_*.csv files.",
            file=sys.stderr,
        )
        return 1

    print(
        f"[SF17] Aggregating ideal HA/CA offsets for n={len(ha_ca_keys)} targets "
        f"from {targets_csv.name} (resolved {len(selected_dir_names)} dirs)"
    )

    ref = Referencing()
    ha_ca_path = figures_dir / "SF17_ideal_offsets_ha_ca.png"
    ha_ca_heatmap = create_grid_heatmap(
        best_ha_offsets,
        best_ca_offsets,
        str(figures_dir),
        h_grid_values=ha_grid_values,
        n_grid_values=ca_grid_values,
        y_label="HA offset (ppm)",
        x_label="CA offset (ppm)",
        title="Best HA/CA Offsets from Grid Search Across Targets",
        y_min=ref.grid_ha_min,
        y_max=ref.grid_ha_max,
        y_step=ref.grid_ha_step,
        x_min=ref.grid_ca_min,
        x_max=ref.grid_ca_max,
        x_step=ref.grid_ca_step,
        output_name=ha_ca_path.name,
    )
    if not ha_ca_heatmap:
        print("Failed to create HA/CA heatmap.", file=sys.stderr)
        return 1
    print(f"SI Fig. S17 HA/CA panel saved to {ha_ca_path.resolve()} (n={len(ha_ca_keys)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
