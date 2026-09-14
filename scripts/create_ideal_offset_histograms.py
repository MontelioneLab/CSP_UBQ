#!/usr/bin/env python3
"""
1D histograms of per-target ideal referencing offsets for three CSP families.

Writes analysis-only PNGs under outputs/ideal_offsets/ (not SI):
  - ideal_offsets_hn.png      (H, N) from 2D HN grid CSVs
  - ideal_offsets_hn_ca.png   (H, N, CA) from csp_table_CA.csv (3D grid)
  - ideal_offsets_ha_ca.png   (HA, CA) from 2D HA/CA grid CSVs

Default target set: CSP_UBQ_ph0.5_temp5C.csv resolved to outputs dirs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt

try:
    from .analyze_offsets import (
        _calculate_bin_edges,
        collect_best_grid_offsets,
        collect_best_grid_offsets_ha_ca,
        collect_hn_ca_offsets_from_ca_tables,
    )
    from .config import Referencing
    from .csp import _frange
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_offsets import (  # type: ignore
        _calculate_bin_edges,
        collect_best_grid_offsets,
        collect_best_grid_offsets_ha_ca,
        collect_hn_ca_offsets_from_ca_tables,
    )
    from scripts.config import Referencing  # type: ignore
    from scripts.csp import _frange  # type: ignore
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def _grid_centers(min_value: float, max_value: float, step: float) -> List[float]:
    decimals = max(0, min(6, len(str(step).split(".")[-1]) if "." in str(step) else 0))
    return _frange(float(min_value), float(max_value), float(step), decimals)


def _write_histogram_figure(
    panels: Sequence[Tuple[str, List[float], List[float]]],
    *,
    title: str,
    output_path: Path,
) -> None:
    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.4 * n_panels, 4.2), squeeze=False)
    cmap = plt.get_cmap("viridis")
    colors = [cmap(0.35), cmap(0.55), cmap(0.75)]

    for i, (label, values, grid_centers) in enumerate(panels):
        ax = axes[0, i]
        edges = _calculate_bin_edges(grid_centers)
        ax.hist(
            values,
            bins=edges,
            color=colors[i % len(colors)],
            edgecolor="black",
            alpha=0.75,
        )
        ax.set_xlabel(f"{label} offset (ppm)", fontsize=12)
        if i == 0:
            ax.set_ylabel("Count", fontsize=12)
        ax.set_title(label, fontsize=13)
        ax.set_xlim(float(edges[0]), float(edges[-1]))
        ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write 1D ideal-offset histograms for H/N, H/N/CA, and HA/CA CSPs."
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs") / "ideal_offsets",
        help="Directory for histogram PNGs (default: outputs/ideal_offsets).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    out_dir = args.out_dir if args.out_dir.is_absolute() else project_root / args.out_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else project_root / args.targets_csv

    if not outputs_dir.exists():
        print(f"Error: outputs directory '{outputs_dir}' does not exist", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    rows = load_target_rows(targets_csv)
    selected_dir_names = {p.name for p in resolve_target_rows(rows, outputs_dir)}
    ref = Referencing()
    written = 0

    (
        _all_h,
        _all_n,
        h_by_target,
        n_by_target,
        _h_grid,
        _n_grid,
    ) = collect_best_grid_offsets(str(outputs_dir))
    hn_keys = sorted(
        k for k in h_by_target if k in n_by_target and k in selected_dir_names
    )
    if hn_keys:
        h_vals = [h_by_target[k] for k in hn_keys]
        n_vals = [n_by_target[k] for k in hn_keys]
        hn_path = out_dir / "ideal_offsets_hn.png"
        _write_histogram_figure(
            [
                ("H", h_vals, _grid_centers(ref.grid_h_min, ref.grid_h_max, ref.grid_h_step)),
                ("N", n_vals, _grid_centers(ref.grid_n_min, ref.grid_n_max, ref.grid_n_step)),
            ],
            title=f"Ideal H/N offsets (n={len(hn_keys)})",
            output_path=hn_path,
        )
        print(f"[ideal-offsets] H/N n={len(hn_keys)} -> {hn_path}")
        written += 1
    else:
        print(
            f"[ideal-offsets] H/N: no matching grid offsets "
            f"(resolved {len(selected_dir_names)} dirs); skipping PNG."
        )

    h_ca, n_ca, ca_ca = collect_hn_ca_offsets_from_ca_tables(str(outputs_dir))
    hnca_keys = sorted(
        k for k in h_ca if k in n_ca and k in ca_ca and k in selected_dir_names
    )
    if hnca_keys:
        hnca_path = out_dir / "ideal_offsets_hn_ca.png"
        _write_histogram_figure(
            [
                ("H", [h_ca[k] for k in hnca_keys], _grid_centers(ref.grid_h_min, ref.grid_h_max, ref.grid_h_step)),
                ("N", [n_ca[k] for k in hnca_keys], _grid_centers(ref.grid_n_min, ref.grid_n_max, ref.grid_n_step)),
                ("CA", [ca_ca[k] for k in hnca_keys], _grid_centers(ref.grid_ca_min, ref.grid_ca_max, ref.grid_ca_step)),
            ],
            title=f"Ideal H/N/CA offsets (n={len(hnca_keys)})",
            output_path=hnca_path,
        )
        print(f"[ideal-offsets] H/N/CA n={len(hnca_keys)} -> {hnca_path}")
        written += 1
    else:
        print(
            f"[ideal-offsets] H/N/CA: no complete triples in csp_table_CA.csv "
            f"(resolved {len(selected_dir_names)} dirs); skipping PNG."
        )

    (
        _all_ha,
        _all_ca,
        ha_by_target,
        ca_by_target,
        _ha_grid,
        _ca_grid,
    ) = collect_best_grid_offsets_ha_ca(str(outputs_dir))
    haca_keys = sorted(
        k for k in ha_by_target if k in ca_by_target and k in selected_dir_names
    )
    if haca_keys:
        haca_path = out_dir / "ideal_offsets_ha_ca.png"
        _write_histogram_figure(
            [
                ("HA", [ha_by_target[k] for k in haca_keys], _grid_centers(ref.grid_ha_min, ref.grid_ha_max, ref.grid_ha_step)),
                ("CA", [ca_by_target[k] for k in haca_keys], _grid_centers(ref.grid_ca_min, ref.grid_ca_max, ref.grid_ca_step)),
            ],
            title=f"Ideal HA/CA offsets (n={len(haca_keys)})",
            output_path=haca_path,
        )
        print(f"[ideal-offsets] HA/CA n={len(haca_keys)} -> {haca_path}")
        written += 1
    else:
        print(
            f"[ideal-offsets] HA/CA: no matching grid offsets "
            f"(resolved {len(selected_dir_names)} dirs); skipping PNG."
        )

    if written == 0:
        print("No ideal-offset histograms written.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
