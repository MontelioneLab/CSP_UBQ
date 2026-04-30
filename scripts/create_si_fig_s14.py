#!/usr/bin/env python3
"""
Generator script filename (s14) ≠ SI index: this produces SI Fig. S11.

SI Fig. S11 — CA-inclusive vs exclusive CSP F1 scores scatterplot.

Reuses existing logic from scripts/analyze_targets_ca.py and writes:
  ./figures/SF11_f1_ca_vs_exclusive.png

Default targets list: CSP_UBQ_ph0.5_temp5C.csv (buffer-filtered subset). Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .analyze_targets_ca import (
        collect_nh_results,
        collect_results,
        render_f1_comparison_scatterplot,
    )
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_targets_ca import (  # type: ignore
        collect_nh_results,
        collect_results,
        render_f1_comparison_scatterplot,
    )
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S11 (CA-inclusive vs exclusive CSP F1 scores)."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root outputs directory with per-target subdirectories.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV file containing holo_pdb targets (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF11_f1_ca_vs_exclusive.png",
        help="Destination for SI Fig. S11.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else project_root / args.targets_csv
    output_image = args.output_image if args.output_image.is_absolute() else project_root / args.output_image

    if not outputs_dir.exists():
        print(f"Error: outputs directory does not exist: {outputs_dir}", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    rows = load_target_rows(targets_csv)
    allowed_targets = {p.name for p in resolve_target_rows(rows, outputs_dir)}
    ca_results, _, _ = collect_results(outputs_dir, allowed_targets)
    nh_results = collect_nh_results(outputs_dir, allowed_targets)

    if not ca_results:
        print("No CA-inclusive results found for selected targets.", file=sys.stderr)
        return 1
    if not nh_results:
        print("No N/H results found for selected targets.", file=sys.stderr)
        return 1

    output_image.parent.mkdir(parents=True, exist_ok=True)
    render_f1_comparison_scatterplot(ca_results, nh_results, output_image)
    print(f"SI Fig. S11 saved to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
