#!/usr/bin/env python3
"""
SI Fig. S17 — F1 vs MCC scatterplot.

Reuses existing logic from scripts/plot_f1_vs_mcc.py and writes:
  ./figures/SF17_f1_vs_mcc.png

By default only systems whose ``system_id`` matches a resolved ``outputs/<dir>``
basename from ``--targets-csv`` rows (``apo_bmrb``/``holo_bmrb`` congruence with
``master_alignment.csv``; see ``scripts.target_resolution``) are plotted.
Override with --targets-csv / --outputs-dir.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .plot_f1_vs_mcc import load_f1_mcc, plot_f1_vs_mcc
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.plot_f1_vs_mcc import load_f1_mcc, plot_f1_vs_mcc  # type: ignore
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S17 (F1 vs MCC scatter plot)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv.",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root outputs directory used to resolve targets-csv rows to system_id (default: outputs).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help=(
            "CSV with holo_pdb plus apo_bmrb/holo_bmrb to filter system_id rows "
            "(default: data/CSP_UBQ_ph0.5_temp5C.csv)."
        ),
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF17_f1_vs_mcc.png",
        help="Destination for SI Fig. S17 image.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    input_csv = args.input if args.input.is_absolute() else project_root / args.input
    output_image = args.output_image if args.output_image.is_absolute() else project_root / args.output_image
    targets_csv = args.targets_csv
    if not targets_csv.is_absolute():
        targets_csv = project_root / targets_csv
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir

    if not input_csv.exists():
        print(f"Error: input CSV does not exist: {input_csv}", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1
    if not outputs_dir.is_dir():
        print(f"Error: outputs directory does not exist: {outputs_dir}", file=sys.stderr)
        return 1

    rows = load_target_rows(targets_csv)
    allowed_system_ids = {p.name for p in resolve_target_rows(rows, outputs_dir)}
    f1_vals, mcc_vals = load_f1_mcc(input_csv, allowed_system_ids=allowed_system_ids)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    plot_f1_vs_mcc(f1_vals, mcc_vals, output_image)
    print(f"SI Fig. S17 saved to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
