#!/usr/bin/env python3
"""
SI Fig. S21 — F1 vs MCC scatterplot.

Reuses existing logic from scripts/plot_f1_vs_mcc.py and writes:
  ./figures/SF21_f1_vs_mcc.png

By default only systems whose system_id matches holo_pdb entries in
CSP_UBQ_ph0.5_temp5C.csv are plotted. Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .plot_f1_vs_mcc import load_f1_mcc, load_holo_pdb_ids_from_targets_csv, plot_f1_vs_mcc
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.plot_f1_vs_mcc import (  # type: ignore
        load_f1_mcc,
        load_holo_pdb_ids_from_targets_csv,
        plot_f1_vs_mcc,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S21 (F1 vs MCC scatter plot)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with holo_pdb column to filter system_id rows (default: CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF21_f1_vs_mcc.png",
        help="Destination for SI Fig. S21 image.",
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

    if not input_csv.exists():
        print(f"Error: input CSV does not exist: {input_csv}", file=sys.stderr)
        return 1
    if not targets_csv.exists():
        print(f"Error: targets CSV does not exist: {targets_csv}", file=sys.stderr)
        return 1

    holo_filter = load_holo_pdb_ids_from_targets_csv(targets_csv)
    f1_vals, mcc_vals = load_f1_mcc(input_csv, holo_pdb_filter=holo_filter)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    plot_f1_vs_mcc(f1_vals, mcc_vals, output_image)
    print(f"SI Fig. S21 saved to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
