#!/usr/bin/env python3
"""
SI Fig. S12 — CA-inclusive vs exclusive CSP F1 scores scatterplot.

Generator script ``create_si_fig_s12.py`` matches output prefix ``SF12_``.

Reuses existing logic from scripts/analyze_targets_ca.py and writes:
  ./figures/SF12_f1_ca_vs_exclusive.png

Targets are trimmed with the **same CA-shift eligibility rule** as SI Fig. S15 /
the standalone 1D F1 boxplot: among resolved pipeline directories,
``target_basenames_passing_ca_shift_coverage`` retains only outputs whose
``1d_analysis.csv`` has both ``CA_apo`` and ``CA_holo`` on strictly more than
``--min-ca-coverage`` of rows (default
``DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE`` = 0.5). Same gate as SI Fig. S11 / S15 /
the standalone 1D F1 boxplot. CA-inclusive vs N/H F1 scores are collected on that subset.

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
    from .analyze_targets_single_atom_shifts import (
        DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        target_basenames_passing_ca_shift_coverage,
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
    from scripts.analyze_targets_single_atom_shifts import (  # type: ignore
        DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        target_basenames_passing_ca_shift_coverage,
    )
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Fig. S12 (CA-inclusive vs exclusive CSP F1 scores)."
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
        default=Path("figures") / "SF12_f1_ca_vs_exclusive.png",
        help="Destination for SI Fig. S12.",
    )
    parser.add_argument(
        "--min-ca-coverage",
        type=float,
        default=DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        help=(
            "Same as SI Fig. S15: require strictly more than this fraction of "
            "1d_analysis.csv rows with both CA_apo and CA_holo (default: %(default)s)."
        ),
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
    if not allowed_targets:
        print("No targets resolved from CSV against outputs/", file=sys.stderr)
        return 1

    min_cov = float(args.min_ca_coverage)
    eligible, coverage_map = target_basenames_passing_ca_shift_coverage(
        outputs_dir,
        min_coverage=min_cov,
        allowed_basenames={k: True for k in allowed_targets},
    )
    print(
        f"{len(eligible)} targets pass CA row coverage > {min_cov:.0%} "
        f"(among {len(allowed_targets)} CSV-resolved; "
        f"{len(coverage_map)} with readable 1d_analysis CA columns)"
    )
    if not eligible:
        print(
            "No targets left after CA shift coverage filter; cannot render SI Fig. S12.",
            file=sys.stderr,
        )
        return 1

    ca_results, _, _ = collect_results(outputs_dir, eligible, "nh_ca")
    nh_results = collect_nh_results(outputs_dir, eligible)

    if not ca_results:
        print("No CA-inclusive results found for selected targets.", file=sys.stderr)
        return 1
    if not nh_results:
        print("No N/H results found for selected targets.", file=sys.stderr)
        return 1

    output_image.parent.mkdir(parents=True, exist_ok=True)
    render_f1_comparison_scatterplot(ca_results, nh_results, output_image, "nh_ca")
    print(f"SI Fig. S12 saved to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
