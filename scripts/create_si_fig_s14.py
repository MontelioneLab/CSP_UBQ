#!/usr/bin/env python3
"""
SI Fig. S14 — Per-target **F1 scores** for 1D H / N / Cα CSPs (boxplots + paired Wilcoxon).

Implements publication SI Fig. S14 by delegating to
:func:`run_f1_1d_boxplot` in ``create_si_fig_f1_1d_boxplot`` (same F1 /
significance logic as ``SF_f1_1d_boxplot.png``). Targets are filtered with
:class:`scripts.target_resolution` so canonical ``outputs/{HOLO}_{apo_bmrb}/``
directories match ``CSP_UBQ``-style rows.

Older versions of this script incorrectly plotted summarized **|1D CSP|**
magnitudes on the *y*-axis rather than classifier **F1** scores derived from the
same 1D significance rules as the rest of the 1D single-atom analysis.

Outputs:

  ./figures/SF14_1d_CSP_boxplot.png
  ./figures/SF14_1d_CSP_boxplot_stats.csv  (Holm-adjusted pairwise *p*-values)

Default targets: ``data/CSP_UBQ_ph0.5_temp5C.csv``. Override via ``--targets-csv``.
CA-shift gating uses :func:`target_basenames_passing_ca_shift_coverage` in
``analyze_targets_single_atom_shifts`` (same as SI Fig. S10 / S11; ``--min-ca-coverage``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .analyze_targets_single_atom_shifts import DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE
    from .create_si_fig_f1_1d_boxplot import run_f1_1d_boxplot
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.analyze_targets_single_atom_shifts import (  # type: ignore
        DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
    )
    from scripts.create_si_fig_f1_1d_boxplot import run_f1_1d_boxplot  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="SI Fig. S14: 1D H/N/Cα F1 score boxplots with paired Wilcoxon (Holm-adjusted)."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root outputs directory with per-target folders.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="Targets CSV (holo_pdb + apo_bmrb); resolved to output dirs (default: buffer-filtered CSP_UBQ).",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF14_1d_CSP_boxplot.png",
        help="Destination PNG for SI Fig. S14.",
    )
    parser.add_argument(
        "--stats-csv",
        type=Path,
        default=None,
        help="Optional Holm stats table path (default: sibling <stem>_stats.csv next to --output-image).",
    )
    parser.add_argument(
        "--min-ca-coverage",
        type=float,
        default=DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        help=(
            "Same as SI Fig. S11: strictly more than this fraction of 1d_analysis.csv "
            "rows must have both CA_apo and CA_holo (default: %(default)s)."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else project_root / args.targets_csv
    output_image = args.output_image if args.output_image.is_absolute() else project_root / args.output_image
    stats_csv = (
        None
        if args.stats_csv is None
        else (args.stats_csv if args.stats_csv.is_absolute() else project_root / args.stats_csv)
    )

    rc = run_f1_1d_boxplot(
        outputs_dir,
        targets_csv=targets_csv,
        output_image=output_image,
        stats_csv=stats_csv,
        min_ca_coverage=float(args.min_ca_coverage),
    )
    if rc == 0:
        print(f"SI Fig. S14 saved to {output_image.resolve()}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
