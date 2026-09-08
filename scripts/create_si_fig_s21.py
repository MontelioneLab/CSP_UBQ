#!/usr/bin/env python3
"""
SI Fig. S21 — Figure 3-style CA-distance histograms under alternate CSP thresholds.

Thin wrapper around ``create_fig_3_thresholds.py`` with SI defaults
(``data/CSP_UBQ_ph0.5_temp5C.csv`` → ``figures/SF21_figure_3_thresholds.png``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts import create_fig_3_thresholds as _impl  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    # Rebuild argv so create_fig_3_thresholds sees SI defaults when omitted.
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs-dir", type=Path, default=_REPO / "outputs")
    ap.add_argument("--figures-dir", type=Path, default=_REPO / "figures")
    ap.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO / "data" / "CSP_UBQ_ph0.5_temp5C.csv",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=_REPO / "figures" / "SF21_figure_3_thresholds.png",
    )
    # Forward unknown flags to the underlying script.
    args, rest = ap.parse_known_args(argv)

    outputs = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO / args.outputs_dir
    figures = args.figures_dir if args.figures_dir.is_absolute() else _REPO / args.figures_dir
    targets = args.targets_csv if args.targets_csv.is_absolute() else _REPO / args.targets_csv
    out = args.output if args.output.is_absolute() else _REPO / args.output

    forwarded = [
        "--outputs-dir",
        str(outputs),
        "--figures-dir",
        str(figures),
        "--targets-csv",
        str(targets),
        "--output",
        str(out),
        *rest,
    ]
    return int(_impl.main(forwarded))


if __name__ == "__main__":
    raise SystemExit(main())
