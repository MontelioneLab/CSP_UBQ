#!/usr/bin/env python3
"""
SI Fig. S23 — Figure 3 histograms for the same-author/same-study subset.

Thin wrapper around ``create_fig_3.py`` with SI defaults
(``data/CSP_UBQ_ph0.5_temp5C_same_author_list.csv`` → SF23 combined PNG).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts import create_fig_3 as _impl  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs-dir", type=Path, default=_REPO / "outputs")
    ap.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO / "data" / "CSP_UBQ_ph0.5_temp5C_same_author_list.csv",
    )
    ap.add_argument(
        "--output-a",
        type=Path,
        default=_REPO / "figures" / "figure_3_a_same_author_list.png",
    )
    ap.add_argument(
        "--output-b",
        type=Path,
        default=_REPO / "figures" / "figure_3_b_same_author_list.png",
    )
    ap.add_argument(
        "--output-combined",
        type=Path,
        default=_REPO / "figures" / "SF23_figure_3_combined_same_author_list.png",
    )
    args, rest = ap.parse_known_args(argv)

    def _abs(p: Path) -> Path:
        return p if p.is_absolute() else _REPO / p

    forwarded = [
        "--outputs-dir",
        str(_abs(args.outputs_dir)),
        "--targets-csv",
        str(_abs(args.targets_csv)),
        "--output-a",
        str(_abs(args.output_a)),
        "--output-b",
        str(_abs(args.output_b)),
        "--output-combined",
        str(_abs(args.output_combined)),
        *rest,
    ]
    return int(_impl.main(forwarded))


if __name__ == "__main__":
    raise SystemExit(main())
