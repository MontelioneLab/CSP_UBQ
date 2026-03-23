#!/usr/bin/env python3
"""
SI Table S1. Data for all receptors in this study.

Creates a LaTeX table (figures/ST1_all_receptors.tex) with study targets (buffer-
filtered CSP_UBQ_ph0.5_temp5C.csv by default), their apo/holo PDB and BMRB IDs,
F1 and MCC scores, and hyperlinks to RCSB/BMRB where available. Format matches
other supplementary tables (ST2–ST13).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .create_csp_latex_table import build_latex, load_and_merge
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.create_csp_latex_table import build_latex, load_and_merge  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create SI Table S1 (Data for all receptors)."
    )
    parser.add_argument(
        "--csp-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="Path to CSP table (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument(
        "--confusion-csv",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv (default: outputs/confusion_matrix_per_system.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figures") / "ST1_all_receptors.tex",
        help="Output .tex file (default: figures/ST1_all_receptors.tex).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    csp_path = args.csp_csv if args.csp_csv.is_absolute() else project_root / args.csp_csv
    confusion_path = (
        args.confusion_csv
        if args.confusion_csv.is_absolute()
        else project_root / args.confusion_csv
    )
    out_path = args.output if args.output.is_absolute() else project_root / args.output

    if not csp_path.exists():
        print(f"Error: CSP CSV not found: {csp_path}", file=sys.stderr)
        return 1
    if not confusion_path.exists():
        print(f"Error: Confusion CSV not found: {confusion_path}", file=sys.stderr)
        return 1

    df = load_and_merge(
        csp_path,
        confusion_path,
        include_all=True,
        allowed_targets=None,
    )

    if df.empty:
        print("No rows to write.", file=sys.stderr)
        return 1

    latex = build_latex(
        df,
        caption=r"\textbf{SI Table S1. Data for all receptors in this study}",
        label="tab:st1",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(latex, encoding="utf-8")
    print(f"SI Table S1 written to {out_path.resolve()} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
