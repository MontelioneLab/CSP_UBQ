#!/usr/bin/env python3
"""
Compute the average % of FP residues across all protein receptors, and report
the minimum % of FP residues across all receptors.

For each receptor: FP % = (number of FP residues / total classified residues) * 100
FP = False Positive: significant CSP but not in binding site (classification == 'FP').

Also reports the mean of FP/(TP+FP) over receptors with at least one significant
residue (TP+FP > 0), and the mean per-receptor F1 from the same definition as
``f1_score_reporter.py`` (``analyze_targets.compute_f1_score``).

By default only targets listed in the given ``--targets-csv`` are included. Each row is
resolved to the pipeline output folder using ``apo_bmrb``, ``holo_bmrb``, and ``holo_pdb``
congruence with the first data row of ``outputs/<dir>/master_alignment.csv`` (including
``{holo_pdb}_1``, ``_2``, … candidates); when multiple dirs match the same BMRB pair, the
first by suffix order is used — same rule as ``scripts.target_resolution`` /
figure-creation scripts.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional

try:
    from .analyze_targets import (
        AlignmentParsingError,
        compute_f1_score,
        load_alignment,
        PREDICTOR_COLUMNS,
    )
    from .target_resolution import load_target_rows, resolve_target_rows
except ImportError:
    from analyze_targets import (
        AlignmentParsingError,
        compute_f1_score,
        load_alignment,
        PREDICTOR_COLUMNS,
    )
    from target_resolution import load_target_rows, resolve_target_rows

CLASSIFICATION_COLUMN = "classification"
VALID_CLASSIFICATIONS = frozenset({"TP", "FP", "TN", "FN"})

_REQUIRED_TARGET_COLS = frozenset({"apo_bmrb", "holo_bmrb", "holo_pdb"})


def _validate_target_columns(fieldnames: Optional[List[str]], csv_path: Path) -> None:
    if not fieldnames:
        raise ValueError(f"{csv_path} has no header row")
    missing = _REQUIRED_TARGET_COLS - {fn.lower() for fn in fieldnames}
    if missing:
        raise ValueError(
            f"{csv_path} must contain columns {sorted(_REQUIRED_TARGET_COLS)}; missing {sorted(missing)}"
        )


def _try_append_receptor_metrics(
    target_dir: Path,
    fp_percents: List[float],
    target_to_pct: List[tuple[str, float, int, int]],
    fp_over_tp_fp: List[float],
    f1_scores: List[float],
) -> bool:
    """Parse ``master_alignment.csv`` and append aggregate metrics; return True if appended."""
    alignment_path = target_dir / "master_alignment.csv"
    if not alignment_path.exists():
        return False
    try:
        df = load_alignment(alignment_path)
    except AlignmentParsingError as exc:
        print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
        return False
    if CLASSIFICATION_COLUMN not in df.columns:
        return False
    cls = df[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper()
    valid_mask = cls.isin(VALID_CLASSIFICATIONS)
    classified = df.loc[valid_mask]
    total = len(classified)
    if total == 0:
        return False
    cls_upper = classified[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper()
    n_fp = int((cls_upper == "FP").sum())
    n_tp = int((cls_upper == "TP").sum())
    pct = 100.0 * n_fp / total
    fp_percents.append(pct)
    target_to_pct.append((target_dir.name, pct, n_fp, total))
    sig_denom = n_tp + n_fp
    if sig_denom > 0:
        fp_over_tp_fp.append(n_fp / sig_denom)
    predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)
    f1_scores.append(compute_f1_score(df, predicted).f1)
    return True


def discover_targets(outputs_dir: Path) -> list[Path]:
    """Return sorted list of target directories under outputs_dir."""
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    return sorted(
        p for p in outputs_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Average and minimum % of FP residues; mean FP/(TP+FP) and mean F1 "
            "(same F1 as f1_score_reporter.py)."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: outputs)",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help=(
            "CSV with apo_bmrb, holo_bmrb, holo_pdb (one output folder per row). "
            "Default: data/CSP_UBQ_ph0.5_temp5C.csv. Ignored with --all-targets."
        ),
    )
    parser.add_argument(
        "--all-targets",
        action="store_true",
        help="Use every subdirectory under --outputs-dir (ignore --targets-csv).",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else repo_root / args.outputs_dir

    fp_percents: list[float] = []
    target_to_pct: list[tuple[str, float, int, int]] = []  # (target, pct, fp_count, total)
    fp_over_tp_fp: list[float] = []
    f1_scores: list[float] = []

    if args.all_targets:
        print("Using all subdirectories under outputs (--all-targets).")
        for target_dir in discover_targets(outputs_dir):
            _try_append_receptor_metrics(
                target_dir, fp_percents, target_to_pct, fp_over_tp_fp, f1_scores
            )
    else:
        targets_csv = args.targets_csv if args.targets_csv.is_absolute() else repo_root / args.targets_csv
        if not targets_csv.is_file():
            print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
            return 1

        try:
            with targets_csv.open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                _validate_target_columns(reader.fieldnames, targets_csv)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        rows = load_target_rows(targets_csv)
        paths = resolve_target_rows(rows, outputs_dir)
        print(
            f"Resolving {len(rows)} rows from {targets_csv.name} to {len(paths)} output dirs "
            "(BMRB-congruent master_alignment match; see scripts/target_resolution)"
        )

        for path in paths:
            _try_append_receptor_metrics(
                path, fp_percents, target_to_pct, fp_over_tp_fp, f1_scores
            )

    if not fp_percents:
        print("No receptors with classification data found.")
        return 0

    avg_pct = sum(fp_percents) / len(fp_percents)
    min_pct = min(fp_percents)
    min_target = min(target_to_pct, key=lambda x: x[1])

    print(f"Receptors analyzed: {len(fp_percents)}")
    print(f"Average % of FP residues: {avg_pct:.2f}%")
    print(f"Minimum % of FP residues: {min_pct:.2f}% (receptor: {min_target[0]}, {min_target[2]} FP / {min_target[3]} total)")
    if fp_over_tp_fp:
        mean_fp_ratio = sum(fp_over_tp_fp) / len(fp_over_tp_fp)
        print(f"Mean FP / (TP+FP): {mean_fp_ratio:.3f} (over {len(fp_over_tp_fp)} receptors with TP+FP > 0)")
    else:
        print("Mean FP / (TP+FP): n/a (no receptors with TP+FP > 0)")
    if f1_scores:
        mean_f1 = sum(f1_scores) / len(f1_scores)
        print(f"Mean F1 score: {mean_f1:.3f} (over {len(f1_scores)} receptors)")
    else:
        print("Mean F1 score: n/a")

    return 0


if __name__ == "__main__":
    sys.exit(main())
