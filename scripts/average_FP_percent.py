#!/usr/bin/env python3
"""
Compute the average % of FP residues across all protein receptors, and report
the minimum % of FP residues across all receptors.

For each receptor: FP % = (number of FP residues / total classified residues) * 100
FP = False Positive: significant CSP but not in binding site (classification == 'FP').

By default only targets listed in CSP_UBQ_ph0.5_temp5C.csv (holo_pdb) are included.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional, Set

import pandas as pd

CLASSIFICATION_COLUMN = "classification"
VALID_CLASSIFICATIONS = frozenset({"TP", "FP", "TN", "FN"})


def discover_targets(outputs_dir: Path) -> list[Path]:
    """Return sorted list of target directories under outputs_dir."""
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    return sorted(
        p for p in outputs_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )


def load_holo_pdb_set(targets_csv: Path) -> Set[str]:
    """Lowercase holo_pdb IDs from CSV (holo_pdb column)."""
    out: Set[str] = set()
    with targets_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "holo_pdb" not in reader.fieldnames:
            raise ValueError(f"{targets_csv} must contain a 'holo_pdb' column")
        for row in reader:
            h = (row.get("holo_pdb") or "").strip().lower()
            if h:
                out.add(h)
    return out


def output_dir_matches_holo_set(dir_name: str, holo_lower: Set[str]) -> bool:
    """True if directory name (or base before '_') is in the allowed holo set."""
    t = dir_name.lower()
    return t in holo_lower or t.split("_", 1)[0] in holo_lower


def load_alignment(path: Path) -> pd.DataFrame:
    """Load master_alignment.csv as DataFrame."""
    return pd.read_csv(path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Average and minimum % of FP residues across protein receptors"
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
            "CSV with holo_pdb column; only matching output subdirectories are used "
            "(default: data/CSP_UBQ_ph0.5_temp5C.csv). Ignored with --all-targets."
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

    holo_filter: Optional[Set[str]] = None
    if not args.all_targets:
        targets_csv = args.targets_csv if args.targets_csv.is_absolute() else repo_root / args.targets_csv
        if not targets_csv.is_file():
            print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
            return 1
        holo_filter = load_holo_pdb_set(targets_csv)
        print(f"Filtering to {len(holo_filter)} holo_pdb IDs from {targets_csv.name}")
    else:
        print("Using all subdirectories under outputs (--all-targets).")

    fp_percents: list[float] = []
    target_to_pct: list[tuple[str, float, int, int]] = []  # (target, pct, fp_count, total)

    for target_dir in discover_targets(outputs_dir):
        if holo_filter is not None and not output_dir_matches_holo_set(target_dir.name, holo_filter):
            continue
        alignment_path = target_dir / "master_alignment.csv"
        if not alignment_path.exists():
            continue

        try:
            df = load_alignment(alignment_path)
        except Exception as exc:
            print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
            continue

        if CLASSIFICATION_COLUMN not in df.columns:
            continue

        cls = df[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper()
        valid_mask = cls.isin(VALID_CLASSIFICATIONS)
        classified = df.loc[valid_mask]

        total = len(classified)
        if total == 0:
            continue

        n_fp = (classified[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper() == "FP").sum()
        pct = 100.0 * n_fp / total
        fp_percents.append(pct)
        target_to_pct.append((target_dir.name, pct, n_fp, total))

    if not fp_percents:
        print("No receptors with classification data found.")
        return 0

    avg_pct = sum(fp_percents) / len(fp_percents)
    min_pct = min(fp_percents)
    min_target = min(target_to_pct, key=lambda x: x[1])

    print(f"Receptors analyzed: {len(fp_percents)}")
    print(f"Average % of FP residues: {avg_pct:.2f}%")
    print(f"Minimum % of FP residues: {min_pct:.2f}% (receptor: {min_target[0]}, {min_target[2]} FP / {min_target[3]} total)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
