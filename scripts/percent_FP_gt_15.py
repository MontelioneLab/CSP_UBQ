#!/usr/bin/env python3
"""
Compute the percent of all FP residues whose smallest inter-chain CA-CA distance
is greater than 15 Angstroms, and the percent whose minimum inter-chain atomic
distance (any atom) is greater than 10 Angstroms.

FP = False Positive: significant CSP but not in binding site (classification == 'FP').
- CA-CA: min distance from receptor CA to any ligand CA (min_ca_distance_distance).
- Any-atom: min distance from any receptor atom to any ligand atom (min_any_atom_distance_any_atom).

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
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
ANY_ATOM_DISTANCE_COLUMN = "min_any_atom_distance_any_atom"
CA_DISTANCE_THRESHOLD = 15.0  # Angstroms
ANY_ATOM_DISTANCE_THRESHOLD = 10.0  # Angstroms


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
        description="Percent of FP residues with min inter-chain distances above thresholds"
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
    parser.add_argument(
        "--ca-threshold",
        type=float,
        default=CA_DISTANCE_THRESHOLD,
        help=f"CA-CA distance threshold in Angstroms (default: {CA_DISTANCE_THRESHOLD})",
    )
    parser.add_argument(
        "--any-atom-threshold",
        type=float,
        default=ANY_ATOM_DISTANCE_THRESHOLD,
        help=f"Any-atom distance threshold in Angstroms (default: {ANY_ATOM_DISTANCE_THRESHOLD})",
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

    total_fp_ca = 0
    fp_gt_ca_threshold = 0
    total_fp_any = 0
    fp_gt_any_threshold = 0

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

        fp_mask = df[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper() == "FP"
        fp_df = df.loc[fp_mask]

        # CA-CA distance
        if CA_DISTANCE_COLUMN in df.columns:
            ca_valid = fp_df[CA_DISTANCE_COLUMN].notna()
            fp_with_ca = fp_df.loc[ca_valid]
            n_fp_ca = len(fp_with_ca)
            if n_fp_ca > 0:
                n_gt_ca = (fp_with_ca[CA_DISTANCE_COLUMN].astype(float) > args.ca_threshold).sum()
                total_fp_ca += n_fp_ca
                fp_gt_ca_threshold += n_gt_ca

        # Any-atom distance
        if ANY_ATOM_DISTANCE_COLUMN in df.columns:
            any_valid = fp_df[ANY_ATOM_DISTANCE_COLUMN].notna()
            fp_with_any = fp_df.loc[any_valid]
            n_fp_any = len(fp_with_any)
            if n_fp_any > 0:
                n_gt_any = (fp_with_any[ANY_ATOM_DISTANCE_COLUMN].astype(float) > args.any_atom_threshold).sum()
                total_fp_any += n_fp_any
                fp_gt_any_threshold += n_gt_any

    # Report CA-CA
    if total_fp_ca > 0:
        pct_ca = 100.0 * fp_gt_ca_threshold / total_fp_ca
        print(f"--- CA-CA distance ---")
        print(f"Total FP residues (with CA distance): {total_fp_ca}")
        print(f"FP residues with min CA-CA distance > {args.ca_threshold} Å: {fp_gt_ca_threshold}")
        print(f"Percent: {pct_ca:.2f}%")
    else:
        print("No FP residues with CA distance data found.")

    # Report any-atom
    if total_fp_any > 0:
        pct_any = 100.0 * fp_gt_any_threshold / total_fp_any
        print(f"\n--- Any-atom distance ---")
        print(f"Total FP residues (with any-atom distance): {total_fp_any}")
        print(f"FP residues with min inter-chain atomic distance > {args.any_atom_threshold} Å: {fp_gt_any_threshold}")
        print(f"Percent: {pct_any:.2f}%")
    else:
        print("\nNo FP residues with any-atom distance data found.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
