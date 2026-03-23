#!/usr/bin/env python3
"""
Update target CSV files with entries from CSP_UBQ.csv.

For scope-based targets (all_alpha, all_beta, alpha_and_beta), filters by
keywords in scope_fold_type. For EC-based targets (hydrolases, isomerases,
oxidoreductases, transferases, translocases), filters by keywords in ec_classes.

Usage:
    python scripts/update_targets_from_csp_ubq.py [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

# Default workspace root (parent of scripts/)
SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent


def load_csp_ubq(csv_path: Path) -> list[dict]:
    """Load CSP_UBQ.csv and return list of row dicts."""
    rows = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            holo = (row.get("holo_pdb") or "").strip()
            if holo:
                rows.append(row)
    return rows


def matches_scope(row: dict, keyword: str) -> bool:
    """Check if scope_fold_type contains the keyword (case-insensitive)."""
    val = (row.get("scope_fold_type") or "").strip()
    return keyword.lower() in val.lower()


def matches_ec(row: dict, keyword: str) -> bool:
    """Check if ec_classes contains the keyword (case-insensitive)."""
    val = (row.get("ec_classes") or "").strip()
    return keyword.lower() in val.lower()


# Target file -> (filter_type, keyword)
TARGET_CONFIG = {
    "targets_all_alpha_proteins.csv": ("scope", "all alpha"),
    "targets_all_beta_proteins.csv": ("scope", "all beta"),
    "targets_alpha_and_beta_proteins_a_plus_b.csv": ("scope", "a+b"),
    "targets_hydrolases.csv": ("ec", "hydrolase"),
    "targets_isomerases.csv": ("ec", "isomerase"),
    "targets_oxidoreductases.csv": ("ec", "oxidoreductase"),
    "targets_transferases.csv": ("ec", "transferase"),
    "targets_translocases.csv": ("ec", "translocase"),
}


def extract_holo_pdbs(rows: list[dict], filter_type: str, keyword: str) -> list[str]:
    """Extract unique holo_pdb values from rows matching the filter."""
    seen: set[str] = set()
    result: list[str] = []
    for row in rows:
        holo = (row.get("holo_pdb") or "").strip()
        if not holo:
            continue
        if filter_type == "scope":
            ok = matches_scope(row, keyword)
        else:
            ok = matches_ec(row, keyword)
        holo_lower = holo.lower()
        if ok and holo_lower not in seen:
            seen.add(holo_lower)
            result.append(holo)
    return sorted(result, key=lambda x: x.lower())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update target CSVs from CSP_UBQ.csv based on scope_fold_type and ec_classes"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be written without modifying files",
    )
    parser.add_argument(
        "--csp-ubq",
        type=Path,
        default=WORKSPACE_ROOT / "CSP_UBQ.csv",
        help="Path to CSP_UBQ.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        default=WORKSPACE_ROOT,
        help="Workspace root directory (default: %(default)s)",
    )
    args = parser.parse_args()

    if not args.csp_ubq.exists():
        print(f"Error: CSP_UBQ.csv not found: {args.csp_ubq}", file=sys.stderr)
        return 1

    rows = load_csp_ubq(args.csp_ubq)
    if not rows:
        print("No rows in CSP_UBQ.csv", file=sys.stderr)
        return 1

    for filename, (filter_type, keyword) in TARGET_CONFIG.items():
        target_path = args.workspace / filename
        holo_pdbs = extract_holo_pdbs(rows, filter_type, keyword)

        if args.dry_run:
            print(f"[DRY-RUN] {filename}: {len(holo_pdbs)} entries")
            for h in holo_pdbs[:5]:
                print(f"  {h}")
            if len(holo_pdbs) > 5:
                print(f"  ... and {len(holo_pdbs) - 5} more")
            continue

        with open(target_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["holo_pdb"])
            for h in holo_pdbs:
                writer.writerow([h])
        print(f"Updated {filename}: {len(holo_pdbs)} entries")

    return 0


if __name__ == "__main__":
    sys.exit(main())
