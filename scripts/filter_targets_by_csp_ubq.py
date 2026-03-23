#!/usr/bin/env python3
"""
Filter targets_*.csv files: remove any row whose holo_pdb is not found in CSP_UBQ.csv.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent


def load_csp_ubq_holo_pdbs(csp_ubq_path: Path) -> set[str]:
    """Load all holo_pdb values from CSP_UBQ.csv (normalized to lowercase)."""
    pdbs: set[str] = set()
    with open(csp_ubq_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "holo_pdb" not in (reader.fieldnames or []):
            return pdbs
        for row in reader:
            val = (row.get("holo_pdb") or "").strip().lower()
            if val:
                pdbs.add(val)
    return pdbs


def main() -> int:
    parser = argparse.ArgumentParser(description="Remove targets not in CSP_UBQ.csv from targets_*.csv")
    parser.add_argument("--csp-ubq", type=Path, default=Path("data/CSP_UBQ.csv"))
    parser.add_argument("--targets-dir", type=Path, default=Path("data"))
    parser.add_argument("--dry-run", action="store_true", help="Report changes without writing")
    args = parser.parse_args()

    targets_dir = args.targets_dir if args.targets_dir.is_absolute() else (_REPO / args.targets_dir)
    csp_path = args.csp_ubq if args.csp_ubq.is_absolute() else (_REPO / args.csp_ubq)
    if not csp_path.exists():
        print(f"Error: {csp_path} not found", file=sys.stderr)
        return 1

    valid_pdbs = load_csp_ubq_holo_pdbs(csp_path)
    print(f"Loaded {len(valid_pdbs)} holo_pdb values from {csp_path}")

    targets_files = sorted(targets_dir.glob("targets_*.csv"))
    if not targets_files:
        print("No targets_*.csv files found", file=sys.stderr)
        return 0

    total_removed = 0
    for path in targets_files:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            if not fieldnames or "holo_pdb" not in fieldnames:
                continue
            rows = list(reader)

        before = len(rows)
        kept = [r for r in rows if (r.get("holo_pdb") or "").strip().lower() in valid_pdbs]
        removed = before - len(kept)

        if removed > 0:
            total_removed += removed
            if args.dry_run:
                removed_ids = [
                    (r.get("holo_pdb") or "").strip()
                    for r in rows
                    if (r.get("holo_pdb") or "").strip().lower() not in valid_pdbs
                ]
                print(f"[DRY-RUN] {path.name}: would remove {removed} row(s): {removed_ids}")
            else:
                with open(path, "w", newline="", encoding="utf-8") as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(kept)
                print(f"{path.name}: removed {removed} row(s)")

    if total_removed == 0 and not args.dry_run:
        print("No entries to remove.")
    elif args.dry_run and total_removed > 0:
        print(f"[DRY-RUN] Would remove {total_removed} total row(s). Run without --dry-run to apply.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
