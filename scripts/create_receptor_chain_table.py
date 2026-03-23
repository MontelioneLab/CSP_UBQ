"""
Create receptor_chain_id.csv from pipeline outputs.

Scans outputs/*/master_alignment.csv (or csp_table.csv) to extract apo_bmrb,
holo_bmrb, holo_pdb, and receptor chain for each target. Outputs a table with
5 columns: apo_pdb, holo_pdb, apo_bmrb, holo_bmrb, receptor_chain_id.

Usage:
  python scripts/create_receptor_chain_table.py [--output data/receptor_chain_id.csv] [--outputs-dir outputs]
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def discover_targets(outputs_dir: Path) -> list[Path]:
    """Return sorted list of target directories under outputs_dir."""
    if not outputs_dir.exists():
        return []
    return sorted(
        p for p in outputs_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )


def find_receptor_chain(target_dir: Path) -> tuple[str, str, str, str] | None:
    """
    Extract apo_bmrb, holo_bmrb, holo_pdb, receptor_chain from the first data row
    of master_alignment.csv or csp_table.csv.

    Returns:
        (apo_bmrb, holo_bmrb, holo_pdb, receptor_chain_id) or None if not found.
    """
    for fname in ("master_alignment.csv", "csp_table.csv"):
        csv_path = target_dir / fname
        if not csv_path.exists():
            continue
        try:
            with open(csv_path, "r", newline="") as f:
                reader = csv.DictReader(f)
                row = next(reader, None)
                if not row:
                    continue
                apo_bmrb = (row.get("apo_bmrb") or "").strip()
                holo_bmrb = (row.get("holo_bmrb") or "").strip()
                holo_pdb = (row.get("holo_pdb") or "").strip()
                chain = (row.get("chain") or "").strip()
                if apo_bmrb and holo_bmrb and holo_pdb and chain:
                    return (apo_bmrb, holo_bmrb, holo_pdb, chain)
        except (csv.Error, OSError):
            continue
    return None


def load_apo_pdb_map(input_csv: Path) -> dict[tuple[str, str, str], str]:
    """
    Load CSP_UBQ.csv and build (apo_bmrb, holo_bmrb, holo_pdb) -> apo_pdb.
    """
    mapping: dict[tuple[str, str, str], str] = {}
    if not input_csv.exists():
        return mapping
    try:
        with open(input_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                apo_bmrb = (row.get("apo_bmrb") or "").strip()
                holo_bmrb = (row.get("holo_bmrb") or "").strip()
                holo_pdb = (row.get("holo_pdb") or "").strip()
                apo_pdb = (row.get("apo_pdb") or "").strip()
                if apo_bmrb and holo_bmrb and holo_pdb:
                    key = (apo_bmrb, holo_bmrb, holo_pdb)
                    mapping[key] = apo_pdb
    except (csv.Error, OSError):
        pass
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create receptor_chain_id.csv from pipeline outputs."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/receptor_chain_id.csv"),
        help="Output CSV path (default: data/receptor_chain_id.csv)",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Outputs directory (default: outputs)",
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/CSP_UBQ.csv"),
        help="Input CSV for apo_pdb lookup (default: data/CSP_UBQ.csv)",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    outputs_dir = project_root / args.outputs_dir if not args.outputs_dir.is_absolute() else args.outputs_dir
    input_csv = project_root / args.input_csv if not args.input_csv.is_absolute() else args.input_csv
    output_path = project_root / args.output if not args.output.is_absolute() else args.output

    apo_pdb_map = load_apo_pdb_map(input_csv)

    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()

    for target_dir in discover_targets(outputs_dir):
        result = find_receptor_chain(target_dir)
        if not result:
            continue
        apo_bmrb, holo_bmrb, holo_pdb, receptor_chain_id = result
        key = (apo_bmrb, holo_bmrb, holo_pdb)
        if key in seen:
            continue
        seen.add(key)
        apo_pdb = apo_pdb_map.get(key, "")
        rows.append({
            "apo_pdb": apo_pdb,
            "holo_pdb": holo_pdb,
            "apo_bmrb": apo_bmrb,
            "holo_bmrb": holo_bmrb,
            "receptor_chain_id": receptor_chain_id,
        })

    rows.sort(key=lambda r: (r["apo_bmrb"], r["holo_bmrb"], r["holo_pdb"]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["apo_pdb", "holo_pdb", "apo_bmrb", "holo_bmrb", "receptor_chain_id"],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
