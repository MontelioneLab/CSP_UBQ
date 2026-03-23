#!/usr/bin/env python3
"""
Get SCOPe fold type classification for each holo PDB ID in CSP_UBQ.csv.

Downloads the SCOPe classification file once, then looks up fold types for all
unique holo_pdb values. Writes results to a CSV (default: outputs/holo_pdb_scope_fold_types.csv).
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import requests

try:
    from .config import paths
except Exception:
    Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.config import paths


SCOPE_URL = "https://zenodo.org/records/5829561/files/dir.cla.scope.2.08-stable.txt?download=1"

CLASS_MAP = {
    "a": "all alpha proteins",
    "b": "all beta proteins",
    "c": "alpha and beta proteins (a/b)",
    "d": "alpha and beta proteins (a+b)",
    "e": "multi-domain proteins (alpha and beta)",
    "f": "membrane and cell surface proteins",
    "g": "small proteins",
    "h": "coiled coil proteins",
    "i": "low resolution protein structures",
    "j": "peptides",
    "k": "designed proteins",
    "l": "artifacts",
}

def download_scope_data(url: str = SCOPE_URL) -> str:
    """Download SCOPe classification file and return its text."""
    response = requests.get(url, timeout=120)
    if not response.ok:
        raise ValueError(f"Failed to download SCOPe classification file: {response.status_code}")
    return response.text


def parse_scope_folds_by_pdb(data: str) -> dict[str, set[str]]:
    """
    Parse SCOPe dir.cla.scope file and return a mapping:
    pdb_id_lower -> set of fold type strings (e.g. 'all alpha proteins').
    """
    pdb_to_folds: dict[str, set[str]] = {}
    for line in data.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        domain_id = parts[0]
        # domain_id is like d1cf4a_, d1sy9a_, etc. (lowercase PDB)
        if not domain_id.startswith("d") or len(domain_id) < 5:
            continue
        # Extract PDB ID: after 'd', next 4 chars are the PDB ID (e.g. 1cf4, 1sy9)
        pdb_id = domain_id[1:5].lower()
        # Current SCOPe 2.08 format: [sid, pdbid, residues, sccs, sunid, hierarchy]
        # Keep compatibility with older assumptions by falling back to column 2 if needed.
        sccs = ""
        if len(parts) >= 4 and "." in parts[3]:
            sccs = parts[3].strip()
        elif len(parts) >= 2:
            sccs = parts[1].strip()
        if not sccs:
            continue
        class_letter = sccs[0].lower()
        fold_type = CLASS_MAP.get(class_letter, "unknown")
        pdb_to_folds.setdefault(pdb_id, set()).add(fold_type)
    return pdb_to_folds


def get_unique_holo_pdbs(csp_csv_path: Path) -> list[str]:
    """Read CSP_UBQ.csv and return unique, non-empty holo_pdb values (normalized to lowercase)."""
    pdbs = []
    with open(csp_csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            h = (row.get("holo_pdb") or "").strip()
            if not h:
                continue
            # Normalize for lookup (SCOPe uses 4-char lowercase)
            h_lo = h.lower()
            if h_lo not in pdbs:
                pdbs.append(h_lo)
    return pdbs


def main() -> int:
    csp_csv = Path(paths.input_csv)
    out_csv = Path(paths.outputs_dir) / "holo_pdb_scope_fold_types.csv"

    if not csp_csv.exists():
        print(f"Error: {csp_csv} not found.", file=sys.stderr)
        return 1

    print("Loading unique holo PDB IDs from CSP_UBQ.csv...")
    holo_pdbs = get_unique_holo_pdbs(csp_csv)
    print(f"Found {len(holo_pdbs)} unique holo PDB IDs.")

    print("Downloading SCOPe classification file...")
    data = download_scope_data()
    pdb_to_folds = parse_scope_folds_by_pdb(data)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb in holo_pdbs:
        folds = pdb_to_folds.get(pdb, set())
        if not folds:
            fold_str = "No SCOPe classification found"
        else:
            fold_str = ", ".join(sorted(folds))
        rows.append({"holo_pdb": pdb.upper(), "scope_fold_type": fold_str})

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["holo_pdb", "scope_fold_type"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_csv}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
