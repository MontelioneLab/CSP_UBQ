#!/usr/bin/env python3
"""
Annotate CSP_UBQ-style CSV files with per-holo-PDB metadata columns:
- ec_classes
- scope_fold_type

SCOPe fold types are resolved using logic from scripts/get_scope_fold_types.py.
EC classes are resolved via PDBe UniProt mappings and UniProt EC annotations.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set

import requests

try:
    from .config import paths
    from .get_scope_fold_types import (
        download_scope_data,
        fold_types_for,
        parse_scope_folds_by_pdb_chain,
    )
except Exception:
    project_root = str(Path(__file__).resolve().parents[1])
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    from scripts.config import paths
    from scripts.get_scope_fold_types import (
        download_scope_data,
        fold_types_for,
        parse_scope_folds_by_pdb_chain,
    )

EC_CLASS_MAP = {
    "1": "oxidoreductases",
    "2": "transferases",
    "3": "hydrolases",
    "4": "lyases",
    "5": "isomerases",
    "6": "ligases",
    "7": "translocases",
}


def _get_ec_classes_for_pdb(pdb_id: str, timeout_s: float = 30.0) -> str:
    """
    Given a PDB ID, retrieve top-level EC class names from UniProt mappings.
    Returns a comma-separated list, `"Non-enzyme protein"` (UniProt lookups
    succeeded and no EC numbers were found), or `"no EC hit"` (EC class
    could not be resolved for any reason: empty/missing PDB, failed PDBe
    request, or no UniProt mapping).
    """
    pdb_id = (pdb_id or "").strip().lower()
    if not pdb_id:
        return "no EC hit"

    pdbe_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    try:
        pdbe_response = requests.get(pdbe_url, timeout=timeout_s)
        pdbe_response.raise_for_status()
        pdbe_data = pdbe_response.json()
    except requests.RequestException:
        return "no EC hit"

    if pdb_id not in pdbe_data or "UniProt" not in pdbe_data[pdb_id]:
        return "no EC hit"

    uniprot_ids = list(pdbe_data[pdb_id]["UniProt"].keys())
    all_classes: Set[str] = set()

    for uniprot in uniprot_ids:
        uni_url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
        try:
            uni_response = requests.get(uni_url, timeout=timeout_s)
            uni_response.raise_for_status()
            uni_data = uni_response.json()
        except requests.RequestException:
            continue

        ec_numbers: List[str] = []
        protein_description = uni_data.get("proteinDescription", {})

        recommended_name = protein_description.get("recommendedName", {})
        if "ecNumbers" in recommended_name:
            ec_numbers.extend(ec.get("value", "") for ec in recommended_name["ecNumbers"])

        for alt in protein_description.get("alternativeNames", []):
            if "ecNumbers" in alt:
                ec_numbers.extend(ec.get("value", "") for ec in alt["ecNumbers"])

        for sub in protein_description.get("submissionNames", []):
            if "ecNumbers" in sub:
                ec_numbers.extend(ec.get("value", "") for ec in sub["ecNumbers"])

        for ec in set(ec_numbers):
            if not ec:
                continue
            top_level = ec.split(".")[0]
            class_name = EC_CLASS_MAP.get(top_level)
            if class_name:
                all_classes.add(class_name)

    if not all_classes:
        return "Non-enzyme protein"

    return ", ".join(sorted(all_classes))


def _read_csv_rows(csv_path: Path) -> tuple[List[Dict[str, str]], List[str]]:
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def _unique_holo_pdbs(rows: Iterable[Dict[str, str]]) -> List[str]:
    unique: List[str] = []
    for row in rows:
        pdb = (row.get("holo_pdb") or "").strip().lower()
        if pdb and pdb not in unique:
            unique.append(pdb)
    return unique


def _load_receptor_chain_map(csv_path: Path) -> Dict[tuple, str]:
    """
    Load receptor chain assignments from `data/receptor_chain_id.csv`.

    Returns a mapping keyed by (apo_bmrb, holo_bmrb, holo_pdb_lower) to the
    receptor chain letter (lowercase). Missing or blank rows are skipped.
    """
    mapping: Dict[tuple, str] = {}
    if not csv_path.exists():
        return mapping
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            apo_bmrb = (row.get("apo_bmrb") or "").strip()
            holo_bmrb = (row.get("holo_bmrb") or "").strip()
            holo_pdb = (row.get("holo_pdb") or "").strip().lower()
            chain = (row.get("receptor_chain_id") or "").strip().lower()
            if not (apo_bmrb and holo_bmrb and holo_pdb and chain):
                continue
            mapping[(apo_bmrb, holo_bmrb, holo_pdb)] = chain
    return mapping


def annotate_csv_with_ec_and_scope(
    csv_path: Path,
    *,
    verbose: bool = False,
    receptor_chain_csv: Path | None = None,
) -> int:
    """
    Add or update `ec_classes` and `scope_fold_type` columns in a CSP_UBQ-style CSV.

    SCOPe fold types are resolved per row using the receptor chain recorded
    in `receptor_chain_csv` (default: `paths.receptor_chain_csv`). EC classes
    are still resolved per unique `holo_pdb`.

    Returns number of rows written.
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    rows, fieldnames = _read_csv_rows(csv_path)
    if not rows:
        return 0

    chain_csv = Path(receptor_chain_csv) if receptor_chain_csv else Path(paths.receptor_chain_csv)
    receptor_chains = _load_receptor_chain_map(chain_csv)
    if verbose:
        print(f"[META] Loaded {len(receptor_chains)} receptor chain entries from {chain_csv}")

    holo_pdbs = _unique_holo_pdbs(rows)
    if verbose:
        print(f"[META] Found {len(holo_pdbs)} unique holo_pdb IDs")

    if verbose:
        print("[META] Downloading SCOPe classification file...")
    scope_data = download_scope_data()
    scope_lookup = parse_scope_folds_by_pdb_chain(scope_data)

    ec_by_pdb: Dict[str, str] = {pdb: _get_ec_classes_for_pdb(pdb) for pdb in holo_pdbs}

    if "ec_classes" not in fieldnames:
        fieldnames.append("ec_classes")
    if "scope_fold_type" not in fieldnames:
        fieldnames.append("scope_fold_type")

    matched_scope = 0
    missing_chain = 0
    for row in rows:
        pdb = (row.get("holo_pdb") or "").strip().lower()
        apo_bmrb = (row.get("apo_bmrb") or "").strip()
        holo_bmrb = (row.get("holo_bmrb") or "").strip()
        chain = receptor_chains.get((apo_bmrb, holo_bmrb, pdb)) if pdb else None

        if not chain:
            missing_chain += 1
            row["scope_fold_type"] = "No receptor chain mapping"
        else:
            folds = fold_types_for(pdb, chain, scope_lookup)
            if folds:
                matched_scope += 1
                row["scope_fold_type"] = ", ".join(sorted(folds))
            else:
                row["scope_fold_type"] = (
                    "No SCOPe classification found (not present in SCOPe release)"
                )

        row["ec_classes"] = ec_by_pdb.get(pdb, "no EC hit")

    if verbose:
        print(
            "[META] SCOPe coverage: "
            f"{matched_scope}/{len(rows)} rows matched by (holo_pdb, chain); "
            f"{missing_chain} rows had no receptor chain mapping."
        )

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if verbose:
        print(f"[META] Updated {len(rows)} rows in {csv_path}")
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Annotate CSP CSV with ec_classes and scope_fold_type columns."
    )
    parser.add_argument(
        "--input",
        default=paths.input_csv,
        help="Path to CSV file (default: data/CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--receptor-chain-csv",
        default=paths.receptor_chain_csv,
        help="Path to receptor chain CSV (default: data/receptor_chain_id.csv).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress logging.",
    )
    args = parser.parse_args()

    csv_path = Path(args.input)
    try:
        annotate_csv_with_ec_and_scope(
            csv_path,
            verbose=not args.quiet,
            receptor_chain_csv=Path(args.receptor_chain_csv),
        )
    except Exception as exc:
        print(f"[META] ERROR: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

