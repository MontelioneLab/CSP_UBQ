#!/usr/bin/env python3
"""Backfill DSSP secondary-structure labels into outputs/*/master_alignment.csv.

Writes a new ``ss`` column with values H / E / C (helix / sheet / coil-loop)
from mkdssp via BioPython. Residues without a DSSP assignment are left empty.

Master ``pdb_residue_number`` is often the CSP sequential index, not the
authentic PDB residue number. This script maps rows to PDB residue numbers by
aligning the master sequence to occlusion / filter CSVs (same approach as
visualize.py), then looks up DSSP labels.
"""

from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.align import align_global
from scripts.secondary_structure import compute_dssp_secondary_structure

SS_COLUMN = "ss"

AA3_TO1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

# Sibling tables with authentic PDB residue_number + residue_name.
_PDB_SEQ_SOURCES = (
    "occlusion_analysis.csv",
    "ca_distance_filter.csv",
    "any_atom_distance_filter.csv",
    "nn_distance_filter.csv",
    "interaction_filter.csv",
)


def _discover_targets(
    outputs_dir: Path,
    ids: Optional[List[str]] = None,
    limit: Optional[int] = None,
) -> List[Path]:
    if not outputs_dir.is_dir():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    targets = sorted(
        p
        for p in outputs_dir.iterdir()
        if p.is_dir()
        and not p.name.startswith(".")
        and (p / "master_alignment.csv").is_file()
    )
    if ids:
        wanted = {x.strip() for x in ids if x and x.strip()}
        targets = [
            t
            for t in targets
            if t.name in wanted
            or t.name.split("_", 1)[0].upper() in {w.upper() for w in wanted}
        ]
    if limit is not None:
        targets = targets[: max(0, int(limit))]
    return targets


def _resolve_pdb_path(holo_pdb: str, pdb_dir: Path) -> Optional[Path]:
    pid = (holo_pdb or "").strip()
    if not pid:
        return None
    for name in (f"{pid.lower()}.pdb", f"{pid.upper()}.pdb", f"{pid}.pdb"):
        cand = pdb_dir / name
        if cand.is_file():
            return cand
    return None


def _parse_int(raw: Any) -> Optional[int]:
    if raw in ("", None):
        return None
    try:
        return int(float(str(raw).strip()))
    except (ValueError, TypeError):
        return None


def _load_pdb_sequence(
    target_dir: Path,
) -> Tuple[str, List[int]]:
    """Load authentic PDB residue sequence from a sibling filter CSV."""
    for fname in _PDB_SEQ_SOURCES:
        path = target_dir / fname
        if not path.is_file():
            continue
        seq: List[str] = []
        positions: List[int] = []
        with path.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                name = (row.get("residue_name") or "").strip().upper()
                resi = _parse_int(row.get("residue_number"))
                if resi is None or not name or name == "PRO":
                    continue
                aa = AA3_TO1.get(name)
                if not aa:
                    continue
                seq.append(aa)
                positions.append(resi)
        if seq:
            return "".join(seq), positions
    return "", []


def _master_sequence(
    rows: Sequence[Dict[str, str]],
) -> Tuple[str, List[int]]:
    """CSP-side sequence keyed by row index into ``rows`` (prolines skipped)."""
    seq: List[str] = []
    row_indices: List[int] = []
    for i, row in enumerate(rows):
        aa = (row.get("holo_aa") or "").strip().upper()
        if not aa or aa == "P":
            continue
        seq.append(aa)
        row_indices.append(i)
    return "".join(seq), row_indices


def _map_rows_to_pdb_resnums(
    rows: Sequence[Dict[str, str]],
    target_dir: Path,
) -> Dict[int, int]:
    """
    Map master row index -> authentic PDB residue number.

    Prefer sequence alignment to occlusion/filter CSVs. Fall back to
    ``pdb_residue_number`` / ``holo_resi`` when no sibling sequence is available
    or when those values already match DSSP keys (identity numbering).
    """
    csp_seq, csp_row_idx = _master_sequence(rows)
    pdb_seq, pdb_positions = _load_pdb_sequence(target_dir)

    row_to_pdb: Dict[int, int] = {}
    if csp_seq and pdb_seq:
        _a, _b, mapping, _score = align_global(csp_seq, pdb_seq)
        for csp_pos, pdb_pos in mapping:
            # align_global returns 1-based indices into the ungapped sequences
            if 1 <= csp_pos <= len(csp_row_idx) and 1 <= pdb_pos <= len(pdb_positions):
                row_to_pdb[csp_row_idx[csp_pos - 1]] = pdb_positions[pdb_pos - 1]

    # Fill any unmapped rows (including prolines) from stored residue columns.
    for i, row in enumerate(rows):
        if i in row_to_pdb:
            continue
        for key in ("pdb_residue_number", "holo_resi"):
            resi = _parse_int(row.get(key))
            if resi is not None:
                row_to_pdb[i] = resi
                break
    return row_to_pdb


def process_target(
    target_dir: str,
    pdb_dir: str,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """Backfill ``ss`` for one target directory. Safe for ProcessPoolExecutor."""
    tgt = Path(target_dir)
    name = tgt.name
    master = tgt / "master_alignment.csv"
    result: Dict[str, Any] = {
        "name": name,
        "ok": False,
        "n_rows": 0,
        "filled": 0,
        "skipped": False,
        "error": "",
    }
    try:
        with master.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
        if not fieldnames:
            result["error"] = "empty or missing header"
            return result
        result["n_rows"] = len(rows)
        if SS_COLUMN in fieldnames and not overwrite:
            n_filled = sum(1 for r in rows if (r.get(SS_COLUMN) or "").strip())
            if n_filled > 0:
                result["ok"] = True
                result["filled"] = n_filled
                result["skipped"] = True
                return result

        holo_pdb = ""
        chain = ""
        if rows:
            holo_pdb = (rows[0].get("holo_pdb") or "").strip() or name.split("_", 1)[0]
            chain = (
                (rows[0].get("chain") or "").strip()
                or (rows[0].get("chain_id") or "").strip()
                or (rows[0].get("chain_id_occlusion") or "").strip()
            )
        else:
            holo_pdb = name.split("_", 1)[0]

        pdb_path = _resolve_pdb_path(holo_pdb, Path(pdb_dir))
        if pdb_path is None:
            result["error"] = f"PDB not found for {holo_pdb!r}"
            return result

        ss_map = compute_dssp_secondary_structure(str(pdb_path), chain or None)
        if not ss_map:
            result["error"] = f"DSSP returned no labels for {pdb_path.name}"
            return result

        row_to_pdb = _map_rows_to_pdb_resnums(rows, tgt)

        filled = 0
        for i, row in enumerate(rows):
            resi = row_to_pdb.get(i)
            label = ss_map.get(resi, "") if resi is not None else ""
            # If identity numbering failed but authentic map missed this row,
            # try direct lookup of stored pdb_residue_number against DSSP.
            if not label:
                for key in ("pdb_residue_number", "holo_resi"):
                    alt = _parse_int(row.get(key))
                    if alt is not None and alt in ss_map:
                        label = ss_map[alt]
                        break
            row[SS_COLUMN] = label
            if label:
                filled += 1

        if SS_COLUMN not in fieldnames:
            fieldnames = fieldnames + [SS_COLUMN]

        tmp = master.with_suffix(".csv.tmp")
        with tmp.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(master)

        result["ok"] = True
        result["filled"] = filled
        return result
    except Exception as exc:
        result["error"] = str(exc)
        return result


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_ROOT / "outputs",
        help="Root with per-target master_alignment.csv "
        "(default: outputs/)",
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=_ROOT / "PDB_FILES",
        help="Directory of PDB files (default: PDB_FILES/)",
    )
    parser.add_argument("--ids", nargs="*", help="Optional target dir names to process")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N targets")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers (default 1)")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute ss even if the column is already populated",
    )
    args = parser.parse_args(argv)

    targets = _discover_targets(args.outputs, args.ids, args.limit)
    print(f"Found {len(targets)} targets under {args.outputs}", flush=True)
    if not targets:
        return 0

    ok = fail = skipped = 0
    failures: List[Tuple[str, str]] = []
    workers = max(1, int(args.workers))

    def _handle(result: Dict[str, Any]) -> None:
        nonlocal ok, fail, skipped
        if result.get("ok"):
            ok += 1
            if result.get("skipped"):
                skipped += 1
            else:
                print(
                    f"OK {result['name']}: filled {result['filled']}/{result['n_rows']}",
                    flush=True,
                )
        else:
            fail += 1
            failures.append((result["name"], str(result.get("error"))))
            print(f"FAIL {result['name']}: {result.get('error')}", flush=True)

    if workers == 1:
        for i, tgt in enumerate(targets, 1):
            _handle(process_target(str(tgt), str(args.pdb_dir), args.overwrite))
            if i % 25 == 0 or i == len(targets):
                print(
                    f"[{i}/{len(targets)}] ok={ok} fail={fail} skipped={skipped}",
                    flush=True,
                )
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {
                ex.submit(process_target, str(t), str(args.pdb_dir), args.overwrite): t
                for t in targets
            }
            for i, fut in enumerate(as_completed(futs), 1):
                _handle(fut.result())
                if i % 25 == 0 or i == len(targets):
                    print(
                        f"[{i}/{len(targets)}] ok={ok} fail={fail} skipped={skipped}",
                        flush=True,
                    )

    print(
        f"Done: ok={ok} fail={fail} skipped_existing={skipped} total={len(targets)}",
        flush=True,
    )
    if failures:
        print("Failures:", flush=True)
        for name, err in failures[:20]:
            print(f"  {name}: {err}", flush=True)
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more", flush=True)
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
