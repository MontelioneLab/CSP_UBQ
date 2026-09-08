#!/usr/bin/env python3
"""Backfill rci_apo / rci_holo into outputs/*/master_alignment.csv via RCI_GUI.

RCI is computed on the same apo/holo saveframes selected by CSP (best global
alignment score). SHIFTY residue numbers use CSP seq_pos indexing so values
join cleanly onto apo_resi / holo_resi.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

RCI_SCRIPT = _ROOT / "RCI_GUI" / "scripts" / "rci_v_1c_PyNMR-STAR.py"
RCI_CACHE_DIR = _ROOT / ".cache" / "rci_shifty"
# GUI defaults minus -mpl/-dynamr. Omit -no_i: the legacy script defaults to
# incomplete_data_use=0 and has no CLI to enable incomplete residues, but most
# CSP BMRB lists lack HA and/or CO, which would yield empty RCI output.
RCI_CLI_FLAGS = ["-r", "4", "-d", "1", "-end_corr2", "-gapfill2"]
_PATCHED_RCI_SCRIPT: Optional[Path] = None


def _patched_rci_script() -> Path:
    """Return a cached copy of the RCI engine with incomplete_data_use=1."""
    global _PATCHED_RCI_SCRIPT
    if _PATCHED_RCI_SCRIPT is not None and _PATCHED_RCI_SCRIPT.is_file():
        return _PATCHED_RCI_SCRIPT
    src = RCI_SCRIPT.read_text(encoding="utf-8", errors="ignore")
    old = (
        "incomplete_data_use=0 # 0 - do not use data for a residue when chemical "
        "shifts are not known for all atoms, 1 - use even if data is incomplete"
    )
    new = (
        "incomplete_data_use=1 # patched by backfill_rci.py for incomplete CSP "
        "shift lists"
    )
    if old not in src:
        raise RuntimeError("Could not patch incomplete_data_use in RCI script")
    RCI_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out = RCI_CACHE_DIR / "rci_v_1c_PyNMR-STAR_incomplete.py"
    out.write_text(src.replace(old, new, 1), encoding="utf-8")
    _PATCHED_RCI_SCRIPT = out
    return out

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
AMIDE_H = ("H", "HN", "H1")
ALPHA_H = ("HA", "HA2", "HA3")
CO_NAMES = ("C", "CO")


def _median(atoms: Dict[str, List[float]], names: Tuple[str, ...]) -> Optional[float]:
    vals: List[float] = []
    for name in names:
        vals.extend(atoms.get(name, []))
    return statistics.median(vals) if vals else None


def _fmt_shift(value: Optional[float]) -> str:
    if value is None:
        return "0.00"
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    if "." not in text:
        text += ".0"
    return text


def find_star(bmrb_id: str, cs_dir: Path) -> Path:
    for suffix in ("_21.str", "_3.str"):
        path = cs_dir / f"{bmrb_id}{suffix}"
        if path.is_file() and path.stat().st_size > 0:
            return path
    raise FileNotFoundError(f"No STAR file for BMRB {bmrb_id} under {cs_dir}")


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)[:120]


def _cache_key(bmrb_id: str, saveframe: str, sequence: str) -> str:
    digest = hashlib.sha1(f"{bmrb_id}|{saveframe}|{sequence}|inc1|sparse3".encode()).hexdigest()[:12]
    return f"{bmrb_id}__{_safe_name(saveframe)}__{digest}"


def _parse_rows_to_banks(
    rows: List[Dict[str, str]],
    *,
    is_v3: bool,
) -> Tuple[str, Dict[int, Dict[str, List[float]]]]:
    residue_to_atoms: Dict[Tuple[int, str], Dict[str, List[float]]] = {}
    for r in rows:
        seq_id: Optional[int] = None
        comp_id = ""
        atom_id = ""
        value: Optional[float] = None
        for tag, val in r.items():
            if tag == "__tags_order__":
                continue
            if is_v3:
                if "_Atom_chem_shift.Seq_ID" in tag and seq_id is None:
                    try:
                        seq_id = int(val)
                    except (TypeError, ValueError):
                        pass
                elif "_Atom_chem_shift.Comp_ID" in tag:
                    comp_id = str(val).strip()
                elif "_Atom_chem_shift.Atom_ID" in tag:
                    atom_id = str(val).strip().upper()
                elif tag == "_Atom_chem_shift.Val":
                    try:
                        value = float(val)
                    except (TypeError, ValueError):
                        pass
            else:
                if "Residue_seq_code" in tag and seq_id is None:
                    try:
                        seq_id = int(val)
                    except (TypeError, ValueError):
                        pass
                elif "Residue_label" in tag:
                    comp_id = str(val).strip()
                elif "Atom_name" in tag:
                    atom_id = str(val).strip().upper()
                elif tag == "_Chem_shift_value":
                    try:
                        value = float(val)
                    except (TypeError, ValueError):
                        pass

        if value is None and atom_id:
            order_s = r.get("__tags_order__", "")
            if order_s:
                order = order_s.split("\t")
                atom_idx = -1
                needle = "_Atom_chem_shift.Atom_ID" if is_v3 else "Atom_name"
                for j, t in enumerate(order):
                    if needle in t:
                        atom_idx = j
                        break
                if atom_idx != -1 and (atom_idx + 1) < len(order):
                    raw_val = r.get(order[atom_idx + 1])
                    if raw_val not in (None, ".", "?"):
                        try:
                            value = float(raw_val)
                        except (TypeError, ValueError):
                            pass

        if seq_id is None or not comp_id or not atom_id or value is None:
            continue
        bank = residue_to_atoms.setdefault((seq_id, comp_id), {})
        bank.setdefault(atom_id, []).append(value)

    sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
    sequence = "".join(AA3_TO1.get(comp.upper(), "X") for _, comp in sorted_keys)
    by_pos: Dict[int, Dict[str, List[float]]] = {}
    for seq_pos, key in enumerate(sorted_keys, 1):
        by_pos[seq_pos] = residue_to_atoms[key]
    return sequence, by_pos


def extract_saveframe_atom_banks(
    star_path: Path,
) -> List[Tuple[str, str, Dict[int, Dict[str, List[float]]]]]:
    """Return (saveframe_name, sequence, atoms_by_seq_pos) matching CSP numbering."""
    from scripts.bmrb_io import (
        _detect_bmrb_format,
        _parse_all_chem_shift_saveframes,
        _parse_chem_shifts_v3,
        _tokenize_star_lines,
    )

    text = star_path.read_text(encoding="utf-8", errors="ignore")
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(str(star_path))

    if fmt == "3":
        saveframes = _parse_chem_shifts_v3(lines)
        all_rows: List[Dict[str, str]] = []
        for _name, rows in saveframes:
            all_rows.extend(rows)
        sequence, by_pos = _parse_rows_to_banks(all_rows, is_v3=True)
        if not sequence:
            return []
        return [("assigned_chemical_shifts_1", sequence, by_pos)]

    results: List[Tuple[str, str, Dict[int, Dict[str, List[float]]]]] = []
    for saveframe_name, rows in _parse_all_chem_shift_saveframes(lines):
        if not rows:
            continue
        sequence, by_pos = _parse_rows_to_banks(rows, is_v3=False)
        # Mirror CSP filter: require at least some amide H and N
        has_h = any(_median(atoms, AMIDE_H) is not None for atoms in by_pos.values())
        has_n = any(_median(atoms, ("N",)) is not None for atoms in by_pos.values())
        if sequence and has_h and has_n:
            results.append((saveframe_name, sequence, by_pos))
    return results


def select_best_saveframe_pair(
    apo_banks: List[Tuple[str, str, Dict[int, Dict[str, List[float]]]]],
    holo_banks: List[Tuple[str, str, Dict[int, Dict[str, List[float]]]]],
) -> Tuple[
    Tuple[str, str, Dict[int, Dict[str, List[float]]]],
    Tuple[str, str, Dict[int, Dict[str, List[float]]]],
    float,
]:
    from scripts.align import align_global

    best_score = float("-inf")
    best: Optional[
        Tuple[
            Tuple[str, str, Dict[int, Dict[str, List[float]]]],
            Tuple[str, str, Dict[int, Dict[str, List[float]]]],
        ]
    ] = None
    for apo in apo_banks:
        for holo in holo_banks:
            _a, _h, _map, score = align_global(apo[1], holo[1])
            if score > best_score:
                best_score = score
                best = (apo, holo)
    if best is None:
        raise RuntimeError("No apo/holo saveframe pair could be aligned")
    return best[0], best[1], best_score


def write_shifty(
    out_path: Path,
    sequence: str,
    atoms_by_pos: Dict[int, Dict[str, List[float]]],
    *,
    min_atom_count: int = 3,
) -> None:
    """Write SHIFTY with CSP seq_pos numbering.

    Atom columns with fewer than ``min_atom_count`` assigned values are written
    as 0.00 so the legacy RCI smoother is not handed 1-residue atom series
    (which raises IndexError).
    """
    # Build per-column values first, then blank sparse columns.
    rows: List[Dict[str, Optional[float]]] = []
    for seq_pos, aa in enumerate(sequence, 1):
        atoms = atoms_by_pos.get(seq_pos, {})
        rows.append(
            {
                "HA": _median(atoms, ALPHA_H),
                "CA": _median(atoms, ("CA",)),
                "CB": _median(atoms, ("CB",)),
                "CO": _median(atoms, CO_NAMES),
                "N": _median(atoms, ("N",)),
                "HN": _median(atoms, AMIDE_H),
            }
        )
    keep = {
        col: sum(1 for r in rows if r[col] is not None) >= min_atom_count
        for col in ("HA", "CA", "CB", "CO", "N", "HN")
    }

    lines = ["#NUM\tAA\tHA\tCA\tCB\tCO\tN\tHN"]
    for seq_pos, aa in enumerate(sequence, 1):
        vals = rows[seq_pos - 1]
        lines.append(
            "\t".join(
                [
                    str(seq_pos),
                    aa,
                    _fmt_shift(vals["HA"] if keep["HA"] else None),
                    _fmt_shift(vals["CA"] if keep["CA"] else None),
                    _fmt_shift(vals["CB"] if keep["CB"] else None),
                    _fmt_shift(vals["CO"] if keep["CO"] else None),
                    _fmt_shift(vals["N"] if keep["N"] else None),
                    _fmt_shift(vals["HN"] if keep["HN"] else None),
                ]
            )
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_rci_txt(path: Path) -> Dict[int, Tuple[float, str]]:
    out: Dict[int, Tuple[float, str]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            resi = int(float(parts[0]))
            rci = float(parts[1])
        except ValueError:
            continue
        aa = parts[2].strip().upper()
        if len(aa) > 1:
            aa = AA3_TO1.get(aa, aa[0])
        out[resi] = (rci, aa)
    return out


def run_rci_on_shifty(shifty_path: Path, work_dir: Path, python_exe: str) -> Dict[int, Tuple[float, str]]:
    work_dir.mkdir(parents=True, exist_ok=True)
    local_input = work_dir / shifty_path.name
    if local_input.resolve() != shifty_path.resolve():
        shutil.copy2(shifty_path, local_input)

    script = _patched_rci_script()
    cmd = [python_exe, str(script), "-b", local_input.name, *RCI_CLI_FLAGS]
    last_err = ""
    for attempt in range(2):
        proc = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
        )
        rci_path = work_dir / f"{local_input.name}.RCI.txt"
        if proc.returncode == 0 and rci_path.is_file() and rci_path.stat().st_size > 0:
            return parse_rci_txt(rci_path)
        last_err = (proc.stderr or proc.stdout or "")[-800:]
        # Clear partial artifacts before retry
        for p in work_dir.glob(f"{local_input.name}*"):
            if p.name == local_input.name:
                continue
            if p.is_file():
                p.unlink()
            elif p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
    raise RuntimeError(
        f"RCI failed for {shifty_path.name} (rc={proc.returncode}): {last_err}"
    )


def get_or_compute_rci(
    bmrb_id: str,
    saveframe: str,
    sequence: str,
    atoms_by_pos: Dict[int, Dict[str, List[float]]],
    dest_dir: Path,
    python_exe: str,
    cache_dir: Path,
    role: str,
) -> Dict[int, Tuple[float, str]]:
    key = _cache_key(bmrb_id, saveframe, sequence)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached_rci = cache_dir / f"{key}.RCI.txt"
    cached_shifty = cache_dir / f"{key}.shifty.txt"

    dest_dir.mkdir(parents=True, exist_ok=True)
    out_shifty = dest_dir / f"{role}.shifty.txt"
    out_rci = dest_dir / f"{role}.RCI.txt"

    if not (cached_rci.is_file() and cached_shifty.is_file()):
        write_shifty(cached_shifty, sequence, atoms_by_pos)
        with tempfile.TemporaryDirectory(prefix=f"rci_{bmrb_id}_", dir=str(cache_dir)) as tmp:
            result = run_rci_on_shifty(cached_shifty, Path(tmp), python_exe)
            produced = Path(tmp) / f"{cached_shifty.name}.RCI.txt"
            if produced.is_file():
                shutil.copy2(produced, cached_rci)
            else:
                lines = [f"{resi} {val} {aa}" for resi, (val, aa) in sorted(result.items())]
                cached_rci.write_text("\n".join(lines) + "\n", encoding="utf-8")

    shutil.copy2(cached_shifty, out_shifty)
    shutil.copy2(cached_rci, out_rci)
    return parse_rci_txt(cached_rci)


def read_master(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    return fieldnames, rows


def write_master(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _insert_rci_fields(fieldnames: List[str]) -> List[str]:
    cleaned = [c for c in fieldnames if c not in ("rci_apo", "rci_holo")]
    if "holo_aa" in cleaned:
        idx = cleaned.index("holo_aa") + 1
        return cleaned[:idx] + ["rci_apo", "rci_holo"] + cleaned[idx:]
    return cleaned + ["rci_apo", "rci_holo"]


def _seq_aa(sequence: str, resi: int) -> Optional[str]:
    if resi < 1 or resi > len(sequence):
        return None
    return sequence[resi - 1]


def process_target(
    tgt_dir: str,
    cs_dir: str,
    python_exe: str,
    cache_dir: str,
) -> Dict[str, Any]:
    target = Path(tgt_dir)
    master_path = target / "master_alignment.csv"
    fieldnames, rows = read_master(master_path)
    if not rows:
        return {"name": target.name, "ok": False, "error": "empty master_alignment"}

    apo_bmrb = str(rows[0].get("apo_bmrb", "")).strip()
    holo_bmrb = str(rows[0].get("holo_bmrb", "")).strip()
    if not apo_bmrb or not holo_bmrb:
        return {"name": target.name, "ok": False, "error": "missing apo/holo BMRB ids"}

    cs = Path(cs_dir)
    apo_star = find_star(apo_bmrb, cs)
    holo_star = find_star(holo_bmrb, cs)
    apo_banks = extract_saveframe_atom_banks(apo_star)
    holo_banks = extract_saveframe_atom_banks(holo_star)
    if not apo_banks or not holo_banks:
        return {
            "name": target.name,
            "ok": False,
            "error": f"no usable saveframes apo={len(apo_banks)} holo={len(holo_banks)}",
        }

    apo_sel, holo_sel, score = select_best_saveframe_pair(apo_banks, holo_banks)
    apo_sf, apo_seq, apo_atoms = apo_sel
    holo_sf, holo_seq, holo_atoms = holo_sel

    # Sanity: master residues must match selected sequences at CSP indices
    mismatches = 0
    checked = 0
    for row in rows:
        try:
            apo_resi = int(float(row["apo_resi"]))
            holo_resi = int(float(row["holo_resi"]))
        except (KeyError, TypeError, ValueError):
            continue
        apo_aa = (row.get("apo_aa") or "").strip().upper()
        holo_aa = (row.get("holo_aa") or "").strip().upper()
        if apo_aa:
            checked += 1
            if _seq_aa(apo_seq, apo_resi) != apo_aa:
                mismatches += 1
        if holo_aa:
            checked += 1
            if _seq_aa(holo_seq, holo_resi) != holo_aa:
                mismatches += 1
    mismatch_frac = (mismatches / checked) if checked else 1.0
    if mismatch_frac > 0.05:
        return {
            "name": target.name,
            "ok": False,
            "error": (
                f"saveframe mismatch vs master: {mismatches}/{checked} "
                f"(apo_sf={apo_sf}, holo_sf={holo_sf}, score={score:.1f})"
            ),
        }

    rci_dir = target / "rci"
    rci_dir.mkdir(parents=True, exist_ok=True)
    (rci_dir / "meta.txt").write_text(
        f"apo_bmrb={apo_bmrb}\napo_saveframe={apo_sf}\n"
        f"holo_bmrb={holo_bmrb}\nholo_saveframe={holo_sf}\n"
        f"alignment_score={score}\n",
        encoding="utf-8",
    )

    apo_rci = get_or_compute_rci(
        apo_bmrb,
        apo_sf,
        apo_seq,
        apo_atoms,
        rci_dir,
        python_exe,
        Path(cache_dir),
        role="apo",
    )
    holo_rci = get_or_compute_rci(
        holo_bmrb,
        holo_sf,
        holo_seq,
        holo_atoms,
        rci_dir,
        python_exe,
        Path(cache_dir),
        role="holo",
    )

    filled_apo = filled_holo = aa_skip = 0
    for row in rows:
        row["rci_apo"] = ""
        row["rci_holo"] = ""
        try:
            apo_resi = int(float(row["apo_resi"]))
            holo_resi = int(float(row["holo_resi"]))
        except (KeyError, TypeError, ValueError):
            continue
        apo_aa = (row.get("apo_aa") or "").strip().upper()
        holo_aa = (row.get("holo_aa") or "").strip().upper()

        if apo_resi in apo_rci:
            val, aa = apo_rci[apo_resi]
            if apo_aa and aa and aa != apo_aa:
                aa_skip += 1
            else:
                row["rci_apo"] = f"{val:.6g}"
                filled_apo += 1
        if holo_resi in holo_rci:
            val, aa = holo_rci[holo_resi]
            if holo_aa and aa and aa != holo_aa:
                aa_skip += 1
            else:
                row["rci_holo"] = f"{val:.6g}"
                filled_holo += 1

    new_fields = _insert_rci_fields(fieldnames)
    write_master(master_path, new_fields, rows)
    return {
        "name": target.name,
        "ok": True,
        "n_rows": len(rows),
        "filled_apo": filled_apo,
        "filled_holo": filled_holo,
        "aa_skip": aa_skip,
        "apo_sf": apo_sf,
        "holo_sf": holo_sf,
        "score": score,
        "mismatch_frac": mismatch_frac,
    }


def _discover_targets(outputs: Path, ids: Optional[List[str]], limit: Optional[int]) -> List[Path]:
    targets = sorted(
        p for p in outputs.iterdir() if p.is_dir() and (p / "master_alignment.csv").is_file()
    )
    if ids:
        want = {x.strip() for x in ids if x.strip()}
        targets = [p for p in targets if p.name in want]
    if limit is not None:
        targets = targets[:limit]
    return targets


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_ROOT / "outputs",
        help="Root with per-target master_alignment.csv (default: outputs/)",
    )
    parser.add_argument(
        "--cs-dir",
        type=Path,
        default=_ROOT / "CS_Lists",
        help="BMRB STAR cache directory",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=RCI_CACHE_DIR,
        help="Shared RCI result cache",
    )
    parser.add_argument("--ids", nargs="*", help="Optional target dir names to process")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N targets")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers (default 1)")
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used to run RCI_GUI (default: current interpreter)",
    )
    args = parser.parse_args(argv)

    if not RCI_SCRIPT.is_file():
        print(f"ERROR: RCI script not found: {RCI_SCRIPT}", file=sys.stderr)
        return 2

    targets = _discover_targets(args.outputs, args.ids, args.limit)
    print(f"Found {len(targets)} targets under {args.outputs}", flush=True)
    if not targets:
        return 0

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    failures: List[Tuple[str, str]] = []
    coverage: List[Tuple[str, int, int, int]] = []

    workers = max(1, int(args.workers))
    if workers == 1:
        results_iter = (
            process_target(str(t), str(args.cs_dir), args.python, str(args.cache_dir))
            for t in targets
        )
        enumerated = enumerate(results_iter, 1)
        for i, result in enumerated:
            if result.get("ok"):
                ok += 1
                coverage.append(
                    (
                        result["name"],
                        result["n_rows"],
                        result["filled_apo"],
                        result["filled_holo"],
                    )
                )
            else:
                fail += 1
                failures.append((result["name"], str(result.get("error"))))
                print(f"FAIL {result['name']}: {result.get('error')}", flush=True)
            if i % 10 == 0 or i == len(targets):
                print(f"[{i}/{len(targets)}] ok={ok} fail={fail}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {
                ex.submit(
                    process_target,
                    str(t),
                    str(args.cs_dir),
                    args.python,
                    str(args.cache_dir),
                ): t
                for t in targets
            }
            for i, fut in enumerate(as_completed(futs), 1):
                try:
                    result = fut.result()
                except Exception as e:
                    fail += 1
                    name = futs[fut].name
                    failures.append((name, str(e)))
                    print(f"FAIL {name}: {e}", flush=True)
                    continue
                if result.get("ok"):
                    ok += 1
                    coverage.append(
                        (
                            result["name"],
                            result["n_rows"],
                            result["filled_apo"],
                            result["filled_holo"],
                        )
                    )
                else:
                    fail += 1
                    failures.append((result["name"], str(result.get("error"))))
                    print(f"FAIL {result['name']}: {result.get('error')}", flush=True)
                if i % 10 == 0 or i == len(targets):
                    print(f"[{i}/{len(targets)}] ok={ok} fail={fail}", flush=True)

    if coverage:
        fracs = [
            (a / n if n else 0.0, h / n if n else 0.0) for _name, n, a, h in coverage
        ]
        mean_apo = sum(a for a, _ in fracs) / len(fracs)
        mean_holo = sum(h for _, h in fracs) / len(fracs)
        print(
            f"Coverage mean filled_apo={mean_apo:.1%} filled_holo={mean_holo:.1%} "
            f"over {len(coverage)} targets",
            flush=True,
        )
        for name, n, a, h in coverage[:5]:
            print(f"  sample {name}: apo {a}/{n} holo {h}/{n}", flush=True)

    print(f"Done. ok={ok} fail={fail}", flush=True)
    for name, err in failures[:40]:
        print(f"  {name}: {err}")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
