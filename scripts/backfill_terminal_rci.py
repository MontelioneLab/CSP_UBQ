#!/usr/bin/env python3
"""Backfill terminal high-RCI boolean columns into master_alignment.csv.

Adds:
  - apo_high_terminal_RCI
  - holo_high_terminal_RCI
  - both_high_terminal_RCI

Terminal stretches are found by walking from the N- and C-termini of each
state's RCI series. Residues with RCI > threshold are high; up to ``max_gap``
consecutive below-threshold residues may be included; two consecutive fails
end the stretch (failing residues excluded).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

BOOL_COLS = (
    "apo_high_terminal_RCI",
    "holo_high_terminal_RCI",
    "both_high_terminal_RCI",
)


def parse_rci_txt(path: Path) -> Dict[int, float]:
    out: Dict[int, float] = {}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            resi = int(float(parts[0]))
            val = float(parts[1])
        except ValueError:
            continue
        out[resi] = val
    return out


def rci_from_master_rows(
    rows: Sequence[Dict[str, str]],
    resi_col: str,
    rci_col: str,
) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for row in rows:
        raw_resi = (row.get(resi_col) or "").strip()
        raw_rci = (row.get(rci_col) or "").strip()
        if not raw_resi or not raw_rci:
            continue
        try:
            resi = int(float(raw_resi))
            out[resi] = float(raw_rci)
        except ValueError:
            continue
    return out


def load_state_rci(tgt_dir: Path, role: str, rows: Sequence[Dict[str, str]]) -> Dict[int, float]:
    path = tgt_dir / "rci" / f"{role}.RCI.txt"
    series = parse_rci_txt(path)
    if series:
        return series
    resi_col = "apo_resi" if role == "apo" else "holo_resi"
    rci_col = "rci_apo" if role == "apo" else "rci_holo"
    return rci_from_master_rows(rows, resi_col, rci_col)


def _is_high(rci_by_resi: Dict[int, float], resi: int, threshold: float) -> bool:
    val = rci_by_resi.get(resi)
    return val is not None and val > threshold


def walk_terminal_stretch(
    ordered_resis: Sequence[int],
    rci_by_resi: Dict[int, float],
    *,
    threshold: float,
    max_gap: int,
) -> Set[int]:
    """Walk ``ordered_resis`` from the terminus; return included residue numbers."""
    included: List[int] = []
    gap_run = 0
    pending_gaps: List[int] = []

    for resi in ordered_resis:
        if _is_high(rci_by_resi, resi, threshold):
            if pending_gaps:
                included.extend(pending_gaps)
                pending_gaps = []
            included.append(resi)
            gap_run = 0
            continue

        # Below threshold (or missing): count toward the allowed gap.
        gap_run += 1
        if gap_run > max_gap:
            break
        pending_gaps.append(resi)

    return set(included)


def terminal_high_residues(
    rci_by_resi: Dict[int, float],
    *,
    threshold: float = 0.3,
    max_gap: int = 1,
) -> Set[int]:
    if not rci_by_resi:
        return set()
    # Contiguous integer span covering all observed residue numbers so gaps in
    # numbering are treated as missing (not high).
    lo = min(rci_by_resi)
    hi = max(rci_by_resi)
    ascending = list(range(lo, hi + 1))
    descending = list(range(hi, lo - 1, -1))
    n_set = walk_terminal_stretch(
        ascending, rci_by_resi, threshold=threshold, max_gap=max_gap
    )
    c_set = walk_terminal_stretch(
        descending, rci_by_resi, threshold=threshold, max_gap=max_gap
    )
    return n_set | c_set


def _insert_bool_fields(fieldnames: List[str]) -> List[str]:
    cleaned = [c for c in fieldnames if c not in BOOL_COLS]
    if "rci_holo" in cleaned:
        idx = cleaned.index("rci_holo") + 1
        return cleaned[:idx] + list(BOOL_COLS) + cleaned[idx:]
    if "holo_aa" in cleaned:
        idx = cleaned.index("holo_aa") + 1
        return cleaned[:idx] + list(BOOL_COLS) + cleaned[idx:]
    return cleaned + list(BOOL_COLS)


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


def process_target(
    tgt_dir: Path,
    *,
    threshold: float,
    max_gap: int,
) -> Dict[str, object]:
    master_path = tgt_dir / "master_alignment.csv"
    fieldnames, rows = read_master(master_path)
    if not rows:
        return {"name": tgt_dir.name, "ok": False, "error": "empty master_alignment"}

    apo_rci = load_state_rci(tgt_dir, "apo", rows)
    holo_rci = load_state_rci(tgt_dir, "holo", rows)
    apo_flags = terminal_high_residues(apo_rci, threshold=threshold, max_gap=max_gap)
    holo_flags = terminal_high_residues(holo_rci, threshold=threshold, max_gap=max_gap)

    n_apo = n_holo = n_both = 0
    for row in rows:
        apo_true = False
        holo_true = False
        try:
            apo_resi = int(float(row["apo_resi"]))
            apo_true = apo_resi in apo_flags
        except (KeyError, TypeError, ValueError):
            pass
        try:
            holo_resi = int(float(row["holo_resi"]))
            holo_true = holo_resi in holo_flags
        except (KeyError, TypeError, ValueError):
            pass
        both_true = apo_true and holo_true
        row["apo_high_terminal_RCI"] = "True" if apo_true else "False"
        row["holo_high_terminal_RCI"] = "True" if holo_true else "False"
        row["both_high_terminal_RCI"] = "True" if both_true else "False"
        n_apo += int(apo_true)
        n_holo += int(holo_true)
        n_both += int(both_true)

    new_fields = _insert_bool_fields(fieldnames)
    write_master(master_path, new_fields, rows)
    return {
        "name": tgt_dir.name,
        "ok": True,
        "n_rows": len(rows),
        "n_apo": n_apo,
        "n_holo": n_holo,
        "n_both": n_both,
        "apo_terminal_n": len(apo_flags),
        "holo_terminal_n": len(holo_flags),
    }


def _discover_targets(
    outputs: Path, ids: Optional[List[str]], limit: Optional[int]
) -> List[Path]:
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
        help="Root with per-target master_alignment.csv",
    )
    parser.add_argument("--ids", nargs="*", help="Optional target dir names")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--threshold", type=float, default=0.3)
    parser.add_argument("--max-gap", type=int, default=1)
    args = parser.parse_args(argv)

    targets = _discover_targets(args.outputs, args.ids, args.limit)
    print(
        f"Found {len(targets)} targets under {args.outputs} "
        f"(threshold={args.threshold}, max_gap={args.max_gap})",
        flush=True,
    )
    if not targets:
        return 0

    ok = fail = 0
    tot_apo = tot_holo = tot_both = tot_rows = 0
    failures: List[Tuple[str, str]] = []

    for i, tgt in enumerate(targets, 1):
        try:
            result = process_target(
                tgt, threshold=args.threshold, max_gap=args.max_gap
            )
        except Exception as e:
            fail += 1
            failures.append((tgt.name, str(e)))
            print(f"FAIL {tgt.name}: {e}", flush=True)
            continue
        if not result.get("ok"):
            fail += 1
            failures.append((str(result["name"]), str(result.get("error"))))
            print(f"FAIL {result['name']}: {result.get('error')}", flush=True)
        else:
            ok += 1
            tot_rows += int(result["n_rows"])
            tot_apo += int(result["n_apo"])
            tot_holo += int(result["n_holo"])
            tot_both += int(result["n_both"])
        if i % 25 == 0 or i == len(targets):
            print(f"[{i}/{len(targets)}] ok={ok} fail={fail}", flush=True)

    print(
        f"Residue flags: apo={tot_apo} holo={tot_holo} both={tot_both} "
        f"over {tot_rows} master rows",
        flush=True,
    )
    print(f"Done. ok={ok} fail={fail}", flush=True)
    for name, err in failures[:40]:
        print(f"  {name}: {err}")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
