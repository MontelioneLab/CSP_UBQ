#!/usr/bin/env python3
"""
Shared helper for resolving --targets-csv rows to outputs/<holo_pdb> directories.

A "target row" carries ``holo_pdb`` plus optionally ``apo_bmrb`` and
``holo_bmrb``. The resolver maps each row to a single ``outputs/<dir>`` by:

1. Finding all candidate dirs whose lowercased basename equals ``holo_pdb`` or
   ``holo_pdb_<n>`` (numeric suffix written by the pipeline for duplicate
   ``holo_pdb`` rows).
2. Filtering candidates to the ones whose ``master_alignment.csv`` first data
   row's ``apo_bmrb``/``holo_bmrb`` pair equals the target row's pair.
3. If more than one candidate matches, returning the first (sorted by suffix
   number, then by directory name).
4. If no candidate matches, returning ``None`` and logging a single ``[WARN]``.

When the target row carries no BMRB info (e.g. a comma-separated ``--targets``
flag) the resolver falls back to "first candidate by suffix order" so legacy
CLI behavior is preserved.

Use ``resolve_target_rows`` from figure-creation scripts: it returns the
deduplicated list of resolved ``Path``s, in input row order.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


_REQUIRED_COL = "holo_pdb"
_OPTIONAL_COLS = ("apo_bmrb", "holo_bmrb")


@dataclass(frozen=True)
class TargetRow:
    """A single targets-CSV row with the columns the resolver cares about."""

    holo_pdb: str
    apo_bmrb: str = ""
    holo_bmrb: str = ""

    @property
    def has_bmrb(self) -> bool:
        return bool(self.apo_bmrb or self.holo_bmrb)


def load_target_rows(
    csv_path: Optional[Path],
    *,
    extra_holo_pdbs: Iterable[str] = (),
) -> List[TargetRow]:
    """Read ``csv_path`` and any extra ``--targets`` holo_pdb strings into ``TargetRow``s.

    The CSV must contain a ``holo_pdb`` column; ``apo_bmrb`` / ``holo_bmrb`` are
    optional but recommended (the resolver requires them to disambiguate
    duplicate-suffix dirs). Extra holo_pdb strings are appended as BMRB-less
    rows for backward compatibility with the comma-separated ``--targets`` flag.
    """
    rows: List[TargetRow] = []
    if csv_path is not None:
        with csv_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            field_lookup = {fn.lower(): fn for fn in (reader.fieldnames or [])}
            if _REQUIRED_COL not in field_lookup:
                raise ValueError(f"Missing 'holo_pdb' column in {csv_path}")
            holo_col = field_lookup[_REQUIRED_COL]
            apo_col = field_lookup.get("apo_bmrb")
            holo_b_col = field_lookup.get("holo_bmrb")
            for raw in reader:
                holo = (raw.get(holo_col) or "").strip()
                if not holo:
                    continue
                rows.append(
                    TargetRow(
                        holo_pdb=holo,
                        apo_bmrb=(raw.get(apo_col) or "").strip() if apo_col else "",
                        holo_bmrb=(raw.get(holo_b_col) or "").strip() if holo_b_col else "",
                    )
                )
    for h in extra_holo_pdbs:
        h = (h or "").strip()
        if h:
            rows.append(TargetRow(holo_pdb=h))
    return rows


def _read_first_bmrb_pair(alignment_path: Path) -> Optional[Tuple[str, str]]:
    """Return the first non-empty (apo_bmrb, holo_bmrb) pair in a master_alignment.csv."""
    if not alignment_path.is_file():
        return None
    try:
        with alignment_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                a = (row.get("apo_bmrb") or "").strip()
                h = (row.get("holo_bmrb") or "").strip()
                if a or h:
                    return a, h
    except Exception:
        return None
    return None


def build_resolution_caches(
    outputs_dir: Path,
) -> Tuple[Dict[str, Path], Dict[Path, Optional[Tuple[str, str]]]]:
    """Index ``outputs/`` subdirs and cache each ``master_alignment.csv`` first BMRB pair.

    Returns ``(outputs_index, bmrb_cache)`` where ``outputs_index`` maps the
    lowercased subdirectory basename to its ``Path`` and ``bmrb_cache`` maps
    that ``Path`` to the cached BMRB pair (or ``None`` if unavailable).
    """
    outputs_index: Dict[str, Path] = {}
    if outputs_dir.is_dir():
        for p in sorted(outputs_dir.iterdir()):
            if p.is_dir() and not p.name.startswith("."):
                outputs_index.setdefault(p.name.lower(), p)
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]] = {
        p: _read_first_bmrb_pair(p / "master_alignment.csv") for p in outputs_index.values()
    }
    return outputs_index, bmrb_cache


def _candidate_dirs_for_holo(
    outputs_index: Dict[str, Path], holo_pdb_lower: str
) -> List[Path]:
    """Dirs whose name is ``holo_pdb`` or ``holo_pdb_<n>`` (sorted by suffix number)."""
    base = holo_pdb_lower
    pref = base + "_"
    scored: List[Tuple[int, str, Path]] = []
    for k, p in outputs_index.items():
        if k == base:
            scored.append((0, k, p))
            continue
        if k.startswith(pref):
            rest = k[len(pref):]
            if rest.isdigit():
                scored.append((int(rest), k, p))
    scored.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in scored]


def resolve_row(
    row: TargetRow,
    outputs_index: Dict[str, Path],
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]],
) -> Optional[Path]:
    """Map a single target row to one outputs/<dir>; ``None`` if no BMRB-congruent dir exists."""
    raw_h = row.holo_pdb.strip().lower()
    if not raw_h:
        return None
    candidates = _candidate_dirs_for_holo(outputs_index, raw_h)
    if not candidates:
        return None
    if not row.has_bmrb:
        return candidates[0]
    want = (row.apo_bmrb, row.holo_bmrb)
    matches = [p for p in candidates if bmrb_cache.get(p) == want]
    if matches:
        return matches[0]
    return None


def resolve_target_rows(
    rows: Sequence[TargetRow],
    outputs_dir: Path,
    *,
    log_warnings: bool = True,
    log_stream=sys.stderr,
) -> List[Path]:
    """Resolve every row to a unique ``outputs/<dir>`` Path (input order, deduped)."""
    outputs_index, bmrb_cache = build_resolution_caches(outputs_dir)
    seen: Set[str] = set()
    resolved: List[Path] = []
    for row in rows:
        path = resolve_row(row, outputs_index, bmrb_cache)
        if path is None:
            if log_warnings:
                print(
                    "[WARN] No outputs subdirectory matches "
                    f"holo_pdb={row.holo_pdb!r} apo_bmrb={row.apo_bmrb!r} "
                    f"holo_bmrb={row.holo_bmrb!r}",
                    file=log_stream,
                )
            continue
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        resolved.append(path)
    return resolved


def resolved_dir_names(
    rows: Sequence[TargetRow],
    outputs_dir: Path,
    *,
    log_warnings: bool = True,
) -> Set[str]:
    """Convenience wrapper returning the set of resolved directory basenames."""
    return {p.name for p in resolve_target_rows(rows, outputs_dir, log_warnings=log_warnings)}


def load_and_resolve(
    csv_path: Optional[Path],
    outputs_dir: Path,
    *,
    extra_holo_pdbs: Iterable[str] = (),
    log_warnings: bool = True,
) -> List[Path]:
    """Convenience: load CSV (+ extra holo_pdbs) and return resolved dir Paths."""
    rows = load_target_rows(csv_path, extra_holo_pdbs=extra_holo_pdbs)
    return resolve_target_rows(rows, outputs_dir, log_warnings=log_warnings)
