#!/usr/bin/env python3
"""
Shared helper for resolving --targets-csv rows to ``outputs/<dir>`` directories.

**Canonical layout (default):** each pipeline target directory basename is
``{HOLO_PDB.upper()}_{apo_bmrb}`` (for example ``1D5G_34688``). Resolution is a
direct lookup on that basename (case-insensitive index key).

**Legacy layout:** older trees used ``holo_pdb`` only or ``holo_pdb_<n>`` for
duplicate PDB codes, matched via ``master_alignment.csv`` BMRB pairs. Set
environment variable ``CSP_LEGACY_OUTPUT_DIRS=1`` (or ``true``/``yes``) to
fall back to that behaviour when the canonical dirname is missing or when a row
has no ``apo_bmrb`` (holo-only CLI rows).

Use ``resolve_target_rows`` from figure-creation scripts: it returns the
deduplicated list of resolved ``Path``s, in input row order.
"""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


_REQUIRED_COL = "holo_pdb"
_OPTIONAL_COLS = ("apo_bmrb", "holo_bmrb")


def canonical_output_dir_name(holo_pdb: str, apo_bmrb: str) -> str:
    """Return ``outputs`` subdirectory basename ``{HOLO_PDB.upper()}_{apo_bmrb}`` (stripped inputs)."""
    h = (holo_pdb or "").strip()
    a = (apo_bmrb or "").strip()
    return f"{h.upper()}_{a}"


def _legacy_output_dirs_enabled() -> bool:
    return os.environ.get("CSP_LEGACY_OUTPUT_DIRS", "").lower() in ("1", "true", "yes")


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
    optional. ``apo_bmrb`` is required for canonical dirname resolution unless
    ``CSP_LEGACY_OUTPUT_DIRS`` is set.
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
            rest = k[len(pref) :]
            if rest.isdigit():
                scored.append((int(rest), k, p))
    scored.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in scored]


def _legacy_resolve_row(
    row: TargetRow,
    outputs_index: Dict[str, Path],
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]],
) -> Optional[Path]:
    """Resolve using pre-canonical directory naming and BMRB pair matching."""
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


def resolve_row(
    row: TargetRow,
    outputs_index: Dict[str, Path],
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]],
) -> Optional[Path]:
    """Map a single target row to one outputs/<dir>; ``None`` if no matching dir exists."""
    raw_h = row.holo_pdb.strip()
    if not raw_h:
        return None
    apo = row.apo_bmrb.strip()
    legacy = _legacy_output_dirs_enabled()

    if apo:
        want_key = canonical_output_dir_name(raw_h, apo).lower()
        hit = outputs_index.get(want_key)
        if hit is not None:
            return hit
        if legacy:
            return _legacy_resolve_row(row, outputs_index, bmrb_cache)
        return None

    if legacy:
        return _legacy_resolve_row(row, outputs_index, bmrb_cache)
    return None


def csv_row_to_target_row(row: Mapping[str, str]) -> TargetRow:
    """Build a ``TargetRow`` from a CSP_UBQ-style dict (supports ``holo_pdb_id``)."""
    holo = (row.get("holo_pdb") or row.get("holo_pdb_id") or "").strip()
    return TargetRow(
        holo_pdb=holo,
        apo_bmrb=(row.get("apo_bmrb") or "").strip(),
        holo_bmrb=(row.get("holo_bmrb") or "").strip(),
    )


def logical_output_dirname_for_manifest(row: Mapping[str, str]) -> str:
    """Canonical ``outputs`` subdirectory basename for logging; empty if holo or apo missing."""
    holo = (row.get("holo_pdb") or row.get("holo_pdb_id") or "").strip()
    apo = (row.get("apo_bmrb") or "").strip()
    if not holo or not apo:
        return ""
    return canonical_output_dir_name(holo, apo)


def resolve_output_dir_from_csv_row(
    row: Mapping[str, str],
    outputs_dir: Optional[Path] = None,
    *,
    caches: Optional[Tuple[Dict[str, Path], Dict[Path, Optional[Tuple[str, str]]]]] = None,
) -> Optional[Path]:
    """Resolve one CSV row to an ``outputs/<dir>`` path (canonical lookup + optional legacy via env).

    Pass ``outputs_dir`` alone, or ``caches`` from :func:`build_resolution_caches` for batch use.
    """
    tr = csv_row_to_target_row(row)
    if caches is None:
        if outputs_dir is None:
            raise ValueError("resolve_output_dir_from_csv_row requires outputs_dir when caches is omitted")
        outputs_index, bmrb_cache = build_resolution_caches(outputs_dir)
    else:
        outputs_index, bmrb_cache = caches
    return resolve_row(tr, outputs_index, bmrb_cache)


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
