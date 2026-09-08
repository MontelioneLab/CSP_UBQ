#!/usr/bin/env python3
"""
SI Fig. S24 — Five selected apo–apo control panels (S24A–S24E).

Copies the precomputed aggregate figures into ``figures/SF24{A–E}_*.png``.

Default sources (in order):
  figures/selected_apo_apo_controls/{query}_{match}.png
  else outputs/apo_apo_matches/{query}_{match}/aggregate_csp_hsqc_grid.png
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
DEFAULT_SELECTED = _REPO / "figures" / "selected_apo_apo_controls"
DEFAULT_PAIRS_ROOT = _REPO / "outputs" / "apo_apo_matches"
DEFAULT_OUT_DIR = _REPO / "figures"

# Curated order used by compose_apo_apo_aggregate_grid.py
PANELS: tuple[tuple[str, str], ...] = (
    ("A", "52080_52079"),
    ("B", "28070_28071"),
    ("C", "34000_34001"),
    ("D", "34394_6354"),
    ("E", "17769_51725"),
)


def _source_for_pair(pair_id: str, selected_dir: Path, pairs_root: Path) -> Path:
    selected = selected_dir / f"{pair_id}.png"
    if selected.is_file():
        return selected
    aggregate = pairs_root / pair_id / "aggregate_csp_hsqc_grid.png"
    if aggregate.is_file():
        return aggregate
    raise FileNotFoundError(
        f"Missing apo–apo panel for {pair_id}: tried {selected} and {aggregate}"
    )


def _output_name(letter: str, pair_id: str) -> str:
    return f"SF24{letter}_apo_apo_{pair_id}.png"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--selected-dir",
        type=Path,
        default=DEFAULT_SELECTED,
        help="Directory of selected pair PNGs (default: figures/selected_apo_apo_controls).",
    )
    ap.add_argument(
        "--pairs-root",
        type=Path,
        default=DEFAULT_PAIRS_ROOT,
        help="Fallback root of apo–apo pair directories.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Destination directory for SF24A–E PNGs (default: figures).",
    )
    args = ap.parse_args(argv)

    selected_dir = args.selected_dir if args.selected_dir.is_absolute() else _REPO / args.selected_dir
    pairs_root = args.pairs_root if args.pairs_root.is_absolute() else _REPO / args.pairs_root
    out_dir = args.output_dir if args.output_dir.is_absolute() else _REPO / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    failed = 0
    for letter, pair_id in PANELS:
        try:
            src = _source_for_pair(pair_id, selected_dir, pairs_root)
        except FileNotFoundError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            failed += 1
            continue
        out = out_dir / _output_name(letter, pair_id)
        shutil.copy2(src, out)
        print(f"[SF24{letter}] Wrote {out} (from {src})")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
