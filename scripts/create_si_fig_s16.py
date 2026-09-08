#!/usr/bin/env python3
"""
SI Fig. S16 — PDB Advanced Search interface screenshot.

Static asset. Installs/verifies ``figures/SF16_pdb_search.png`` (preferred) or
leaves ``figures/SF16_pdb_search.pdf`` for ``build_si_merged_pdf.py`` to splice.

Usage:
  python scripts/create_si_fig_s16.py
  python scripts/create_si_fig_s16.py --source path/to/screenshot.png
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
DEFAULT_PNG = _REPO / "figures" / "SF16_pdb_search.png"
DEFAULT_PDF = _REPO / "figures" / "SF16_pdb_search.pdf"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--figures-dir",
        type=Path,
        default=_REPO / "figures",
        help="Figures directory (default: figures/).",
    )
    ap.add_argument(
        "--source",
        type=Path,
        default=None,
        help="Optional PNG/PDF to copy into figures/ as SF16_pdb_search.*",
    )
    ap.add_argument(
        "--require-png",
        action="store_true",
        help="Fail if SF16_pdb_search.png is missing (default: PNG or PDF is OK).",
    )
    args = ap.parse_args(argv)

    figures = args.figures_dir if args.figures_dir.is_absolute() else _REPO / args.figures_dir
    figures.mkdir(parents=True, exist_ok=True)
    png = figures / "SF16_pdb_search.png"
    pdf = figures / "SF16_pdb_search.pdf"

    if args.source is not None:
        src = args.source if args.source.is_absolute() else _REPO / args.source
        if not src.is_file():
            print(f"Error: --source not found: {src}", file=sys.stderr)
            return 1
        dest = png if src.suffix.lower() == ".png" else pdf if src.suffix.lower() == ".pdf" else png
        if src.suffix.lower() not in {".png", ".pdf"}:
            dest = png
        shutil.copy2(src, dest)
        print(f"[SF15] Installed {dest}")

    if png.is_file():
        print(f"[SF15] OK: {png}")
        return 0
    if pdf.is_file() and not args.require_png:
        print(
            f"[SF15] PNG missing; PDF present at {pdf}. "
            "build_si_merged_pdf.py will insert this page.",
            file=sys.stderr,
        )
        return 0

    print(
        "Error: SI Fig. S16 asset missing. Provide figures/SF16_pdb_search.png "
        "(or .pdf) or pass --source.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
