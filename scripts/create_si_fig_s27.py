#!/usr/bin/env python3
"""
SI Fig. S27 — Case study panel for holo 2FIN / apo BMRB 6809.

Copies the precomputed case-study z panel into ``figures/SF27_2FIN_case_study_z.png``.

Default source: outputs/2FIN_6809/2FIN_case_study_z.png
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
DEFAULT_SRC = _REPO / "outputs" / "2FIN_6809" / "2FIN_case_study_z.png"
DEFAULT_OUT = _REPO / "figures" / "SF27_2FIN_case_study_z.png"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    src = args.source if args.source.is_absolute() else _REPO / args.source
    out = args.output if args.output.is_absolute() else _REPO / args.output
    if not src.is_file():
        print(f"Error: missing 2FIN case-study panel: {src}", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, out)
    print(f"[SF27] Wrote {out} (from {src})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
