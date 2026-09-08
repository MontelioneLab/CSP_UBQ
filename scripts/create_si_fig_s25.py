#!/usr/bin/env python3
"""
SI Fig. S25 — Terminal-anchor vs global offsets (HSQC overlays + CSP bars).

By default copies the precomputed four-target panel into
``figures/SF25_terminal_anchor_vs_global.png``.

Default source:
  outputs/hsqc_overlay_anchor_vs_global_all.png

Pass ``--regenerate`` to rebuild via ``scripts/compare_terminal_anchor_offsets.py``
(requires matching per-target outputs under ``--outputs``).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
DEFAULT_SRC = (
    _REPO / "outputs" / "hsqc_overlay_anchor_vs_global_all.png"
)
DEFAULT_OUT = _REPO / "figures" / "SF25_terminal_anchor_vs_global.png"
DEFAULT_OUTPUTS = _REPO / "outputs"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--regenerate",
        action="store_true",
        help="Run compare_terminal_anchor_offsets.py before copying.",
    )
    ap.add_argument(
        "--outputs",
        type=Path,
        default=DEFAULT_OUTPUTS,
        help="Per-target root for --regenerate (default: outputs).",
    )
    args = ap.parse_args(argv)

    src = args.source if args.source.is_absolute() else _REPO / args.source
    out = args.output if args.output.is_absolute() else _REPO / args.output
    outputs = args.outputs if args.outputs.is_absolute() else _REPO / args.outputs

    if args.regenerate:
        cmd = [
            sys.executable,
            str(_REPO / "scripts" / "compare_terminal_anchor_offsets.py"),
            "--outputs",
            str(outputs),
        ]
        print(f"[SF24] Regenerating via: {' '.join(cmd)}")
        r = subprocess.run(cmd, cwd=str(_REPO))
        if r.returncode != 0:
            return r.returncode
        src = outputs / "hsqc_overlay_anchor_vs_global_all.png"

    if not src.is_file():
        print(f"Error: missing terminal-anchor panel: {src}", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, out)
    print(f"[SF24] Wrote {out} (from {src})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
