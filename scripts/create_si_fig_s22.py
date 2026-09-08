#!/usr/bin/env python3
"""SI Fig. S22 — FP/(TP+FP) distributions; thin wrapper for create_si_fig_s22_fp_percent.py."""

from __future__ import annotations

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.create_si_fig_s22_fp_percent import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
