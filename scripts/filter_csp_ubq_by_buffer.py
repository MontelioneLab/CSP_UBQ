#!/usr/bin/env python3
"""
Write a subset of CSP_UBQ.csv whose apo/holo rows satisfy buffer criteria from
apo_holo_exp_conditions.csv (paired by row index with CSP_UBQ).

Default: |ΔpH| ≤ 0.5 and |ΔT| ≤ 5 °C, requiring pH and temperature on both sides.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.config import Paths  # noqa: E402


def row_meets(
    exp_row: dict[str, str], *, ph_max: float, temp_max_c: float
) -> bool:
    try:
        aph = float(exp_row["apo_pH"]) if exp_row.get("apo_pH", "").strip() else None
        hph = float(exp_row["holo_pH"]) if exp_row.get("holo_pH", "").strip() else None
        at = float(exp_row["apo_temperature_C"]) if exp_row.get("apo_temperature_C", "").strip() else None
        ht = float(exp_row["holo_temperature_C"]) if exp_row.get("holo_temperature_C", "").strip() else None
    except ValueError:
        return False
    if aph is None or hph is None or at is None or ht is None:
        return False
    return abs(aph - hph) <= ph_max and abs(at - ht) <= temp_max_c


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csp", type=Path, default=None, help="Input CSP_UBQ CSV")
    ap.add_argument(
        "--exp",
        type=Path,
        default=None,
        help="apo_holo_exp_conditions.csv",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: CSP_UBQ_ph{ph}_temp{C}.csv in repo root)",
    )
    ap.add_argument("--ph-max-diff", type=float, default=0.5)
    ap.add_argument("--temp-max-diff-c", type=float, default=5.0)
    args = ap.parse_args()

    cfg = Paths()
    csp_path = args.csp or Path(cfg.input_csv)
    exp_path = args.exp or (_REPO / "apo_holo_exp_conditions.csv")
    ph = float(args.ph_max_diff)
    tc = float(args.temp_max_diff_c)
    out = args.output
    if out is None:
        # Filename encodes tolerances (avoid ugly floats)
        out = _REPO / f"CSP_UBQ_ph{ph:g}_temp{tc:g}C.csv"

    if not csp_path.is_file() or not exp_path.is_file():
        print("Missing input CSV.", file=sys.stderr)
        return 1

    with csp_path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        fieldnames = r.fieldnames
        csp_rows = list(r)
    with exp_path.open(newline="", encoding="utf-8") as f:
        exp_rows = list(csv.DictReader(f))

    if fieldnames is None:
        print("CSP CSV has no header.", file=sys.stderr)
        return 1

    n = min(len(csp_rows), len(exp_rows))
    if len(csp_rows) != len(exp_rows):
        print(
            f"WARN: row counts differ (csp={len(csp_rows)}, exp={len(exp_rows)}); using first {n} pairs.",
            file=sys.stderr,
        )

    kept = [
        csp_rows[i]
        for i in range(n)
        if row_meets(exp_rows[i], ph_max=ph, temp_max_c=tc)
    ]

    out = out if out.is_absolute() else _REPO / out
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)

    print(f"Wrote {len(kept)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
