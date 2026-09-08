#!/usr/bin/env python3
"""
Filter the same-author apo/holo list to pairs with matched ionic strength and
spectrometer field strength (plus pH / temperature), using apo_holo_exp_conditions.csv.

Default criteria (both sides must report each quantity):
  |ΔpH| ≤ 0.5, |ΔT| ≤ 5 °C, |ΔI| ≤ 50 mM, exact primary field-strength MHz match.

Writes:
  - refined target CSV (same schema as input)
  - audit CSV with per-pair deltas and pass/fail flags
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.config import Paths  # noqa: E402


def _f(val: Optional[str]) -> Optional[float]:
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _bool_flag(ok: bool) -> str:
    return "True" if ok else "False"


def evaluate_pair(
    exp: dict[str, str],
    *,
    ph_max: float,
    temp_max_c: float,
    ionic_max_mm: float,
) -> dict[str, object]:
    aph = _f(exp.get("apo_pH"))
    hph = _f(exp.get("holo_pH"))
    at = _f(exp.get("apo_temperature_C"))
    ht = _f(exp.get("holo_temperature_C"))
    ai = _f(exp.get("apo_ionic_strength_mM"))
    hi = _f(exp.get("holo_ionic_strength_mM"))
    af = _f(exp.get("apo_field_strength_MHz"))
    hf = _f(exp.get("holo_field_strength_MHz"))

    d_ph = abs(aph - hph) if aph is not None and hph is not None else None
    d_t = abs(at - ht) if at is not None and ht is not None else None
    d_i = abs(ai - hi) if ai is not None and hi is not None else None

    ph_ok = d_ph is not None and d_ph <= ph_max
    temp_ok = d_t is not None and d_t <= temp_max_c
    ionic_ok = d_i is not None and d_i <= ionic_max_mm
    spec_ok = af is not None and hf is not None and af == hf
    keep = ph_ok and temp_ok and ionic_ok and spec_ok

    return {
        "apo_pH": aph,
        "holo_pH": hph,
        "delta_pH": d_ph,
        "apo_temperature_C": at,
        "holo_temperature_C": ht,
        "delta_T_C": d_t,
        "apo_ionic_strength_mM": ai,
        "holo_ionic_strength_mM": hi,
        "delta_ionic_mM": d_i,
        "apo_field_strength_MHz": af,
        "holo_field_strength_MHz": hf,
        "ph_ok": ph_ok,
        "temp_ok": temp_ok,
        "ionic_ok": ionic_ok,
        "spectrometer_ok": spec_ok,
        "keep": keep,
        "missing_ionic": ai is None or hi is None,
        "missing_field": af is None or hf is None,
        "mismatched_field": (
            af is not None and hf is not None and af != hf
        ),
    }


def _fmt_num(x: object, nd: int = 4) -> str:
    if x is None:
        return ""
    try:
        v = float(x)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return ""
    r = round(v, nd)
    if nd <= 1 and r == int(r):
        return str(int(r))
    return str(r)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    cfg = Paths()
    ap.add_argument(
        "--input",
        type=Path,
        default=Path(cfg.data_dir) / "CSP_UBQ_ph0.5_temp5C_same_author_list.csv",
        help="Same-author target list CSV",
    )
    ap.add_argument(
        "--exp",
        type=Path,
        default=Path(cfg.exp_conditions_csv),
        help="apo_holo_exp_conditions CSV",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(cfg.data_dir)
        / "CSP_UBQ_ph0.5_temp5C_same_author_matched_conditions.csv",
        help="Refined target CSV",
    )
    ap.add_argument(
        "--audit",
        type=Path,
        default=Path(cfg.data_dir)
        / "CSP_UBQ_ph0.5_temp5C_same_author_condition_audit.csv",
        help="Per-pair diagnostic CSV",
    )
    ap.add_argument("--ph-max-diff", type=float, default=0.5)
    ap.add_argument("--temp-max-diff-c", type=float, default=5.0)
    ap.add_argument("--ionic-max-diff-mm", type=float, default=50.0)
    args = ap.parse_args()

    in_path = args.input if args.input.is_absolute() else _REPO / args.input
    exp_path = args.exp if args.exp.is_absolute() else _REPO / args.exp
    out_path = args.output if args.output.is_absolute() else _REPO / args.output
    audit_path = args.audit if args.audit.is_absolute() else _REPO / args.audit

    if not in_path.is_file() or not exp_path.is_file():
        print(f"Missing input: {in_path} or {exp_path}", file=sys.stderr)
        return 1

    with in_path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        fieldnames = r.fieldnames
        pairs = list(r)
    with exp_path.open(newline="", encoding="utf-8") as f:
        exp_by_bmrb = {
            (row["apo_bmrb"].strip(), row["holo_bmrb"].strip()): row
            for row in csv.DictReader(f)
        }

    if fieldnames is None:
        print("Input CSV has no header.", file=sys.stderr)
        return 1

    kept: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    n_missing_exp = 0
    n_missing_ionic = 0
    n_mismatched_field = 0
    n_missing_field = 0

    for pair in pairs:
        key = (pair["apo_bmrb"].strip(), pair["holo_bmrb"].strip())
        exp = exp_by_bmrb.get(key)
        if exp is None:
            n_missing_exp += 1
            audit_rows.append(
                {
                    "apo_bmrb": pair["apo_bmrb"],
                    "holo_bmrb": pair["holo_bmrb"],
                    "apo_pdb": pair.get("apo_pdb", ""),
                    "holo_pdb": pair.get("holo_pdb", ""),
                    "keep": "False",
                    "ph_ok": "False",
                    "temp_ok": "False",
                    "ionic_ok": "False",
                    "spectrometer_ok": "False",
                    "missing_ionic": "True",
                    "missing_field": "True",
                    "mismatched_field": "False",
                    "note": "no_exp_conditions_row",
                }
            )
            continue

        ev = evaluate_pair(
            exp,
            ph_max=args.ph_max_diff,
            temp_max_c=args.temp_max_diff_c,
            ionic_max_mm=args.ionic_max_diff_mm,
        )
        if ev["missing_ionic"]:
            n_missing_ionic += 1
        if ev["missing_field"]:
            n_missing_field += 1
        if ev["mismatched_field"]:
            n_mismatched_field += 1
        if ev["keep"]:
            kept.append(pair)

        audit_rows.append(
            {
                "apo_bmrb": pair["apo_bmrb"],
                "holo_bmrb": pair["holo_bmrb"],
                "apo_pdb": pair.get("apo_pdb", ""),
                "holo_pdb": pair.get("holo_pdb", ""),
                "apo_pH": _fmt_num(ev["apo_pH"], 3),
                "holo_pH": _fmt_num(ev["holo_pH"], 3),
                "delta_pH": _fmt_num(ev["delta_pH"], 3),
                "apo_temperature_C": _fmt_num(ev["apo_temperature_C"], 2),
                "holo_temperature_C": _fmt_num(ev["holo_temperature_C"], 2),
                "delta_T_C": _fmt_num(ev["delta_T_C"], 2),
                "apo_ionic_strength_mM": _fmt_num(ev["apo_ionic_strength_mM"], 4),
                "holo_ionic_strength_mM": _fmt_num(ev["holo_ionic_strength_mM"], 4),
                "delta_ionic_mM": _fmt_num(ev["delta_ionic_mM"], 4),
                "apo_field_strength_MHz": _fmt_num(ev["apo_field_strength_MHz"], 1),
                "holo_field_strength_MHz": _fmt_num(ev["holo_field_strength_MHz"], 1),
                "ph_ok": _bool_flag(bool(ev["ph_ok"])),
                "temp_ok": _bool_flag(bool(ev["temp_ok"])),
                "ionic_ok": _bool_flag(bool(ev["ionic_ok"])),
                "spectrometer_ok": _bool_flag(bool(ev["spectrometer_ok"])),
                "keep": _bool_flag(bool(ev["keep"])),
                "missing_ionic": _bool_flag(bool(ev["missing_ionic"])),
                "missing_field": _bool_flag(bool(ev["missing_field"])),
                "mismatched_field": _bool_flag(bool(ev["mismatched_field"])),
                "note": "",
            }
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)

    audit_fields = [
        "apo_bmrb",
        "holo_bmrb",
        "apo_pdb",
        "holo_pdb",
        "apo_pH",
        "holo_pH",
        "delta_pH",
        "apo_temperature_C",
        "holo_temperature_C",
        "delta_T_C",
        "apo_ionic_strength_mM",
        "holo_ionic_strength_mM",
        "delta_ionic_mM",
        "apo_field_strength_MHz",
        "holo_field_strength_MHz",
        "ph_ok",
        "temp_ok",
        "ionic_ok",
        "spectrometer_ok",
        "keep",
        "missing_ionic",
        "missing_field",
        "mismatched_field",
        "note",
    ]
    with audit_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=audit_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(audit_rows)

    print(
        f"Kept {len(kept)}/{len(pairs)} pairs → {out_path}\n"
        f"  missing exp row: {n_missing_exp}\n"
        f"  missing ionic (either side): {n_missing_ionic}\n"
        f"  missing field (either side): {n_missing_field}\n"
        f"  mismatched field: {n_mismatched_field}\n"
        f"Audit → {audit_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
