#!/usr/bin/env python3
"""
Build apo/holo experimental-condition summary from BMRB NMR-STAR files.

Reads CSP_UBQ.csv (apo_bmrb, holo_bmrb, apo_pdb, holo_pdb), loads cached STAR
files from CS_Lists/{id}_21.str or {id}_3.str, extracts pH, temperature (°C),
ionic strength (mM), and explicit NaCl (mM, NMR-STAR 3 sample components only).

Output: apo_holo_exp_conditions.csv with conditions_similar True only when
both sides report NaCl and |ΔNaCl| <= NACL_TOLERANCE_MM, |ΔpH| <= 0.5, and
|ΔT_C| <= 5. Missing pH, temperature, or NaCl on either side yields False.

Many legacy (_21) entries lack NaCl component rows; conditions_similar is then
usually False even when ionic strength is present.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

# --- Tunable comparison thresholds ---
NACL_TOLERANCE_MM = 50.0
PH_TOLERANCE = 0.5
TEMP_TOLERANCE_C = 5.0

# Repo root (parent of scripts/)
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.config import Paths  # noqa: E402


# STAR helpers (aligned with scripts/bmrb_io.py; local to avoid importing requests).
def _tokenize_star_lines(text: str) -> List[str]:
    return text.splitlines()


def _split_star_line(line: str) -> List[str]:
    tokens: List[str] = []
    buf: List[str] = []
    in_single = False
    in_double = False
    i = 0
    while i < len(line):
        c = line[i]
        if c == "'" and not in_double:
            in_single = not in_single
            i += 1
            continue
        if c == '"' and not in_single:
            in_double = not in_double
            i += 1
            continue
        if not in_single and not in_double and c.isspace():
            if buf:
                tokens.append("".join(buf))
                buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    if buf:
        tokens.append("".join(buf))
    return tokens


def _detect_bmrb_format(star_path: str) -> str:
    filename = os.path.basename(star_path)
    if filename.endswith("_21.str"):
        return "21"
    if filename.endswith("_3.str"):
        return "3"
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    if "_Atom_chem_shift.Val" in content:
        return "3"
    if "_Chem_shift_value" in content:
        return "21"
    return "21"


@dataclass
class ExtractedConditions:
    pH: Optional[float] = None
    temperature_C: Optional[float] = None
    ionic_strength_mM: Optional[float] = None
    nacl_mM: Optional[float] = None


def _parse_float_token(tok: str) -> Optional[float]:
    t = tok.strip().strip("'\"")
    if t in (".", "", "?"):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def _normalize_condition_key(type_str: str) -> Optional[str]:
    s = type_str.strip().strip("'\"").lower()
    if s in ("ph", "p_d", "pd"):
        return "ph"
    if "temperature" in s:
        return "temperature"
    if "ionic" in s and "strength" in s:
        return "ionic_strength"
    return None


def _temp_to_celsius(val: float, units: str) -> Optional[float]:
    u = units.strip().strip("'\"").lower()
    if u in ("k", "kelvin"):
        return val - 273.15
    if u in ("c", "celsius", "°c", "degc"):
        return val
    if u in ("f", "fahrenheit"):
        return (val - 32.0) * 5.0 / 9.0
    return None


def _ionic_to_mM(val: float, units: str) -> Optional[float]:
    u = units.strip().strip("'\"").lower()
    if u in ("mm", "millimolar"):
        return val
    if u in ("m", "molar"):
        return val * 1000.0
    return None


def _iter_saveframes(lines: List[str]) -> Iterator[Tuple[int, int, str]]:
    """Yield (start_line_index, end_line_index_of_save_marker, saveframe_name)."""
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        s = raw.strip()
        if s.startswith("save_") and s != "save_":
            name = s
            start = i
            i += 1
            while i < n and lines[i].strip() != "save_":
                i += 1
            yield (start, i, name)
            i += 1
        else:
            i += 1


def _parse_loop_at(
    lines: List[str], i: int, end_bound: int
) -> Optional[Tuple[List[str], List[List[str]], int]]:
    """
    If lines[i] is loop_, parse tags and rows until stop_.
    Returns (tags, list of token rows, index after stop_) or None.
    """
    if i >= end_bound or lines[i].strip() != "loop_":
        return None
    i += 1
    tags: List[str] = []
    while i < end_bound:
        line = lines[i].strip()
        if line.startswith("_"):
            tags.append(line)
            i += 1
        else:
            break
    if not tags:
        return None
    rows: List[List[str]] = []
    while i < end_bound:
        line = lines[i].strip()
        if line == "stop_":
            return tags, rows, i + 1
        if line == "save_":
            break
        if not line:
            i += 1
            continue
        tokens = _split_star_line(line)
        if len(tokens) == len(tags):
            rows.append(tokens)
            i += 1
        elif len(tokens) > 0:
            buf = tokens[:]
            i += 1
            while i < end_bound and len(buf) < len(tags):
                nl = lines[i].strip()
                if nl in ("stop_", "save_"):
                    break
                if nl:
                    buf.extend(_split_star_line(nl))
                i += 1
            if len(buf) >= len(tags):
                rows.append(buf[: len(tags)])
            else:
                i += 1
        else:
            i += 1
    return None


def _saveframe_body(lines: List[str], start: int, save_end: int) -> str:
    return "\n".join(lines[start + 1 : save_end])


def _is_sample_conditions_saveframe_v21(body: str) -> bool:
    return bool(
        re.search(r"^\s*_Saveframe_category\s+sample_conditions\s*$", body, re.MULTILINE)
    )


def _is_sample_conditions_saveframe_v3(body: str) -> bool:
    return "_Sample_condition_list.Sf_category" in body and re.search(
        r"^\s*_Sample_condition_list\.Sf_category\s+sample_conditions\s*$", body, re.MULTILINE
    )


def _is_sample_saveframe_v3(body: str) -> bool:
    return bool(
        re.search(r"^\s*_Sample\.Sf_category\s+sample\s*$", body, re.MULTILINE)
    )


def _merge_condition_row(
    acc: Dict[str, Tuple[float, str]], row: Dict[str, str], tag_type: str, tag_val: str, tag_units: str
) -> None:
    tkey = _normalize_condition_key(row.get(tag_type, ""))
    if not tkey:
        return
    val = _parse_float_token(row.get(tag_val, ""))
    if val is None:
        return
    units = row.get(tag_units, ".")
    if tkey not in acc:
        acc[tkey] = (val, units)


def extract_conditions_from_star(star_path: str) -> ExtractedConditions:
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(star_path)

    merged: Dict[str, Tuple[float, str]] = {}

    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if fmt == "3":
            if _is_sample_conditions_saveframe_v3(body):
                i = start + 1
                while i < save_end:
                    parsed = _parse_loop_at(lines, i, save_end)
                    if not parsed:
                        i += 1
                        continue
                    tags, rows, ni = parsed
                    i = ni
                    tag_names = [t.strip() for t in tags]
                    if not any("Sample_condition_variable.Type" in t for t in tag_names):
                        continue
                    t_type = next(
                        (t for t in tags if t.strip().endswith("Sample_condition_variable.Type")),
                        "",
                    )
                    t_val = next(
                        (t for t in tags if t.strip().endswith("Sample_condition_variable.Val")),
                        "",
                    )
                    t_units = next(
                        (t for t in tags if "Sample_condition_variable.Val_units" in t),
                        "",
                    )
                    if not (t_type and t_val and t_units):
                        continue
                    for toks in rows:
                        row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                        _merge_condition_row(merged, row, t_type, t_val, t_units)
        else:
            if _is_sample_conditions_saveframe_v21(body):
                i = start + 1
                while i < save_end:
                    parsed = _parse_loop_at(lines, i, save_end)
                    if not parsed:
                        i += 1
                        continue
                    tags, rows, ni = parsed
                    i = ni
                    if not any("_Variable_type" == t.strip() for t in tags):
                        continue
                    t_type = next(t for t in tags if t.strip() == "_Variable_type")
                    t_val = next(t for t in tags if t.strip() == "_Variable_value")
                    t_units = next(t for t in tags if t.strip() == "_Variable_value_units")
                    for toks in rows:
                        row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                        _merge_condition_row(merged, row, t_type, t_val, t_units)

    out = ExtractedConditions()
    if "ph" in merged:
        out.pH = merged["ph"][0]
    if "temperature" in merged:
        val, units = merged["temperature"]
        tc = _temp_to_celsius(val, units)
        if tc is not None:
            out.temperature_C = tc
    if "ionic_strength" in merged:
        val, units = merged["ionic_strength"]
        im = _ionic_to_mM(val, units)
        if im is not None:
            out.ionic_strength_mM = im

    # NaCl: NMR-STAR 3 sample components only
    if fmt == "3":
        nacl_done = False
        for start, save_end, _name in _iter_saveframes(lines):
            if nacl_done:
                break
            body = _saveframe_body(lines, start, save_end)
            if not _is_sample_saveframe_v3(body):
                continue
            i = start + 1
            while i < save_end:
                parsed = _parse_loop_at(lines, i, save_end)
                if not parsed:
                    i += 1
                    continue
                tags, rows, ni = parsed
                i = ni
                tag_stripped = [t.strip() for t in tags]
                if not any("Sample_component.Mol_common_name" in t for t in tag_stripped):
                    continue
                t_name = next(
                    (t for t in tags if "Sample_component.Mol_common_name" in t), ""
                )
                t_conc = next(
                    (
                        t
                        for t in tags
                        if "Sample_component.Concentration_val" in t
                        and "min" not in t.lower()
                        and "max" not in t.lower()
                        and "err" not in t.lower()
                    ),
                    "",
                )
                t_units = next(
                    (t for t in tags if "Sample_component.Concentration_val_units" in t), ""
                )
                if not (t_name and t_conc and t_units):
                    continue
                for toks in rows:
                    row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                    name_raw = row.get(t_name, "").strip().strip("'\"")
                    if "sodium chloride" not in name_raw.lower():
                        continue
                    cv = _parse_float_token(row.get(t_conc, ""))
                    if cv is None:
                        continue
                    u_raw = row.get(t_units, "mM").strip().strip("'\"")
                    c_mm = _ionic_to_mM(cv, u_raw)
                    if c_mm is not None:
                        out.nacl_mM = c_mm
                        nacl_done = True
                        break

    return out


def resolve_star_path(bmrb_id: str, cs_dir: Path) -> Optional[Path]:
    for suffix in ("_21.str", "_3.str"):
        p = cs_dir / f"{bmrb_id}{suffix}"
        if p.is_file() and p.stat().st_size > 0:
            return p
    return None


def ensure_star_path(bmrb_id: str, cs_dir: Path, do_fetch: bool) -> Optional[Path]:
    p = resolve_star_path(bmrb_id, cs_dir)
    if p is not None:
        return p
    if do_fetch:
        try:
            from scripts.bmrb_io import fetch_bmrb  # noqa: WPS433

            got = fetch_bmrb(bmrb_id, cache_dir=str(cs_dir), force=False)
            return Path(got) if got else None
        except Exception:
            return None
    return None


def _similar_salt(a: Optional[float], b: Optional[float]) -> bool:
    if a is None or b is None:
        return False
    return abs(a - b) <= NACL_TOLERANCE_MM


def _similar_ph(a: Optional[float], b: Optional[float]) -> bool:
    if a is None or b is None:
        return False
    return abs(a - b) <= PH_TOLERANCE


def _similar_temp(a: Optional[float], b: Optional[float]) -> bool:
    if a is None or b is None:
        return False
    return abs(a - b) <= TEMP_TOLERANCE_C


def row_to_csv_values(ec: ExtractedConditions) -> Dict[str, Any]:
    def fmt(x: Optional[float], nd: int = 4) -> str:
        if x is None:
            return ""
        r = round(float(x), nd)
        if nd <= 2 and r == int(r):
            return str(int(r))
        return str(r)

    return {
        "pH": fmt(ec.pH, 3),
        "temperature_C": fmt(ec.temperature_C, 2),
        "NaCl_mM": fmt(ec.nacl_mM, 4),
        "ionic_strength_mM": fmt(ec.ionic_strength_mM, 4),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Apo/holo BMRB experimental conditions table.")
    ap.add_argument(
        "--input",
        type=Path,
        default=Path(Paths().input_csv),
        help="Input CSV (default: data/CSP_UBQ.csv)",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(Paths().exp_conditions_csv),
        help="Output CSV path",
    )
    ap.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(Paths().cs_cache_dir),
        help="Directory with {bmrb}_21.str / {bmrb}_3.str",
    )
    ap.add_argument(
        "--fetch",
        action="store_true",
        help="Download missing STAR files via bmrb_io.fetch_bmrb",
    )
    args = ap.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames_in = reader.fieldnames or []
        rows = list(reader)

    out_fields = [
        "apo_bmrb",
        "holo_bmrb",
        "apo_pdb",
        "holo_pdb",
        "apo_star_path",
        "holo_star_path",
        "apo_pH",
        "apo_temperature_C",
        "apo_NaCl_mM",
        "apo_ionic_strength_mM",
        "holo_pH",
        "holo_temperature_C",
        "holo_NaCl_mM",
        "holo_ionic_strength_mM",
        "salt_similar",
        "ph_similar",
        "temp_similar",
        "conditions_similar",
    ]

    out_rows: List[Dict[str, Any]] = []
    for row in rows:
        apo_id = (row.get("apo_bmrb") or "").strip()
        holo_id = (row.get("holo_bmrb") or "").strip()
        apo_pdb = (row.get("apo_pdb") or "").strip()
        holo_pdb = (row.get("holo_pdb") or "").strip()

        apo_path = ensure_star_path(apo_id, args.cs_dir, args.fetch) if apo_id else None
        holo_path = ensure_star_path(holo_id, args.cs_dir, args.fetch) if holo_id else None

        apo_ec = extract_conditions_from_star(str(apo_path)) if apo_path else ExtractedConditions()
        holo_ec = extract_conditions_from_star(str(holo_path)) if holo_path else ExtractedConditions()

        salt_ok = _similar_salt(apo_ec.nacl_mM, holo_ec.nacl_mM)
        ph_ok = _similar_ph(apo_ec.pH, holo_ec.pH)
        temp_ok = _similar_temp(apo_ec.temperature_C, holo_ec.temperature_C)
        all_ok = salt_ok and ph_ok and temp_ok

        av = row_to_csv_values(apo_ec)
        hv = row_to_csv_values(holo_ec)
        out_rows.append(
            {
                "apo_bmrb": apo_id,
                "holo_bmrb": holo_id,
                "apo_pdb": apo_pdb,
                "holo_pdb": holo_pdb,
                "apo_star_path": str(apo_path) if apo_path else "",
                "holo_star_path": str(holo_path) if holo_path else "",
                "apo_pH": av["pH"],
                "apo_temperature_C": av["temperature_C"],
                "apo_NaCl_mM": av["NaCl_mM"],
                "apo_ionic_strength_mM": av["ionic_strength_mM"],
                "holo_pH": hv["pH"],
                "holo_temperature_C": hv["temperature_C"],
                "holo_NaCl_mM": hv["NaCl_mM"],
                "holo_ionic_strength_mM": hv["ionic_strength_mM"],
                "salt_similar": salt_ok,
                "ph_similar": ph_ok,
                "temp_similar": temp_ok,
                "conditions_similar": all_ok,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
