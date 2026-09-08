#!/usr/bin/env python3
"""
Find sequence-similar single-polypeptide apo BMRB entries for holo IDs in
CSP_UBQ.csv that are absent from CSP_UBQ_ph0.5_temp5C.csv.

Searches BMRB with the holo *receptor* sequence (longest polypeptide polymer
in the holo STAR; assigned-shift fallback). Keeps candidates that:
  - cover 100% of the query receptor sequence
  - have >= 95% sequence identity
  - contain exactly one polypeptide polymer

Writes:
  data/CSP_UBQ_excluded_holo_apo_sequence_matches.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.apo_holo_exp_conditions import (  # noqa: E402
    PH_TOLERANCE,
    TEMP_TOLERANCE_C,
    _detect_bmrb_format,
    _iter_saveframes,
    _saveframe_body,
    _similar_ph,
    _similar_temp,
    _tokenize_star_lines,
    ensure_star_path,
    extract_conditions_from_star,
)
from scripts.config import Paths  # noqa: E402
from scripts.find_apo_apo_matches import (  # noqa: E402
    count_polypeptide_polymers,
    filter_fasta_hits,
    search_bmrb_fasta,
    write_csv,
)

PAIR_FIELDS = [
    "holo_bmrb",
    "apo_bmrb",
    "percent_id",
    "query_coverage",
    "alignment_length",
    "holo_seq",
    "apo_seq",
    "holo_seq_len",
    "apo_seq_len",
    "current_apo_bmrb",
    "holo_pdb",
    "candidate_is_existing_holo",
    "holo_star_path",
    "apo_star_path",
    "apo_pH",
    "apo_temperature_C",
    "apo_pressure_atm",
    "holo_pH",
    "holo_temperature_C",
    "holo_pressure_atm",
    "delta_pH",
    "delta_T_C",
    "ph_similar",
    "temp_similar",
    "ph_temp_match",
]


def _clean_seq(raw: str) -> str:
    return re.sub(r"[^A-Za-z]", "", raw or "").upper()


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_REPO_ROOT))
    except ValueError:
        return str(path)


def _fmt_cond(x: Optional[float], nd: int) -> str:
    if x is None:
        return ""
    r = round(float(x), nd)
    if nd <= 2 and r == int(r):
        return str(int(r))
    return str(r)


def pair_condition_fields(apo_star: Path, holo_star: Path) -> Dict[str, str]:
    """pH, °C temperature, pressure in atm, and ±0.5 / ±5 °C match flags."""
    apo_ec = extract_conditions_from_star(str(apo_star))
    holo_ec = extract_conditions_from_star(str(holo_star))
    ph_ok = _similar_ph(apo_ec.pH, holo_ec.pH)
    temp_ok = _similar_temp(apo_ec.temperature_C, holo_ec.temperature_C)
    d_ph = ""
    if apo_ec.pH is not None and holo_ec.pH is not None:
        d_ph = _fmt_cond(abs(apo_ec.pH - holo_ec.pH), 3)
    d_t = ""
    if apo_ec.temperature_C is not None and holo_ec.temperature_C is not None:
        d_t = _fmt_cond(abs(apo_ec.temperature_C - holo_ec.temperature_C), 2)
    return {
        "apo_pH": _fmt_cond(apo_ec.pH, 3),
        "apo_temperature_C": _fmt_cond(apo_ec.temperature_C, 2),
        "apo_pressure_atm": _fmt_cond(apo_ec.pressure_atm, 4),
        "holo_pH": _fmt_cond(holo_ec.pH, 3),
        "holo_temperature_C": _fmt_cond(holo_ec.temperature_C, 2),
        "holo_pressure_atm": _fmt_cond(holo_ec.pressure_atm, 4),
        "delta_pH": d_ph,
        "delta_T_C": d_t,
        "ph_similar": str(ph_ok),
        "temp_similar": str(temp_ok),
        "ph_temp_match": str(ph_ok and temp_ok),
    }


def annotate_pairs_with_conditions(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Fill condition columns on existing pair rows using cached STAR paths."""
    cache: Dict[str, Any] = {}
    out: List[Dict[str, Any]] = []
    for row in rows:
        apo_path = (row.get("apo_star_path") or "").strip()
        holo_path = (row.get("holo_star_path") or "").strip()
        updated = dict(row)
        if apo_path and holo_path:
            key = (apo_path, holo_path)
            if key not in cache:
                apo_p = Path(apo_path)
                holo_p = Path(holo_path)
                if not apo_p.is_absolute():
                    apo_p = _REPO_ROOT / apo_p
                if not holo_p.is_absolute():
                    holo_p = _REPO_ROOT / holo_p
                cache[key] = pair_condition_fields(apo_p, holo_p)
            updated.update(cache[key])
        out.append(updated)
    return out


def _join_unique(values: Sequence[str]) -> str:
    seen: Set[str] = set()
    out: List[str] = []
    for v in values:
        s = (v or "").strip()
        if not s or s in seen:
            continue
        seen.add(s)
        out.append(s)
    return ";".join(out)


def _one_letter_seq_from_body(body: str) -> str:
    """Extract a polymer one-letter sequence from a single saveframe body."""
    tag = (
        r"(?:_Entity\.Polymer_seq_one_letter_code|"
        r"_Mol_residue_sequence|"
        r"_Polymer_seq_one_letter_code)"
    )
    m = re.search(rf"{tag}\s+(\S+)", body)
    if m:
        tok = m.group(1).strip().strip("'\"")
        if tok and tok not in (".", "?") and not tok.startswith(";"):
            seq = _clean_seq(tok)
            if len(seq) >= 5:
                return seq
    m = re.search(rf"{tag}\s*\n\s*;\s*\n?(.*?)\n\s*;", body, re.DOTALL)
    if m:
        seq = _clean_seq(m.group(1))
        if len(seq) >= 5:
            return seq
    return ""


def _is_v3_polypeptide_entity(body: str) -> bool:
    if not re.search(r"^\s*_Entity\.Type\s+polymer\s*$", body, re.MULTILINE):
        return False
    return bool(
        re.search(r"^\s*_Entity\.Polymer_type\s+polypeptide", body, re.MULTILINE | re.IGNORECASE)
    )


def polypeptide_sequences(star_path: Path) -> List[str]:
    """One-letter sequences for every polypeptide polymer entity in a STAR file."""
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(str(star_path))
    seqs: List[str] = []

    if fmt == "3":
        for start, save_end, _name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not _is_v3_polypeptide_entity(body):
                continue
            seq = _one_letter_seq_from_body(body)
            if seq:
                seqs.append(seq)
        return seqs

    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if not re.search(r"^\s*_Saveframe_category\s+monomeric_polymer\s*$", body, re.MULTILINE):
            continue
        seq = _one_letter_seq_from_body(body)
        if seq:
            seqs.append(seq)
    return seqs


def _assigned_shift_sequences(star_path: Path) -> List[str]:
    try:
        from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # noqa: WPS433

        frames = parse_sequence_and_shifts_from_saveframes(str(star_path))
    except Exception:
        return []
    seqs: List[str] = []
    for entry in frames:
        if not entry:
            continue
        seq = _clean_seq(str(entry[0] or ""))
        if seq:
            seqs.append(seq)
    return seqs


def extract_receptor_sequence(star_path: Path) -> str:
    """Longest polypeptide polymer sequence; assigned-shift fallback."""
    polymer = polypeptide_sequences(star_path)
    if polymer:
        return max(polymer, key=len)
    assigned = _assigned_shift_sequences(star_path)
    if assigned:
        return max(assigned, key=len)
    return ""


def extract_apo_sequence(star_path: Path) -> str:
    """Sequence of the single polypeptide (or longest if parsing is ambiguous)."""
    return extract_receptor_sequence(star_path)


def load_excluded_holos(
    input_csv: Path,
    exclude_csv: Path,
    id_filter: Optional[Set[str]],
) -> Tuple[List[str], Dict[str, Dict[str, str]], Set[str]]:
    """Return excluded holo IDs, per-holo context, and all CSP_UBQ holo IDs."""
    all_holos: Set[str] = set()
    context: Dict[str, Dict[str, List[str]]] = {}
    with open(input_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            holo = (row.get("holo_bmrb") or "").strip()
            if not holo:
                continue
            all_holos.add(holo)
            bucket = context.setdefault(holo, {"apo": [], "pdb": []})
            bucket["apo"].append((row.get("apo_bmrb") or "").strip())
            bucket["pdb"].append((row.get("holo_pdb") or "").strip())

    excluded_holos: Set[str] = set()
    with open(exclude_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            holo = (row.get("holo_bmrb") or "").strip()
            if holo:
                excluded_holos.add(holo)

    holos = [h for h in all_holos if h not in excluded_holos]
    if id_filter:
        holos = [h for h in holos if h in id_filter]
    holos.sort(key=lambda x: (len(x), x))

    joined: Dict[str, Dict[str, str]] = {
        h: {
            "current_apo_bmrb": _join_unique(context[h]["apo"]),
            "holo_pdb": _join_unique(context[h]["pdb"]),
        }
        for h in holos
    }
    return holos, joined, all_holos


def process_holo(
    holo_id: str,
    *,
    ctx: Dict[str, str],
    cs_dir: Path,
    do_fetch: bool,
    min_identity: float,
    min_coverage: float,
    all_holo_ids: Set[str],
    max_candidates: Optional[int],
) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    """Return (accepted pair rows, unmatched_reason or None)."""
    holo_star = ensure_star_path(holo_id, cs_dir, do_fetch)
    if holo_star is None:
        return [], "query_star_missing"

    holo_seq = extract_receptor_sequence(holo_star)
    if not holo_seq:
        return [], "query_sequence_empty"

    hits = search_bmrb_fasta(holo_seq)
    filtered = filter_fasta_hits(
        hits,
        query_id=holo_id,
        query_len=len(holo_seq),
        min_identity=min_identity,
        min_coverage=min_coverage,
    )
    if max_candidates is not None:
        filtered = filtered[:max_candidates]

    print(
        f"[HOLO] {holo_id}: seq_len={len(holo_seq)} fasta_hits={len(hits)} "
        f"after_filter={len(filtered)}",
        flush=True,
    )
    if not filtered:
        return [], "no_fasta_hits_passing_sequence_filter"

    pairs: List[Dict[str, Any]] = []
    for hit in filtered:
        cand_star = ensure_star_path(hit.entry_id, cs_dir, do_fetch)
        if cand_star is None:
            print(
                f"[HOLO] {holo_id}: skip {hit.entry_id} (star missing)",
                flush=True,
            )
            continue
        n_poly = count_polypeptide_polymers(cand_star)
        if n_poly != 1:
            print(
                f"[HOLO] {holo_id}: skip {hit.entry_id} (n_polypeptide={n_poly})",
                flush=True,
            )
            continue
        apo_seq = extract_apo_sequence(cand_star)
        if not apo_seq:
            print(
                f"[HOLO] {holo_id}: skip {hit.entry_id} (apo sequence empty)",
                flush=True,
            )
            continue
        cov = hit.query_coverage(len(holo_seq))
        row = {
            "holo_bmrb": holo_id,
            "apo_bmrb": hit.entry_id,
            "percent_id": f"{hit.percent_id:.2f}",
            "query_coverage": f"{cov:.4f}",
            "alignment_length": str(hit.alignment_length),
            "holo_seq": holo_seq,
            "apo_seq": apo_seq,
            "holo_seq_len": str(len(holo_seq)),
            "apo_seq_len": str(len(apo_seq)),
            "current_apo_bmrb": ctx.get("current_apo_bmrb", ""),
            "holo_pdb": ctx.get("holo_pdb", ""),
            "candidate_is_existing_holo": (
                "true" if hit.entry_id in all_holo_ids else "false"
            ),
            "holo_star_path": _relpath(holo_star),
            "apo_star_path": _relpath(cand_star),
        }
        row.update(pair_condition_fields(cand_star, holo_star))
        pairs.append(row)

    if not pairs:
        return [], "no_single_polypeptide_candidates"
    return pairs, None


def main(argv: Optional[Sequence[str]] = None) -> int:
    cfg = Paths()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input",
        type=Path,
        default=Path(cfg.input_csv),
        help="Full CSP_UBQ CSV (default: data/CSP_UBQ.csv)",
    )
    ap.add_argument(
        "--exclude",
        type=Path,
        default=Path(cfg.filtered_input_csv),
        help="CSV whose holo_bmrb IDs are skipped (default: data/CSP_UBQ_ph0.5_temp5C.csv)",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(cfg.data_dir) / "CSP_UBQ_excluded_holo_apo_sequence_matches.csv",
        help="Output pairs CSV",
    )
    ap.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(cfg.cs_cache_dir),
        help="STAR cache directory (default: CS_Lists)",
    )
    ap.add_argument(
        "--fetch",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Download missing STAR files (default: on)",
    )
    ap.add_argument(
        "--ids",
        type=str,
        default="",
        help="Comma-separated holo_bmrb IDs to process (default: all excluded)",
    )
    ap.add_argument("--min-identity", type=float, default=95.0, help="Min FASTA %% identity")
    ap.add_argument(
        "--min-coverage",
        type=float,
        default=1.0,
        help="Min query coverage (aligned query span / query length)",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for per-holo processing (default: 1)",
    )
    ap.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="Cap FASTA candidates evaluated per query (after sequence filter)",
    )
    ap.add_argument(
        "--annotate-only",
        action="store_true",
        help="Skip FASTA search; add pH/temp/pressure columns to an existing pairs CSV",
    )
    ap.add_argument(
        "--match-output",
        type=Path,
        default=None,
        help="CSV of pairs passing |ΔpH|≤0.5 and |ΔT|≤5 °C "
        "(default: data/CSP_UBQ_excluded_holo_apo_ph0.5_temp5C_matches.csv)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    match_out = args.match_output
    if match_out is None:
        match_out = Path(cfg.data_dir) / "CSP_UBQ_excluded_holo_apo_ph0.5_temp5C_matches.csv"

    if args.annotate_only:
        with open(args.output, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
        if not existing:
            print(f"No rows in {args.output}", file=sys.stderr)
            return 1
        annotated = annotate_pairs_with_conditions(existing)
        write_csv(args.output, PAIR_FIELDS, annotated)
        matches = [r for r in annotated if r.get("ph_temp_match") == "True"]
        write_csv(match_out, PAIR_FIELDS, matches)
        print(
            f"Annotated {len(annotated)} pairs -> {args.output}; "
            f"{len(matches)} ph/temp matches -> {match_out}",
            flush=True,
        )
        return 0

    id_filter: Optional[Set[str]] = None
    if args.ids.strip():
        id_filter = {x.strip() for x in args.ids.split(",") if x.strip()}

    holos, context, all_holo_ids = load_excluded_holos(args.input, args.exclude, id_filter)
    if not holos:
        print("No excluded holo_bmrb IDs to process.", file=sys.stderr)
        return 1

    print(f"Processing {len(holos)} excluded holo_bmrb IDs", flush=True)

    all_pairs: List[Dict[str, Any]] = []
    unmatched: List[Tuple[str, str]] = []

    def _run_one(holo_id: str) -> Tuple[str, List[Dict[str, Any]], Optional[str]]:
        pairs, reason = process_holo(
            holo_id,
            ctx=context.get(holo_id, {}),
            cs_dir=args.cs_dir,
            do_fetch=args.fetch,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            all_holo_ids=all_holo_ids,
            max_candidates=args.max_candidates,
        )
        return holo_id, pairs, reason

    if args.workers <= 1:
        results = [_run_one(h) for h in holos]
    else:
        results = []
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_run_one, h): h for h in holos}
            for fut in as_completed(futs):
                holo_id = futs[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:
                    print(f"[ERROR] {holo_id}: {exc}", file=sys.stderr)
                    unmatched.append((holo_id, f"exception:{exc}"))

    for holo_id, pairs, reason in results:
        all_pairs.extend(pairs)
        if reason:
            unmatched.append((holo_id, reason))

    all_pairs.sort(key=lambda r: (r.get("holo_bmrb", ""), r.get("apo_bmrb", "")))
    write_csv(args.output, PAIR_FIELDS, all_pairs)
    matches = [r for r in all_pairs if r.get("ph_temp_match") == "True"]
    write_csv(match_out, PAIR_FIELDS, matches)

    matched_holos = {r["holo_bmrb"] for r in all_pairs}
    print(
        f"Wrote {len(all_pairs)} pairs for {len(matched_holos)} holos -> {args.output}",
        flush=True,
    )
    print(
        f"Wrote {len(matches)} pairs with |ΔpH|≤{PH_TOLERANCE} and "
        f"|ΔT|≤{TEMP_TOLERANCE_C} °C -> {match_out}",
        flush=True,
    )
    if unmatched:
        print(f"Unmatched holos: {len(unmatched)}", flush=True)
        for holo_id, reason in unmatched:
            print(f"  {holo_id}: {reason}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
