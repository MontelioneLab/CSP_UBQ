#!/usr/bin/env python3
"""
Find equivalent apo BMRB entries for each apo in CSP_UBQ.csv.

For each unique apo_bmrb, searches BMRB via the FASTA API for sequence-similar
entries, then keeps candidates that:
  - have exactly one polypeptide polymer
  - have PubMed ID or DOI provenance
  - match experimental conditions (pH, temperature, ionic strength, spectrometer)

Writes:
  outputs/apo_apo_matches/apo_apo_pairs.csv
  outputs/apo_apo_matches/apo_apo_candidates_audit.csv
  outputs/apo_apo_matches/shifts/{id}_*.str (+ parsed residue tables when possible)

Re-extract pH, temperature, sample buffer, and entity names onto an existing
pairs CSV with ``--enrich path/to/apo_apo_pairs.csv``.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.apo_holo_exp_conditions import (  # noqa: E402
    IONIC_TOLERANCE_MM,
    PH_TOLERANCE,
    TEMP_TOLERANCE_C,
    ExtractedConditions,
    _detect_bmrb_format,
    _ionic_to_mM,
    _iter_saveframes,
    _parse_float_token,
    _parse_loop_at,
    _saveframe_body,
    _similar_spectrometer,
    _tokenize_star_lines,
    ensure_star_path,
    extract_conditions_from_star,
    row_to_csv_values,
)
from scripts.config import Paths  # noqa: E402

BMRB_FASTA_URL = "https://api.bmrb.io/v2/search/fasta/{sequence}"

PAIR_FIELDS = [
    "query_apo_bmrb",
    "match_apo_bmrb",
    "percent_id",
    "query_coverage",
    "alignment_length",
    "match_is_csp_holo",
    "query_star_path",
    "match_star_path",
    "query_pubmed",
    "query_doi",
    "match_pubmed",
    "match_doi",
    "query_buffer",
    "match_buffer",
    "query_entities",
    "match_entities",
    "query_pH",
    "match_pH",
    "delta_pH",
    "query_temperature_C",
    "match_temperature_C",
    "delta_T_C",
    "query_ionic_strength_mM",
    "match_ionic_strength_mM",
    "delta_ionic_mM",
    "query_field_strength_MHz",
    "match_field_strength_MHz",
    "query_field_strength_Hz",
    "match_field_strength_Hz",
    "ph_ok",
    "temp_ok",
    "ionic_ok",
    "spectrometer_ok",
]

AUDIT_FIELDS = [
    "query_apo_bmrb",
    "candidate_bmrb",
    "percent_id",
    "query_coverage",
    "status",
    "reason",
]


@dataclass
class FastaHit:
    entry_id: str
    entity_id: str
    percent_id: float
    alignment_length: int
    q_start: int
    q_end: int
    e_value: Optional[float] = None
    bit_score: Optional[float] = None

    def query_coverage(self, query_len: int) -> float:
        if query_len <= 0:
            return 0.0
        return (self.q_end - self.q_start + 1) / float(query_len)


@dataclass
class Provenance:
    pubmed: str = ""
    doi: str = ""

    @property
    def adequate(self) -> bool:
        return bool(self.pubmed) or bool(self.doi)


def _fmt_num(x: Optional[float], nd: int = 4) -> str:
    if x is None:
        return ""
    r = round(float(x), nd)
    if nd <= 1 and r == int(r):
        return str(int(r))
    return str(r)


def _mhz_to_hz(mhz: Optional[float]) -> Optional[float]:
    if mhz is None:
        return None
    return float(mhz) * 1_000_000.0


def _strip_star_token(tok: str) -> str:
    t = tok.strip().strip("'\"")
    if t in (".", "?", ""):
        return ""
    return t


def load_unique_apo_ids(input_csv: Path, id_filter: Optional[Set[str]]) -> Tuple[List[str], Set[str]]:
    """Return sorted unique apo_bmrb IDs and the set of holo_bmrb IDs."""
    apos: Set[str] = set()
    holos: Set[str] = set()
    with open(input_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            apo = (row.get("apo_bmrb") or "").strip()
            holo = (row.get("holo_bmrb") or "").strip()
            if apo:
                apos.add(apo)
            if holo:
                holos.add(holo)
    if id_filter:
        apos = {a for a in apos if a in id_filter}
    return sorted(apos, key=lambda x: (len(x), x)), holos


def _polymer_residue_sequence(star_path: Path) -> str:
    """Extract polymer residue sequence from monomeric_polymer / Entity polymer tags."""
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)

    # Prefer _Mol_residue_sequence / _Entity.Polymer_seq_one_letter_code style tags
    for pat in (
        r"_Mol_residue_sequence\s+(\S+)",
        r"_Entity\.Polymer_seq_one_letter_code\s+(\S+)",
        r"_Polymer_seq_one_letter_code\s+(\S+)",
    ):
        m = re.search(pat, text)
        if m:
            tok = _strip_star_token(m.group(1))
            if tok and not tok.startswith(";"):
                seq = re.sub(r"[^A-Za-z]", "", tok).upper()
                if len(seq) >= 5:
                    return seq

    # Multi-line semicolon blocks for polymer sequence
    for pat in (
        r"_Mol_residue_sequence\s*\n\s*;\s*\n?(.*?)\n\s*;",
        r"_Entity\.Polymer_seq_one_letter_code\s*\n\s*;\s*\n?(.*?)\n\s*;",
    ):
        m = re.search(pat, text, re.DOTALL)
        if m:
            seq = re.sub(r"[^A-Za-z]", "", m.group(1)).upper()
            if len(seq) >= 5:
                return seq

    try:
        from scripts.bmrb_io import _extract_sequence_from_saveframe  # noqa: WPS433

        seq = (_extract_sequence_from_saveframe(lines) or "").strip().upper()
        return re.sub(r"[^A-Z]", "", seq)
    except Exception:
        return ""


def extract_query_sequence(star_path: Path) -> str:
    """Prefer longest H+N assigned-shift sequence; fall back to polymer residue sequence."""
    best = ""
    try:
        from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # noqa: WPS433

        frames = parse_sequence_and_shifts_from_saveframes(str(star_path))
        for seq, _h, _n, _ca, _ha, _name in frames:
            s = re.sub(r"[^A-Za-z]", "", (seq or "")).upper()
            if len(s) > len(best):
                best = s
    except Exception:
        pass
    if best:
        return best
    return _polymer_residue_sequence(star_path)


def search_bmrb_fasta(sequence: str, timeout: float = 60.0) -> List[FastaHit]:
    """Query BMRB FASTA search API; return parsed hits."""
    import json

    seq = re.sub(r"[^A-Za-z]", "", sequence)
    if not seq:
        return []
    url = BMRB_FASTA_URL.format(sequence=urllib.parse.quote(seq, safe=""))
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "CSP_UBQ-apo-apo/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
        print(f"[WARN] FASTA search failed: {exc}", file=sys.stderr)
        return []

    hits: List[FastaHit] = []
    if not isinstance(payload, list):
        return hits
    for item in payload:
        if not isinstance(item, dict):
            continue
        try:
            entry_id = str(item.get("entry_id", "")).strip()
            if not entry_id:
                continue
            percent_id = float(item.get("percent_id"))
            aln = int(float(item.get("alignment_length")))
            q_start = int(float(item.get("q.start")))
            q_end = int(float(item.get("q.end")))
        except (TypeError, ValueError, KeyError):
            continue
        e_val = None
        bit = None
        try:
            if item.get("e-value") is not None:
                e_val = float(item["e-value"])
        except (TypeError, ValueError):
            pass
        try:
            if item.get("bit_score") is not None:
                bit = float(item["bit_score"])
        except (TypeError, ValueError):
            pass
        hits.append(
            FastaHit(
                entry_id=entry_id,
                entity_id=str(item.get("entity_id", "")).strip(),
                percent_id=percent_id,
                alignment_length=aln,
                q_start=q_start,
                q_end=q_end,
                e_value=e_val,
                bit_score=bit,
            )
        )
    return hits


def filter_fasta_hits(
    hits: Sequence[FastaHit],
    *,
    query_id: str,
    query_len: int,
    min_identity: float,
    min_coverage: float,
) -> List[FastaHit]:
    """Deduplicate by entry_id (best percent_id), drop self, apply thresholds."""
    best: Dict[str, FastaHit] = {}
    for h in hits:
        if h.entry_id == query_id:
            continue
        if h.percent_id < min_identity:
            continue
        if h.query_coverage(query_len) < min_coverage:
            continue
        prev = best.get(h.entry_id)
        if prev is None or h.percent_id > prev.percent_id or (
            h.percent_id == prev.percent_id and h.alignment_length > prev.alignment_length
        ):
            best[h.entry_id] = h
    return sorted(best.values(), key=lambda x: (-x.percent_id, -x.alignment_length, x.entry_id))


def count_polypeptide_polymers(star_path: Path) -> int:
    """Count polypeptide polymer entities in an NMR-STAR file."""
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(str(star_path))

    if fmt == "3":
        n = 0
        for start, save_end, _name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not re.search(
                r"^\s*_Entity\.Sf_category\s+entity\s*$", body, re.MULTILINE
            ):
                continue
            if not re.search(r"^\s*_Entity\.Type\s+polymer\s*$", body, re.MULTILINE):
                continue
            if re.search(r"^\s*_Entity\.Polymer_type\s+polypeptide", body, re.MULTILINE | re.IGNORECASE):
                n += 1
        return n

    # NMR-STAR 2.1: monomeric_polymer saveframes (protein polymers)
    n = 0
    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if re.search(
            r"^\s*_Saveframe_category\s+monomeric_polymer\s*$", body, re.MULTILINE
        ):
            n += 1
    return n


def extract_provenance(star_path: Path) -> Provenance:
    """Extract first usable PubMed ID and DOI from citation / entry tags."""
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    pubmed = ""
    doi = ""

    for pat in (
        r"(?:_Citation\.PubMed_ID|_PubMed_ID)\s+(\S+)",
        r"(?:_Citation\.PubMed_ID|_PubMed_ID)\s+\n\s*;\s*([^\n;]+)\s*;",
    ):
        for m in re.finditer(pat, text):
            val = _strip_star_token(m.group(1))
            if val and val.lower() not in ("none", "null", "n/a"):
                # numeric PubMed preferred
                if val.isdigit() or re.match(r"^\d+$", val):
                    pubmed = val
                    break
                if not pubmed:
                    pubmed = val
        if pubmed:
            break

    for pat in (
        r"(?:_Citation\.DOI|_DOI)\s+(\S+)",
        r"(?:_Citation\.DOI|_DOI)\s+\n\s*;\s*([^\n;]+)\s*;",
    ):
        for m in re.finditer(pat, text):
            val = _strip_star_token(m.group(1))
            if val and val.lower() not in ("none", "null", "n/a"):
                doi = val
                break
        if doi:
            break

    return Provenance(pubmed=pubmed, doi=doi)


def _mol_label_display_map(lines: List[str]) -> Dict[str, str]:
    """Map v2.1 $Mol_label tokens to molecular_system component names."""
    mapping: Dict[str, str] = {}
    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if not re.search(
            r"^\s*_Saveframe_category\s+molecular_system\s*$", body, re.MULTILINE
        ):
            continue
        i = start + 1
        while i < save_end:
            parsed = _parse_loop_at(lines, i, save_end)
            if not parsed:
                i += 1
                continue
            tags, rows, ni = parsed
            i = ni
            stripped = [t.strip() for t in tags]
            if "_Mol_system_component_name" not in stripped or "_Mol_label" not in stripped:
                continue
            t_name = next(t for t in tags if t.strip() == "_Mol_system_component_name")
            t_lab = next(t for t in tags if t.strip() == "_Mol_label")
            for toks in rows:
                row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                name = _strip_star_token(row.get(t_name, ""))
                lab = _strip_star_token(row.get(t_lab, "")).lstrip("$")
                if lab and name:
                    mapping[lab] = name
    return mapping


def _format_buffer_part(name: str, conc_tok: str, units_tok: str) -> str:
    name = _strip_star_token(name)
    if not name:
        return ""
    cv = _parse_float_token(conc_tok)
    u_raw = _strip_star_token(units_tok)
    if cv is None:
        return name
    mm = _ionic_to_mM(cv, u_raw) if u_raw else None
    if mm is not None:
        return f"{name}={_fmt_num(mm, 4)} mM"
    return f"{name}={_fmt_num(cv, 4)} {u_raw}".strip()


def extract_entity_names(star_path: Path) -> str:
    """Join polymer then ligand entity display names (STAR inventory order)."""
    # Lazy import: apo_polymer_inventory imports count_polypeptide_polymers from here.
    from scripts.apo_polymer_inventory import inventory_entities  # noqa: WPS433

    _fmt, entities = inventory_entities(star_path)
    if any((e.entity_id or "").strip() for e in entities):

        def _sort_key(e: Any) -> Tuple[int, int]:
            try:
                eid = int(float(e.entity_id)) if e.entity_id else 10**9
            except (TypeError, ValueError):
                eid = 10**9
            return (eid, 0)

        ordered = sorted(enumerate(entities), key=lambda iv: (_sort_key(iv[1])[0], iv[0]))
        entities = [e for _, e in ordered]

    names: List[str] = []
    seen = set()
    for e in entities:
        name = (e.name or "").strip()
        if not name or name in (".", "?"):
            continue
        if name in seen:
            continue
        seen.add(name)
        names.append(name)
    return "; ".join(names)


def extract_buffer_text(star_path: Path) -> str:
    """Serialize sample-component buffer conditions as name=conc units; ..."""
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(str(star_path))
    parts: List[str] = []

    if fmt == "3":
        for start, save_end, _name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not re.search(r"^\s*_Sample\.Sf_category\s+sample\s*$", body, re.MULTILINE):
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
                t_name = next((t for t in tags if "Sample_component.Mol_common_name" in t), "")
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
                if not t_name:
                    continue
                for toks in rows:
                    row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                    part = _format_buffer_part(
                        row.get(t_name, ""),
                        row.get(t_conc, "") if t_conc else "",
                        row.get(t_units, "") if t_units else "",
                    )
                    if part:
                        parts.append(part)
        return "; ".join(parts)

    label_map = _mol_label_display_map(lines)
    # v2.1: sample loops use _Mol_label + _Concentration_value (not _Mol_common_name).
    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if not re.search(r"^\s*_Saveframe_category\s+sample\s*$", body, re.MULTILINE):
            continue
        i = start + 1
        while i < save_end:
            parsed = _parse_loop_at(lines, i, save_end)
            if not parsed:
                i += 1
                continue
            tags, rows, ni = parsed
            i = ni
            name_tag = next(
                (
                    t
                    for t in tags
                    if t.strip()
                    in (
                        "_Mol_label",
                        "_Mol_common_name",
                        "_Molecule_name",
                        "_Chemical_name",
                    )
                    or "Mol_common_name" in t
                ),
                "",
            )
            if not name_tag:
                continue
            conc_tag = next(
                (
                    t
                    for t in tags
                    if t.strip()
                    in (
                        "_Concentration_value",
                        "_Concentration_val",
                        "_Concentration",
                    )
                    or (
                        "Concentration_val" in t
                        and "min" not in t.lower()
                        and "max" not in t.lower()
                        and "err" not in t.lower()
                    )
                ),
                "",
            )
            units_tag = next(
                (
                    t
                    for t in tags
                    if t.strip()
                    in (
                        "_Concentration_value_units",
                        "_Concentration_val_units",
                        "_Concentration_units",
                    )
                    or "Concentration_val_units" in t
                    or t.strip() == "_Concentration_value_units"
                ),
                "",
            )
            for toks in rows:
                row = {tags[j]: toks[j] for j in range(min(len(tags), len(toks)))}
                raw_name = _strip_star_token(row.get(name_tag, ""))
                if not raw_name:
                    continue
                lab = raw_name.lstrip("$")
                name = label_map.get(lab, raw_name if not raw_name.startswith("$") else lab)
                part = _format_buffer_part(
                    name,
                    row.get(conc_tag, "") if conc_tag else "",
                    row.get(units_tag, "") if units_tag else "",
                )
                if part:
                    parts.append(part)
    return "; ".join(parts)


def conditions_match(
    query: ExtractedConditions,
    cand: ExtractedConditions,
    *,
    ph_tol: float,
    temp_tol: float,
    ionic_tol: float,
) -> Tuple[bool, Dict[str, Any]]:
    """Require pH, T, ionic, and field present and within tolerances."""
    d_ph = (
        abs(query.pH - cand.pH)
        if query.pH is not None and cand.pH is not None
        else None
    )
    d_t = (
        abs(query.temperature_C - cand.temperature_C)
        if query.temperature_C is not None and cand.temperature_C is not None
        else None
    )
    d_i = (
        abs(query.ionic_strength_mM - cand.ionic_strength_mM)
        if query.ionic_strength_mM is not None and cand.ionic_strength_mM is not None
        else None
    )
    ph_ok = d_ph is not None and d_ph <= ph_tol
    temp_ok = d_t is not None and d_t <= temp_tol
    ionic_ok = d_i is not None and d_i <= ionic_tol
    spec_ok = _similar_spectrometer(query.field_strength_MHz, cand.field_strength_MHz)
    keep = ph_ok and temp_ok and ionic_ok and spec_ok
    detail = {
        "delta_pH": d_ph,
        "delta_T_C": d_t,
        "delta_ionic_mM": d_i,
        "ph_ok": ph_ok,
        "temp_ok": temp_ok,
        "ionic_ok": ionic_ok,
        "spectrometer_ok": spec_ok,
        "missing_ph": query.pH is None or cand.pH is None,
        "missing_temp": query.temperature_C is None or cand.temperature_C is None,
        "missing_ionic": query.ionic_strength_mM is None or cand.ionic_strength_mM is None,
        "missing_field": query.field_strength_MHz is None or cand.field_strength_MHz is None,
    }
    return keep, detail


def copy_star_to_shifts(src: Path, shifts_dir: Path) -> Path:
    """Copy STAR into shifts_dir; return destination path."""
    shifts_dir.mkdir(parents=True, exist_ok=True)
    dest = shifts_dir / src.name
    if dest.resolve() != src.resolve():
        shutil.copy2(src, dest)
    return dest


def export_parsed_shifts(star_path: Path) -> None:
    """Best-effort residue-shift CSV export next to STAR (CS_Lists/parsed or shifts/parsed)."""
    try:
        from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # noqa: WPS433

        parse_sequence_and_shifts_from_saveframes(str(star_path))
    except Exception as exc:
        print(f"[WARN] parsed export failed for {star_path.name}: {exc}", file=sys.stderr)


def ensure_candidate_star(
    bmrb_id: str,
    cs_dir: Path,
    do_fetch: bool,
) -> Optional[Path]:
    """Resolve/fetch STAR into CS_Lists cache (not yet copied to shifts/)."""
    return ensure_star_path(bmrb_id, cs_dir, do_fetch)


def build_pair_row(
    *,
    query_id: str,
    hit: FastaHit,
    query_len: int,
    query_star: Path,
    match_star: Path,
    query_prov: Provenance,
    match_prov: Provenance,
    query_ec: ExtractedConditions,
    match_ec: ExtractedConditions,
    query_buffer: str,
    match_buffer: str,
    cond_detail: Dict[str, Any],
    holo_ids: Set[str],
) -> Dict[str, Any]:
    qv = row_to_csv_values(query_ec)
    mv = row_to_csv_values(match_ec)
    q_hz = _mhz_to_hz(query_ec.field_strength_MHz)
    m_hz = _mhz_to_hz(match_ec.field_strength_MHz)
    return {
        "query_apo_bmrb": query_id,
        "match_apo_bmrb": hit.entry_id,
        "percent_id": _fmt_num(hit.percent_id, 2),
        "query_coverage": _fmt_num(hit.query_coverage(query_len), 4),
        "alignment_length": str(hit.alignment_length),
        "match_is_csp_holo": "True" if hit.entry_id in holo_ids else "False",
        "query_star_path": str(query_star),
        "match_star_path": str(match_star),
        "query_pubmed": query_prov.pubmed,
        "query_doi": query_prov.doi,
        "match_pubmed": match_prov.pubmed,
        "match_doi": match_prov.doi,
        "query_buffer": query_buffer,
        "match_buffer": match_buffer,
        "query_entities": extract_entity_names(query_star),
        "match_entities": extract_entity_names(match_star),
        "query_pH": qv["pH"],
        "match_pH": mv["pH"],
        "delta_pH": _fmt_num(cond_detail.get("delta_pH"), 3),
        "query_temperature_C": qv["temperature_C"],
        "match_temperature_C": mv["temperature_C"],
        "delta_T_C": _fmt_num(cond_detail.get("delta_T_C"), 2),
        "query_ionic_strength_mM": qv["ionic_strength_mM"],
        "match_ionic_strength_mM": mv["ionic_strength_mM"],
        "delta_ionic_mM": _fmt_num(cond_detail.get("delta_ionic_mM"), 4),
        "query_field_strength_MHz": qv["field_strength_MHz"],
        "match_field_strength_MHz": mv["field_strength_MHz"],
        "query_field_strength_Hz": _fmt_num(q_hz, 1),
        "match_field_strength_Hz": _fmt_num(m_hz, 1),
        "ph_ok": str(bool(cond_detail.get("ph_ok"))),
        "temp_ok": str(bool(cond_detail.get("temp_ok"))),
        "ionic_ok": str(bool(cond_detail.get("ionic_ok"))),
        "spectrometer_ok": str(bool(cond_detail.get("spectrometer_ok"))),
    }


def evaluate_candidate(
    *,
    query_id: str,
    hit: FastaHit,
    query_len: int,
    query_star: Path,
    query_ec: ExtractedConditions,
    query_prov: Provenance,
    query_buffer: str,
    cs_dir: Path,
    shifts_dir: Path,
    do_fetch: bool,
    export_parsed: bool,
    ph_tol: float,
    temp_tol: float,
    ionic_tol: float,
    holo_ids: Set[str],
) -> Tuple[Optional[Dict[str, Any]], Dict[str, Any]]:
    """Return (pair_row_or_None, audit_row)."""
    cov = hit.query_coverage(query_len)
    audit_base = {
        "query_apo_bmrb": query_id,
        "candidate_bmrb": hit.entry_id,
        "percent_id": _fmt_num(hit.percent_id, 2),
        "query_coverage": _fmt_num(cov, 4),
    }

    cand_cached = ensure_candidate_star(hit.entry_id, cs_dir, do_fetch)
    if cand_cached is None:
        return None, {**audit_base, "status": "reject", "reason": "star_missing"}

    n_poly = count_polypeptide_polymers(cand_cached)
    if n_poly != 1:
        return None, {
            **audit_base,
            "status": "reject",
            "reason": f"n_polypeptide={n_poly}",
        }

    match_prov = extract_provenance(cand_cached)
    if not match_prov.adequate:
        return None, {**audit_base, "status": "reject", "reason": "no_pubmed_or_doi"}

    match_ec = extract_conditions_from_star(str(cand_cached))
    ok, detail = conditions_match(
        query_ec, match_ec, ph_tol=ph_tol, temp_tol=temp_tol, ionic_tol=ionic_tol
    )
    if not ok:
        reasons = []
        if detail["missing_ph"]:
            reasons.append("missing_ph")
        elif not detail["ph_ok"]:
            reasons.append("ph_mismatch")
        if detail["missing_temp"]:
            reasons.append("missing_temp")
        elif not detail["temp_ok"]:
            reasons.append("temp_mismatch")
        if detail["missing_ionic"]:
            reasons.append("missing_ionic")
        elif not detail["ionic_ok"]:
            reasons.append("ionic_mismatch")
        if detail["missing_field"]:
            reasons.append("missing_field")
        elif not detail["spectrometer_ok"]:
            reasons.append("field_mismatch")
        return None, {
            **audit_base,
            "status": "reject",
            "reason": ",".join(reasons) or "conditions_fail",
        }

    # Accepted: copy query + match into shifts/ and optionally export parsed tables
    query_copy = copy_star_to_shifts(query_star, shifts_dir)
    match_star = copy_star_to_shifts(cand_cached, shifts_dir)
    match_buffer = extract_buffer_text(match_star)
    if export_parsed:
        export_parsed_shifts(match_star)

    pair = build_pair_row(
        query_id=query_id,
        hit=hit,
        query_len=query_len,
        query_star=query_copy,
        match_star=match_star,
        query_prov=query_prov,
        match_prov=match_prov,
        query_ec=query_ec,
        match_ec=match_ec,
        query_buffer=query_buffer,
        match_buffer=match_buffer,
        cond_detail=detail,
        holo_ids=holo_ids,
    )
    return pair, {**audit_base, "status": "accept", "reason": "ok"}


def process_query_apo(
    apo_id: str,
    *,
    cs_dir: Path,
    shifts_dir: Path,
    do_fetch: bool,
    export_parsed: bool,
    min_identity: float,
    min_coverage: float,
    ph_tol: float,
    temp_tol: float,
    ionic_tol: float,
    holo_ids: Set[str],
    max_candidates: Optional[int],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    pairs: List[Dict[str, Any]] = []
    audits: List[Dict[str, Any]] = []

    query_star = ensure_star_path(apo_id, cs_dir, do_fetch)
    if query_star is None:
        audits.append(
            {
                "query_apo_bmrb": apo_id,
                "candidate_bmrb": "",
                "percent_id": "",
                "query_coverage": "",
                "status": "reject",
                "reason": "query_star_missing",
            }
        )
        return pairs, audits

    seq = extract_query_sequence(query_star)
    if not seq:
        audits.append(
            {
                "query_apo_bmrb": apo_id,
                "candidate_bmrb": "",
                "percent_id": "",
                "query_coverage": "",
                "status": "reject",
                "reason": "query_sequence_empty",
            }
        )
        return pairs, audits

    query_ec = extract_conditions_from_star(str(query_star))
    query_prov = extract_provenance(query_star)
    query_buffer = extract_buffer_text(query_star)

    hits = search_bmrb_fasta(seq)
    filtered = filter_fasta_hits(
        hits,
        query_id=apo_id,
        query_len=len(seq),
        min_identity=min_identity,
        min_coverage=min_coverage,
    )
    if max_candidates is not None:
        filtered = filtered[:max_candidates]

    print(
        f"[APO] {apo_id}: seq_len={len(seq)} fasta_hits={len(hits)} "
        f"after_filter={len(filtered)}",
        flush=True,
    )

    if not filtered:
        audits.append(
            {
                "query_apo_bmrb": apo_id,
                "candidate_bmrb": "",
                "percent_id": "",
                "query_coverage": "",
                "status": "reject",
                "reason": "no_fasta_hits_passing_sequence_filter",
            }
        )
        return pairs, audits

    for hit in filtered:
        pair, audit = evaluate_candidate(
            query_id=apo_id,
            hit=hit,
            query_len=len(seq),
            query_star=query_star,
            query_ec=query_ec,
            query_prov=query_prov,
            query_buffer=query_buffer,
            cs_dir=cs_dir,
            shifts_dir=shifts_dir,
            do_fetch=do_fetch,
            export_parsed=export_parsed,
            ph_tol=ph_tol,
            temp_tol=temp_tol,
            ionic_tol=ionic_tol,
            holo_ids=holo_ids,
        )
        audits.append(audit)
        if pair is not None:
            pairs.append(pair)

    return pairs, audits


def write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fieldnames), extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _resolve_star_for_row(
    bmrb_id: str,
    csv_rel: str,
    *,
    pairs_csv: Path,
    cs_dir: Path,
) -> Optional[Path]:
    """Resolve a STAR path from an existing pairs-CSV row."""
    raw = (csv_rel or "").strip()
    candidates: List[Path] = []
    if raw:
        p = Path(raw)
        candidates.append(p if p.is_absolute() else _REPO_ROOT / p)
        candidates.append(pairs_csv.parent / p)
        candidates.append(pairs_csv.parent / p.name)
        candidates.append(pairs_csv.parent / "shifts" / p.name)
    for suffix in ("_21.str", "_3.str"):
        candidates.append(pairs_csv.parent / "shifts" / f"{bmrb_id}{suffix}")
        candidates.append(cs_dir / f"{bmrb_id}{suffix}")
    seen = set()
    for p in candidates:
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        if p.is_file() and p.stat().st_size > 0:
            return p
    return resolve_star_path(bmrb_id, cs_dir)


def enrich_pairs_csv(
    pairs_csv: Path,
    *,
    cs_dir: Path,
    out_csv: Optional[Path] = None,
) -> Path:
    """Re-extract pH, temperature, buffer, and entities onto an existing pairs CSV."""
    with open(pairs_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        dest = out_csv or pairs_csv
        write_csv(dest, PAIR_FIELDS, [])
        return dest

    cache: Dict[str, Tuple[str, str, ExtractedConditions]] = {}

    def _lookup(bmrb_id: str, csv_rel: str) -> Tuple[str, str, ExtractedConditions]:
        star = _resolve_star_for_row(bmrb_id, csv_rel, pairs_csv=pairs_csv, cs_dir=cs_dir)
        key = str(star) if star else f"missing:{bmrb_id}"
        if key not in cache:
            if star is None:
                cache[key] = ("", "", ExtractedConditions())
            else:
                cache[key] = (
                    extract_buffer_text(star),
                    extract_entity_names(star),
                    extract_conditions_from_star(str(star)),
                )
        return cache[key]

    enriched: List[Dict[str, Any]] = []
    for row in rows:
        qid = (row.get("query_apo_bmrb") or "").strip()
        mid = (row.get("match_apo_bmrb") or "").strip()
        q_buf, q_ent, q_ec = _lookup(qid, row.get("query_star_path") or "")
        m_buf, m_ent, m_ec = _lookup(mid, row.get("match_star_path") or "")
        qv = row_to_csv_values(q_ec)
        mv = row_to_csv_values(m_ec)
        out = dict(row)
        out["query_buffer"] = q_buf
        out["match_buffer"] = m_buf
        out["query_entities"] = q_ent
        out["match_entities"] = m_ent
        if qv["pH"]:
            out["query_pH"] = qv["pH"]
        if mv["pH"]:
            out["match_pH"] = mv["pH"]
        if qv["temperature_C"]:
            out["query_temperature_C"] = qv["temperature_C"]
        if mv["temperature_C"]:
            out["match_temperature_C"] = mv["temperature_C"]
        if qv["ionic_strength_mM"]:
            out["query_ionic_strength_mM"] = qv["ionic_strength_mM"]
        if mv["ionic_strength_mM"]:
            out["match_ionic_strength_mM"] = mv["ionic_strength_mM"]
        enriched.append(out)

    dest = out_csv or pairs_csv
    existing = list(rows[0].keys())
    fields: List[str] = []
    for col in PAIR_FIELDS:
        if col not in fields:
            fields.append(col)
    for col in existing:
        if col not in fields:
            fields.append(col)
    write_csv(dest, fields, enriched)
    return dest


def main(argv: Optional[Sequence[str]] = None) -> int:
    cfg = Paths()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input",
        type=Path,
        default=Path(cfg.input_csv),
        help="Input CSV with apo_bmrb column (default: data/CSP_UBQ.csv)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(cfg.outputs_dir) / "apo_apo_matches",
        help="Output directory (default: outputs/apo_apo_matches)",
    )
    ap.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(cfg.cs_cache_dir),
        help="STAR cache directory (default: CS_Lists)",
    )
    ap.add_argument(
        "--fetch",
        action="store_true",
        help="Download missing STAR files via bmrb_io.fetch_bmrb",
    )
    ap.add_argument(
        "--ids",
        type=str,
        default="",
        help="Comma-separated apo_bmrb IDs to process (default: all unique)",
    )
    ap.add_argument("--min-identity", type=float, default=95.0, help="Min FASTA %% identity")
    ap.add_argument(
        "--min-coverage",
        type=float,
        default=0.9,
        help="Min query coverage (aligned query span / query length)",
    )
    ap.add_argument("--ph-tol", type=float, default=PH_TOLERANCE)
    ap.add_argument("--temp-tol", type=float, default=TEMP_TOLERANCE_C)
    ap.add_argument("--ionic-tol", type=float, default=IONIC_TOLERANCE_MM)
    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for per-query processing (default: 1)",
    )
    ap.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="Cap FASTA candidates evaluated per query (after sequence filter)",
    )
    ap.add_argument(
        "--no-export-parsed",
        action="store_true",
        help="Skip residue-shift CSV export for accepted match STARs",
    )
    ap.add_argument(
        "--enrich",
        type=Path,
        default=None,
        help=(
            "Re-extract pH, temperature, buffer, and entities onto an existing "
            "apo_apo_pairs.csv and write it back (skips FASTA search)."
        ),
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    if args.enrich is not None:
        dest = enrich_pairs_csv(args.enrich, cs_dir=args.cs_dir)
        print(f"Enriched {args.enrich} → {dest}")
        return 0

    id_filter: Optional[Set[str]] = None
    if args.ids.strip():
        id_filter = {x.strip() for x in args.ids.split(",") if x.strip()}

    apo_ids, holo_ids = load_unique_apo_ids(args.input, id_filter)
    if not apo_ids:
        print("No apo_bmrb IDs to process.", file=sys.stderr)
        return 1

    out_dir: Path = args.out
    shifts_dir = out_dir / "shifts"
    shifts_dir.mkdir(parents=True, exist_ok=True)
    # Ensure parsed/ exists under shifts for any exports that land beside STAR
    (shifts_dir / "parsed").mkdir(parents=True, exist_ok=True)

    all_pairs: List[Dict[str, Any]] = []
    all_audits: List[Dict[str, Any]] = []
    export_parsed = not args.no_export_parsed

    def _run_one(apo_id: str) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        return process_query_apo(
            apo_id,
            cs_dir=args.cs_dir,
            shifts_dir=shifts_dir,
            do_fetch=args.fetch,
            export_parsed=export_parsed,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            ph_tol=args.ph_tol,
            temp_tol=args.temp_tol,
            ionic_tol=args.ionic_tol,
            holo_ids=holo_ids,
            max_candidates=args.max_candidates,
        )

    if args.workers <= 1:
        for apo_id in apo_ids:
            pairs, audits = _run_one(apo_id)
            all_pairs.extend(pairs)
            all_audits.extend(audits)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_run_one, apo_id): apo_id for apo_id in apo_ids}
            for fut in as_completed(futs):
                apo_id = futs[fut]
                try:
                    pairs, audits = fut.result()
                except Exception as exc:
                    print(f"[ERROR] {apo_id}: {exc}", file=sys.stderr)
                    all_audits.append(
                        {
                            "query_apo_bmrb": apo_id,
                            "candidate_bmrb": "",
                            "percent_id": "",
                            "query_coverage": "",
                            "status": "reject",
                            "reason": f"exception:{exc}",
                        }
                    )
                    continue
                all_pairs.extend(pairs)
                all_audits.extend(audits)

    # Stable sort
    all_pairs.sort(
        key=lambda r: (r.get("query_apo_bmrb", ""), r.get("match_apo_bmrb", ""))
    )
    all_audits.sort(
        key=lambda r: (
            r.get("query_apo_bmrb", ""),
            r.get("candidate_bmrb", ""),
            r.get("status", ""),
        )
    )

    pairs_path = out_dir / "apo_apo_pairs.csv"
    audit_path = out_dir / "apo_apo_candidates_audit.csv"
    write_csv(pairs_path, PAIR_FIELDS, all_pairs)
    write_csv(audit_path, AUDIT_FIELDS, all_audits)

    n_accept = sum(1 for a in all_audits if a.get("status") == "accept")
    print(
        f"Wrote {len(all_pairs)} accepted pairs to {pairs_path} "
        f"({n_accept} accepts in audit of {len(all_audits)} rows) → {audit_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
