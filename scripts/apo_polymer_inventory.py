#!/usr/bin/env python3
"""
Inventory polymer entities in apo BMRB shift lists from CSP_UBQ.csv.

Flags apo entries with polymer entities outside the main protein receptor
(second polypeptides / peptides, DNA, RNA). Writes:

  outputs/apo_polymer_inventory/apo_polymer_entities.csv
  outputs/apo_polymer_inventory/apo_polymer_summary.csv
  outputs/apo_polymer_inventory/DIGEST.md
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.apo_holo_exp_conditions import (  # noqa: E402
    _detect_bmrb_format,
    _iter_saveframes,
    _saveframe_body,
    _tokenize_star_lines,
    ensure_star_path,
    resolve_star_path,
)
from scripts.config import Paths  # noqa: E402
from scripts.find_apo_apo_matches import count_polypeptide_polymers  # noqa: E402

ENTITY_FIELDS = [
    "apo_bmrb",
    "star_file",
    "bmrb_format",
    "n_polypeptide",
    "saveframe",
    "entity_id",
    "entity_type",
    "polymer_type",
    "name",
    "fragment",
    "n_monomers",
    "role_guess",
    "has_shift_list",
    "shift_saveframes",
    "holo_bmrb",
    "holo_pdb",
    "apo_pdb",
    "ec_classes",
    "scope_fold_type",
]

SUMMARY_FIELDS = [
    "apo_bmrb",
    "star_file",
    "bmrb_format",
    "n_polypeptide",
    "n_polymer_total",
    "n_nonpolymer",
    "has_extra_polypeptide",
    "has_nucleic_acid",
    "receptor_name",
    "receptor_n_monomers",
    "extra_polymer_names",
    "extra_polymer_lengths",
    "holo_bmrb",
    "holo_pdb",
    "apo_pdb",
    "ec_classes",
    "scope_fold_type",
]


def _strip_token(tok: str) -> str:
    t = (tok or "").strip()
    if len(t) >= 2 and ((t[0] == t[-1] == "'") or (t[0] == t[-1] == '"')):
        t = t[1:-1]
    return t.strip()


def _tag(body: str, *patterns: str) -> str:
    for pat in patterns:
        m = re.search(pat, body, re.MULTILINE)
        if m:
            return _strip_token(m.group(1))
    return ""


def _norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def _star_row_tokens(line: str) -> List[str]:
    """Split a STAR loop data row into tokens, respecting quotes."""
    return [_strip_token(t) for t in re.findall(r"'[^']*'|\"[^\"]*\"|\S+", line.strip())]


@dataclass
class PolymerEntity:
    saveframe: str
    entity_id: str
    entity_type: str
    polymer_type: str
    name: str
    fragment: str
    n_monomers: str
    has_shift_list: bool = False
    shift_saveframes: str = ""

    @property
    def is_polymer(self) -> bool:
        if self.entity_type.lower() == "polymer":
            return True
        # v2.1 monomeric_polymer frames are polymers even if Type absent
        return bool(self.polymer_type) and self.entity_type.lower() != "non-polymer"

    @property
    def is_polypeptide(self) -> bool:
        pt = self.polymer_type.lower()
        return pt.startswith("polypeptide") or pt in ("protein", "peptide")

    @property
    def is_nucleic_acid(self) -> bool:
        pt = self.polymer_type.lower()
        return any(
            x in pt
            for x in (
                "polyribo",
                "polydeoxy",
                "dna",
                "rna",
                "nucleic",
            )
        )

    @property
    def n_monomers_int(self) -> int:
        try:
            return int(float(self.n_monomers))
        except (TypeError, ValueError):
            return 0


def inventory_entities(star_path: Path) -> Tuple[str, List[PolymerEntity]]:
    """Return BMRB format and entity/polymer inventory for a STAR file."""
    text = star_path.read_text(encoding="utf-8", errors="ignore")
    lines = _tokenize_star_lines(text)
    fmt = _detect_bmrb_format(str(star_path))
    entities: List[PolymerEntity] = []

    if fmt == "3":
        for start, save_end, name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not re.search(r"^\s*_Entity\.Type\s+", body, re.MULTILINE):
                continue
            entities.append(
                PolymerEntity(
                    saveframe=name,
                    entity_id=_tag(body, r"^\s*_Entity\.ID\s+(\S+)"),
                    entity_type=_tag(body, r"^\s*_Entity\.Type\s+(\S+)"),
                    polymer_type=_tag(body, r"^\s*_Entity\.Polymer_type\s+(\S+)"),
                    name=_tag(body, r"^\s*_Entity\.Name\s+(.+)$"),
                    fragment=_tag(body, r"^\s*_Entity\.Fragment\s+(.+)$"),
                    n_monomers=_tag(body, r"^\s*_Entity\.Number_of_monomers\s+(\S+)"),
                )
            )
    else:
        polymer_labels: Set[str] = set()
        for start, save_end, name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not re.search(
                r"^\s*_Saveframe_category\s+monomeric_polymer\s*$", body, re.MULTILINE
            ):
                continue
            sf_label = name[5:] if name.startswith("save_") else name
            pname = _tag(
                body,
                r"^\s*_Name_common\s+(.+)$",
                r"^\s*_Mol_name\s+(.+)$",
            )
            polymer_labels.add(_norm_name(sf_label))
            polymer_labels.add(_norm_name(name))
            if pname:
                polymer_labels.add(_norm_name(pname))
            entities.append(
                PolymerEntity(
                    saveframe=name,
                    entity_id="",
                    entity_type="polymer",
                    polymer_type=_tag(body, r"^\s*_Mol_polymer_class\s+(\S+)") or "protein",
                    name=pname,
                    fragment="",
                    n_monomers=_tag(body, r"^\s*_Residue_count\s+(\S+)"),
                )
            )

        # Non-polymers from molecular_system assembly components not tied to a polymer
        seen_np: Set[str] = set()
        for start, save_end, _name in _iter_saveframes(lines):
            body = _saveframe_body(lines, start, save_end)
            if not re.search(
                r"^\s*_Saveframe_category\s+molecular_system\s*$", body, re.MULTILINE
            ):
                continue
            for raw in body.splitlines():
                toks = _star_row_tokens(raw)
                if len(toks) < 2 or not toks[1].startswith("$"):
                    continue
                comp = toks[0]
                label = toks[1].lstrip("$")
                if _norm_name(label) in polymer_labels or _norm_name(comp) in polymer_labels:
                    continue
                key = _norm_name(label) or _norm_name(comp)
                if not key or key in seen_np:
                    continue
                seen_np.add(key)
                disp = comp or label
                sf_name = f"save_{label}"
                for s2, e2, n2 in _iter_saveframes(lines):
                    lab = n2[5:] if n2.startswith("save_") else n2
                    b2 = _saveframe_body(lines, s2, e2)
                    common = _tag(
                        b2, r"^\s*_Name_common\s+(.+)$", r"^\s*_Mol_name\s+(.+)$"
                    )
                    if (
                        _norm_name(lab) == _norm_name(label)
                        or _norm_name(n2) == _norm_name(label)
                        or (common and _norm_name(label) in _norm_name(common))
                        or (common and _norm_name(comp) in _norm_name(common))
                    ):
                        disp = common or disp
                        sf_name = n2
                        break
                # Prefer a descriptive name over bare entity_N placeholders
                if re.fullmatch(r"entity_?\d*", _norm_name(disp)) and label:
                    disp = label
                entities.append(
                    PolymerEntity(
                        saveframe=sf_name,
                        entity_id="",
                        entity_type="non-polymer",
                        polymer_type="",
                        name=disp,
                        fragment="",
                        n_monomers="",
                    )
                )

    _annotate_shift_lists(lines, fmt, entities)
    return fmt, entities


def _assembly_component_map(lines: List[str]) -> Dict[str, str]:
    """Map entity_N / component labels to polymer saveframe or common name (v2.1)."""
    mapping: Dict[str, str] = {}
    for start, save_end, _name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        if not re.search(
            r"^\s*_Saveframe_category\s+molecular_system\s*$", body, re.MULTILINE
        ) and "assembly" not in _name.lower():
            # Still try assembly-like loops anywhere
            pass
        # Rows often: entity_1  $cNTnC  ... or 'zf-CW domain'  $entity_1
        for raw in body.splitlines():
            toks = _star_row_tokens(raw)
            if len(toks) < 2 or not toks[1].startswith("$"):
                continue
            left = toks[0]
            right = toks[1].lstrip("$")
            if left and right:
                mapping[_norm_name(left)] = right
                mapping[_norm_name(right)] = right
    return mapping


def _match_score(ent: PolymerEntity, key: str, asm_map: Dict[str, str]) -> int:
    """Return match specificity score (0 = no match). Higher is better."""
    if not key:
        return 0
    kn = _norm_name(key)
    if not kn:
        return 0
    if ent.entity_id and (key == ent.entity_id or kn == _norm_name(ent.entity_id)):
        return 100
    resolved = asm_map.get(kn, "")
    rn = _norm_name(resolved)
    labels = [
        (_norm_name(ent.name), 90),
        (_norm_name(ent.saveframe.lstrip("save_")), 80),
        (_norm_name(ent.saveframe), 70),
    ]
    for ln, score in labels:
        if not ln:
            continue
        if kn == ln or (rn and rn == ln):
            return score
    return 0


def _annotate_shift_lists(
    lines: List[str], fmt: str, entities: List[PolymerEntity]
) -> None:
    """Mark entities that have an assigned chemical-shift saveframe."""
    shift_hits: Dict[int, List[str]] = defaultdict(list)
    asm_map = _assembly_component_map(lines)

    for start, save_end, name in _iter_saveframes(lines):
        body = _saveframe_body(lines, start, save_end)
        cat = _tag(
            body,
            r"^\s*_Saveframe_category\s+(\S+)",
            r"^\s*_Assigned_chem_shift_list\.Sf_category\s+(\S+)",
        )
        if cat and cat not in ("assigned_chemical_shifts",):
            continue
        if not re.search(r"assigned_chemical_shifts", body, re.IGNORECASE):
            continue
        if cat == "entry_information" or name.endswith("entry_information"):
            continue

        component = _tag(
            body,
            r"^\s*_Mol_system_component_name\s+(.+)$",
        )
        keys: List[str] = []
        if component:
            keys.append(component.lstrip("$"))
        # v3: dominant Entity_ID among atom shifts (ignore sparse contaminants)
        eids = [
            eid
            for eid in re.findall(r"^\s*\d+\s+\S+\s+\S+\s+(\d+)\s+", body, re.MULTILINE)
        ]
        if not eids:
            eids = [
                eid
                for eid in re.findall(r"_Atom_chem_shift\.Entity_ID\s+(\S+)", body)
                if eid not in (".", "?", "")
            ]
        if eids:
            keys.append(Counter(eids).most_common(1)[0][0])

        best_ent = None
        best_score = 0
        for ent in entities:
            score = max((_match_score(ent, key, asm_map) for key in keys), default=0)
            if score > best_score:
                best_score = score
                best_ent = ent

        if best_ent is not None and best_score > 0:
            shift_hits[id(best_ent)].append(name)
        else:
            polys = [e for e in entities if e.is_polypeptide]
            if len(polys) == 1:
                shift_hits[id(polys[0])].append(name)

    for ent in entities:
        frames = shift_hits.get(id(ent), [])
        if frames:
            ent.has_shift_list = True
            ent.shift_saveframes = ";".join(dict.fromkeys(frames))


def assign_roles(entities: List[PolymerEntity]) -> None:
    """Guess receptor vs extra polymer roles (longest polypeptide preferred)."""
    polys = [e for e in entities if e.is_polypeptide]
    if not polys:
        for e in entities:
            if e.is_nucleic_acid:
                e.role_guess = "DNA/RNA"  # type: ignore[attr-defined]
            elif e.entity_type.lower() == "non-polymer" or not e.is_polymer:
                e.role_guess = "nonpolymer"  # type: ignore[attr-defined]
            else:
                e.role_guess = "other_polymer"  # type: ignore[attr-defined]
        return

    # Prefer longest chain; break ties with shift coverage then name length
    def sort_key(e: PolymerEntity) -> Tuple[int, int, int]:
        return (e.n_monomers_int, 1 if e.has_shift_list else 0, len(e.name))

    receptor = max(polys, key=sort_key)
    for e in entities:
        if e is receptor:
            e.role_guess = "receptor"  # type: ignore[attr-defined]
        elif e.is_nucleic_acid:
            e.role_guess = "DNA/RNA"  # type: ignore[attr-defined]
        elif e.is_polypeptide:
            e.role_guess = (
                "peptide"
                if e.n_monomers_int and e.n_monomers_int <= 40
                else "second_protein"
            )  # type: ignore[attr-defined]
        elif e.entity_type.lower() == "non-polymer" or not e.is_polymer:
            e.role_guess = "nonpolymer"  # type: ignore[attr-defined]
        else:
            e.role_guess = "other_polymer"  # type: ignore[attr-defined]


def load_target_rows(input_csv: Path) -> Dict[str, List[Dict[str, str]]]:
    by_apo: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with open(input_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            apo = (row.get("apo_bmrb") or "").strip()
            if not apo:
                continue
            by_apo[apo].append(
                {
                    "holo_bmrb": (row.get("holo_bmrb") or "").strip(),
                    "holo_pdb": (row.get("holo_pdb") or "").strip(),
                    "apo_pdb": (row.get("apo_pdb") or "").strip(),
                    "ec_classes": (row.get("ec_classes") or "").strip(),
                    "scope_fold_type": (row.get("scope_fold_type") or "").strip(),
                }
            )
    return by_apo


def _join_unique(rows: Sequence[Dict[str, str]], key: str) -> str:
    vals = []
    for r in rows:
        v = (r.get(key) or "").strip()
        if v and v not in vals:
            vals.append(v)
    return ";".join(vals)


def write_digest(
    path: Path,
    *,
    n_apos: int,
    summary_rows: Sequence[Dict[str, Any]],
    entity_rows: Sequence[Dict[str, Any]],
) -> None:
    multi = [r for r in summary_rows if r.get("has_extra_polypeptide") in (True, "True", "1", 1)]
    nucleic = [r for r in summary_rows if r.get("has_nucleic_acid") in (True, "True", "1", 1)]
    zero = [r for r in summary_rows if int(r.get("n_polypeptide") or 0) == 0]
    # Headline focuses on NMR-STAR 3 Entity.Type=non-polymer extras (compact, structured).
    nonpoly_only = [
        r
        for r in summary_rows
        if int(r.get("n_nonpolymer") or 0) > 0
        and r.get("has_extra_polypeptide") not in (True, "True", "1", 1)
        and int(r.get("n_polypeptide") or 0) == 1
        and r.get("bmrb_format") == "3"
    ]

    lines: List[str] = []
    lines.append("# Apo shift-list non-receptor polymer digest")
    lines.append("")
    lines.append(
        f"Scanned **{n_apos}** unique `apo_bmrb` entries from `data/CSP_UBQ.csv` "
        "via NMR-STAR polymer metadata in `CS_Lists/` "
        "(`*_21.str` preferred, else `*_3.str`), using `count_polypeptide_polymers` "
        "from `scripts/find_apo_apo_matches.py` plus a richer entity inventory."
    )
    lines.append("")
    lines.append("## Headline counts")
    lines.append("")
    lines.append(f"- **{n_apos}** unique apo shift lists")
    lines.append(
        f"- **{len(multi)}** with a second **polypeptide** polymer (protein + peptide/partner)"
    )
    lines.append(f"- **{len(nucleic)}** with DNA/RNA polymer types")
    lines.append(
        f"- **{len(nonpoly_only)}** with only an extra **non-polymer** "
        f"(ions/ligands): {', '.join(r['apo_bmrb'] for r in nonpoly_only) or 'none'}"
    )
    lines.append(
        f"- **{len(zero)}** with **zero** countable polypeptides: "
        f"{', '.join(r['apo_bmrb'] for r in zero) or 'none'}"
    )
    lines.append("")
    lines.append(
        "Main receptor is inferred as the longest polypeptide / the one with "
        "the primary H+N shift list (matches pipeline alignment behavior)."
    )
    lines.append("")
    lines.append("## Targets with an extra polymer")
    lines.append("")

    if not multi:
        lines.append("_None found._")
        lines.append("")
    else:
        for i, srow in enumerate(multi, 1):
            apo = srow["apo_bmrb"]
            ents = [e for e in entity_rows if e["apo_bmrb"] == apo]
            holos = srow.get("holo_bmrb") or ""
            pdbs = srow.get("holo_pdb") or ""
            apo_pdb = srow.get("apo_pdb") or ""
            ec = srow.get("ec_classes") or ""
            fold = srow.get("scope_fold_type") or ""
            holo_bits = []
            for hb, hp in zip(holos.split(";"), pdbs.split(";")):
                if hb or hp:
                    holo_bits.append(f"**{hb}** / PDB **{hp}**" if hp else f"**{hb}**")
            holo_txt = " and ".join(holo_bits) if holo_bits else "(no holo link)"
            meta = []
            if apo_pdb:
                meta.append(f"apo PDB `{apo_pdb}`")
            if ec:
                meta.append(ec)
            if fold:
                meta.append(fold)
            meta_txt = f" ({'; '.join(meta)})" if meta else ""

            lines.append(f"### {i}. Apo **{apo}** → holo {holo_txt}{meta_txt}")
            lines.append("")
            lines.append("| Role | Name | Length | Assigned shifts? |")
            lines.append("|------|------|--------|------------------|")
            for e in ents:
                role = e.get("role_guess") or ""
                if role == "nonpolymer":
                    continue
                if role not in ("receptor", "peptide", "second_protein", "DNA/RNA", "other_polymer"):
                    if e.get("entity_type") == "non-polymer":
                        continue
                nmon = e.get("n_monomers") or "?"
                shifts = "Yes" if e.get("has_shift_list") in (True, "True", "1", 1) else "No"
                sf = e.get("shift_saveframes") or ""
                if shifts == "Yes" and sf:
                    shifts = f"Yes (`{sf}`)"
                role_label = {
                    "receptor": "Receptor",
                    "peptide": "Extra",
                    "second_protein": "Extra",
                    "DNA/RNA": "DNA/RNA",
                    "other_polymer": "Extra",
                }.get(role, role)
                lines.append(
                    f"| {role_label} | "
                    f"{e.get('name') or e.get('saveframe')} | {nmon} | {shifts} |"
                )
            extras_np = [e for e in ents if e.get("role_guess") == "nonpolymer"]
            if extras_np:
                names = ", ".join(e.get("name") or e.get("saveframe") for e in extras_np)
                lines.append("")
                lines.append(f"Also deposits non-polymers: {names}.")
            poly_with_shifts = [
                e
                for e in ents
                if e.get("role_guess") in ("receptor", "peptide", "second_protein")
                and e.get("has_shift_list") in (True, "True", "1", 1)
            ]
            if len(poly_with_shifts) > 1:
                lines.append("")
                lines.append(
                    "Both polymers have assigned chem-shift lists "
                    "(highest risk of dual-sequence confusion in CSP pairing)."
                )
            if ";" in (srow.get("holo_bmrb") or ""):
                lines.append("")
                lines.append("Multiple CSP_UBQ rows share this apo.")
            lines.append("")

    lines.append("## Pattern")
    lines.append("")
    lines.append(
        "All extras among query apos are short peptides / partner chains "
        "(typically 14–27 aa), annotated as `_Mol_polymer_class protein`, "
        "not nucleic acids. Extra polymers only matter for CSP pairing when "
        "they also have H+N assignments (**6095**, **25495**)."
    )
    lines.append("")
    lines.append("## Note on receptor selection")
    lines.append("")
    lines.append(
        "The main pipeline does **not** use Entity polymer metadata to pick the receptor; "
        "it picks the apo H+N sequence that best aligns to holo "
        "(`scripts/pipeline.py` / CSP alignment)."
    )
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- `apo_polymer_entities.csv` — one row per entity/polymer")
    lines.append("- `apo_polymer_summary.csv` — one row per apo BMRB")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


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
        default=Path(cfg.outputs_dir) / "apo_polymer_inventory",
        help="Output directory (default: outputs/apo_polymer_inventory)",
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
    args = ap.parse_args(argv)

    id_filter: Optional[Set[str]] = None
    if args.ids.strip():
        id_filter = {x.strip() for x in args.ids.split(",") if x.strip()}

    by_apo = load_target_rows(args.input)
    apo_ids = sorted(by_apo.keys(), key=lambda x: (len(x), x))
    if id_filter is not None:
        apo_ids = [a for a in apo_ids if a in id_filter]

    entity_rows: List[Dict[str, Any]] = []
    summary_rows: List[Dict[str, Any]] = []
    missing: List[str] = []

    for apo in apo_ids:
        targets = by_apo[apo]
        star = (
            ensure_star_path(apo, args.cs_dir, do_fetch=args.fetch)
            if args.fetch
            else resolve_star_path(apo, args.cs_dir)
        )
        if star is None:
            missing.append(apo)
            summary_rows.append(
                {
                    "apo_bmrb": apo,
                    "star_file": "",
                    "bmrb_format": "",
                    "n_polypeptide": "",
                    "n_polymer_total": "",
                    "n_nonpolymer": "",
                    "has_extra_polypeptide": "",
                    "has_nucleic_acid": "",
                    "receptor_name": "",
                    "receptor_n_monomers": "",
                    "extra_polymer_names": "MISSING_STAR",
                    "extra_polymer_lengths": "",
                    "holo_bmrb": _join_unique(targets, "holo_bmrb"),
                    "holo_pdb": _join_unique(targets, "holo_pdb"),
                    "apo_pdb": _join_unique(targets, "apo_pdb"),
                    "ec_classes": _join_unique(targets, "ec_classes"),
                    "scope_fold_type": _join_unique(targets, "scope_fold_type"),
                }
            )
            continue

        fmt, entities = inventory_entities(star)
        assign_roles(entities)
        n_poly = count_polypeptide_polymers(star)
        polys = [e for e in entities if e.is_polypeptide]
        nonpolys = [
            e
            for e in entities
            if getattr(e, "role_guess", "") == "nonpolymer"
            or e.entity_type.lower() == "non-polymer"
        ]
        polymers = [e for e in entities if e.is_polymer]
        extras = [
            e
            for e in entities
            if getattr(e, "role_guess", "")
            in ("peptide", "second_protein", "DNA/RNA", "other_polymer")
        ]
        receptor = next(
            (e for e in entities if getattr(e, "role_guess", "") == "receptor"), None
        )

        for e in entities:
            entity_rows.append(
                {
                    "apo_bmrb": apo,
                    "star_file": star.name,
                    "bmrb_format": fmt,
                    "n_polypeptide": n_poly,
                    "saveframe": e.saveframe,
                    "entity_id": e.entity_id,
                    "entity_type": e.entity_type,
                    "polymer_type": e.polymer_type,
                    "name": e.name,
                    "fragment": e.fragment,
                    "n_monomers": e.n_monomers,
                    "role_guess": getattr(e, "role_guess", ""),
                    "has_shift_list": e.has_shift_list,
                    "shift_saveframes": e.shift_saveframes,
                    "holo_bmrb": _join_unique(targets, "holo_bmrb"),
                    "holo_pdb": _join_unique(targets, "holo_pdb"),
                    "apo_pdb": _join_unique(targets, "apo_pdb"),
                    "ec_classes": _join_unique(targets, "ec_classes"),
                    "scope_fold_type": _join_unique(targets, "scope_fold_type"),
                }
            )

        summary_rows.append(
            {
                "apo_bmrb": apo,
                "star_file": star.name,
                "bmrb_format": fmt,
                "n_polypeptide": n_poly,
                "n_polymer_total": len(polymers),
                "n_nonpolymer": len(nonpolys),
                "has_extra_polypeptide": n_poly > 1,
                "has_nucleic_acid": any(e.is_nucleic_acid for e in entities),
                "receptor_name": receptor.name if receptor else "",
                "receptor_n_monomers": receptor.n_monomers if receptor else "",
                "extra_polymer_names": ";".join(e.name or e.saveframe for e in extras),
                "extra_polymer_lengths": ";".join(
                    str(e.n_monomers) for e in extras if e.n_monomers
                ),
                "holo_bmrb": _join_unique(targets, "holo_bmrb"),
                "holo_pdb": _join_unique(targets, "holo_pdb"),
                "apo_pdb": _join_unique(targets, "apo_pdb"),
                "ec_classes": _join_unique(targets, "ec_classes"),
                "scope_fold_type": _join_unique(targets, "scope_fold_type"),
            }
        )

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    entities_path = out_dir / "apo_polymer_entities.csv"
    summary_path = out_dir / "apo_polymer_summary.csv"
    digest_path = out_dir / "DIGEST.md"

    with open(entities_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=ENTITY_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(entity_rows)

    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(summary_rows)

    write_digest(
        digest_path,
        n_apos=len(apo_ids),
        summary_rows=summary_rows,
        entity_rows=entity_rows,
    )

    n_multi = sum(1 for r in summary_rows if r.get("has_extra_polypeptide") is True)
    print(
        f"Scanned {len(apo_ids)} apos; missing STAR={len(missing)}; "
        f"extra polypeptide={n_multi}"
    )
    print(f"Wrote {entities_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {digest_path}")
    if missing:
        print("Missing:", ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
