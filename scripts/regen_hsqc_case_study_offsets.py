#!/usr/bin/env python3
"""Regenerate hsqc_scatter.png and case-study v1 figures with applied-offset labels.

Reads existing per-target csp_table / binding tables; does not recompute CSPs.
Reuses cached PyMOL views (skips targets with no saved view).
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.case_study import (  # noqa: E402
    _case_study_view_load_candidates,
    generate_case_study_figure,
)
from scripts.config import paths  # noqa: E402
from scripts.csp import CSPResult  # noqa: E402
from scripts.HSQC_visualize import plot_hsqc_variants  # noqa: E402
from scripts.merge_csv import parse_optional_bool  # noqa: E402


def _f(val: Any) -> Optional[float]:
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _b(val: Any) -> Optional[bool]:
    parsed = parse_optional_bool(val)
    if parsed is None and val is not None and str(val).strip() != "":
        # parse_optional_bool returns None for empty; keep False/True for 0/1
        return None
    return parsed


def load_csp_results(csp_table_path: Path) -> List[CSPResult]:
    results: List[CSPResult] = []
    with open(csp_table_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            apo_resi = _f(row.get("apo_resi"))
            holo_resi = _f(row.get("holo_resi"))
            if apo_resi is None or holo_resi is None:
                continue
            results.append(
                CSPResult(
                    apo_index=int(apo_resi),
                    holo_index=int(holo_resi),
                    apo_aa=(row.get("apo_aa") or "").strip() or "X",
                    holo_aa=(row.get("holo_aa") or "").strip() or "X",
                    H_apo=_f(row.get("H_apo")),
                    N_apo=_f(row.get("N_apo")),
                    H_holo=_f(row.get("H_holo")),
                    N_holo=_f(row.get("N_holo")),
                    dH=_f(row.get("dH")),
                    dN=_f(row.get("dN")),
                    csp_A=_f(row.get("csp_A")),
                    significant=_b(row.get("significant")),
                    H_holo_original=_f(row.get("H_holo_original")),
                    N_holo_original=_f(row.get("N_holo_original")),
                    H_offset=_f(row.get("H_offset")),
                    N_offset=_f(row.get("N_offset")),
                    z_score=_f(row.get("csp_z")),
                    CA_apo=_f(row.get("CA_apo")),
                    CA_holo=_f(row.get("CA_holo")),
                    HA_apo=_f(row.get("HA_apo")),
                    HA_holo=_f(row.get("HA_holo")),
                    CA_offset=_f(row.get("CA_offset")),
                    HA_offset=_f(row.get("HA_offset")),
                )
            )
    return results


def load_binding_results(target_dir: Path) -> Optional[dict]:
    """Build a union-style binding_results dict from master_alignment when present."""
    master_path = target_dir / "master_alignment.csv"
    if not master_path.exists():
        return None
    residue_info: List[Dict[str, Any]] = []
    with open(master_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pdb_res = _f(row.get("pdb_residue_number"))
            if pdb_res is None:
                continue
            is_binding = bool(_b(row.get("is_binding_site")))
            aa = (row.get("holo_aa") or row.get("apo_aa") or "").strip()
            residue_info.append(
                {
                    "residue_number": int(pdb_res),
                    "residue_name": aa,
                    "has_hbond": is_binding,
                    "has_charge_complement": False,
                    "has_pi_contact": False,
                    "has_sasa_occlusion": False,
                    "has_ca_distance": False,
                    "has_any_atom_sub_2A": False,
                }
            )
    if not residue_info:
        return None
    n_union = sum(1 for info in residue_info if info["has_hbond"])
    return {
        "residue_info": residue_info,
        "dataset_type": "union",
        "n_union_residues": n_union,
    }


def load_metadata_index(csv_path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """Map (holo_pdb_upper, apo_bmrb) -> row fields."""
    index: Dict[Tuple[str, str], Dict[str, str]] = {}
    if not csv_path.exists():
        return index
    with open(csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            holo = (row.get("holo_pdb") or row.get("holo_pdb_id") or "").strip().upper()
            apo = (row.get("apo_bmrb") or "").strip()
            if holo and apo:
                index[(holo, apo)] = row
    return index


def has_saved_view(pdb_id: str, view_key: str) -> bool:
    for cand in _case_study_view_load_candidates(paths.pymol_views_dir, pdb_id, view_key):
        if os.path.exists(cand):
            return True
    return False


def process_target(
    target_dir: Path,
    meta_index: Dict[Tuple[str, str], Dict[str, str]],
    *,
    skip_case_study: bool = False,
) -> Tuple[str, Optional[Tuple[float, float]]]:
    csp_path = target_dir / "csp_table.csv"
    results = load_csp_results(csp_path)
    if not results:
        raise RuntimeError("no CSP rows")

    parts = target_dir.name.split("_", 1)
    pdb_id = parts[0]
    apo_bmrb = parts[1] if len(parts) > 1 else ""
    meta = meta_index.get((pdb_id.upper(), apo_bmrb), {})
    holo_bmrb = (meta.get("holo_bmrb") or "").strip() or None
    apo_pdb = (meta.get("apo_pdb") or "").strip() or None
    if holo_bmrb is None:
        # Fall back to first csp_table row
        with open(csp_path, newline="", encoding="utf-8") as handle:
            row = next(csv.DictReader(handle), None)
            if row:
                holo_bmrb = (row.get("holo_bmrb") or "").strip() or None

    binding = load_binding_results(target_dir)
    hsqc_path = target_dir / "hsqc_scatter.png"
    offsets = plot_hsqc_variants(
        results,
        str(hsqc_path),
        title=f"{pdb_id} HSQC comparison",
        binding_results=binding,
    )

    if not skip_case_study:
        if not has_saved_view(pdb_id, target_dir.name):
            raise RuntimeError("no saved PyMOL view")
        generate_case_study_figure(
            str(target_dir),
            pdb_id,
            apo_bmrb=apo_bmrb or None,
            holo_bmrb=holo_bmrb,
            apo_pdb=apo_pdb,
            force_view_reset=False,
            view_key=target_dir.name,
        )
    return target_dir.name, offsets


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=PROJECT_ROOT / "outputs",
        help="Outputs root containing per-target folders",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / "data" / "CSP_UBQ.csv",
        help="CSV used for apo_pdb / holo_bmrb metadata headers",
    )
    parser.add_argument(
        "--ids",
        type=str,
        default="",
        help="Optional comma-separated target folder names (e.g. 1CF4_18251)",
    )
    parser.add_argument(
        "--hsqc-only",
        action="store_true",
        help="Only regenerate hsqc_scatter.png (skip case-study PyMOL compose)",
    )
    args = parser.parse_args()

    meta_index = load_metadata_index(args.input)
    if args.ids.strip():
        wanted = {x.strip() for x in args.ids.split(",") if x.strip()}
        targets = sorted(
            p
            for p in args.out.iterdir()
            if p.is_dir() and p.name in wanted and (p / "csp_table.csv").exists()
        )
    else:
        targets = sorted(
            p
            for p in args.out.iterdir()
            if p.is_dir() and (p / "csp_table.csv").exists() and (p / "hsqc_scatter.png").exists()
        )

    ok = 0
    fail = 0
    print(f"Regenerating {len(targets)} targets under {args.out}")
    for i, target_dir in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] {target_dir.name}", flush=True)
        try:
            name, offsets = process_target(
                target_dir,
                meta_index,
                skip_case_study=args.hsqc_only,
            )
            print(f"  OK {name} offsets={offsets}", flush=True)
            ok += 1
        except Exception as exc:
            print(f"  FAIL {target_dir.name}: {exc}", flush=True)
            fail += 1

    print("=" * 60)
    print(f"DONE ok={ok} fail={fail} targets={len(targets)}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
