#!/usr/bin/env python3
"""
Create Supplementary Figure 4 (`figures/suppl_fig_4.png`) by reproducing the per-atom
classification panel plot used in `outputs/<target>/per_atom_classification_panels.png`.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from csp import CSPResult
from visualize import plot_per_atom_classification_panels


_TRUE_VALUES = {"1", "true", "t", "yes", "y"}


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in _TRUE_VALUES


def _as_int(value: object) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _as_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_csp_results(path: Path) -> List[CSPResult]:
    rows = _read_csv_rows(path)
    results: List[CSPResult] = []
    for row in rows:
        apo_index = _as_int(row.get("apo_resi"))
        holo_index = _as_int(row.get("holo_resi")) or _as_int(row.get("residue_number"))
        if holo_index is None:
            continue
        if apo_index is None:
            apo_index = holo_index

        csp_result = CSPResult(
            apo_index=apo_index,
            holo_index=holo_index,
            apo_aa=(row.get("apo_aa") or "").strip(),
            holo_aa=(row.get("holo_aa") or "").strip(),
            H_apo=_as_float(row.get("H_apo")),
            N_apo=_as_float(row.get("N_apo")),
            H_holo=_as_float(row.get("H_holo")),
            N_holo=_as_float(row.get("N_holo")),
            dH=_as_float(row.get("dH")),
            dN=_as_float(row.get("dN")),
            csp_A=_as_float(row.get("csp_A")) or _as_float(row.get("csp_CA")),
            significant=_as_bool(row.get("significant")) or _as_bool(row.get("csp_CA_significant")),
            significant_1sd=_as_bool(row.get("significant_1sd")) or _as_bool(row.get("csp_CA_significant_1sd")),
            significant_2sd=_as_bool(row.get("significant_2sd")) or _as_bool(row.get("csp_CA_significant_2sd")),
            H_holo_original=_as_float(row.get("H_holo_original")) or _as_float(row.get("H_holo")),
            N_holo_original=_as_float(row.get("N_holo_original")) or _as_float(row.get("N_holo")),
            H_offset=_as_float(row.get("H_offset")),
            N_offset=_as_float(row.get("N_offset")),
            z_score=_as_float(row.get("csp_z")) or _as_float(row.get("csp_CA_z")),
            CA_apo=_as_float(row.get("CA_apo")),
            CA_holo=_as_float(row.get("CA_holo")),
            CA_holo_original=_as_float(row.get("CA_holo_original")) or _as_float(row.get("CA_holo")),
            CA_offset=_as_float(row.get("CA_offset")),
            dCA=_as_float(row.get("dCA")),
        )
        results.append(csp_result)
    return results


def _build_union_binding_results(
    interaction_csv: Path,
    occlusion_csv: Path,
    ca_distance_csv: Path,
) -> Dict[str, object]:
    interaction_rows = _read_csv_rows(interaction_csv)
    occlusion_rows = _read_csv_rows(occlusion_csv)
    ca_distance_rows = _read_csv_rows(ca_distance_csv)

    union_map: Dict[int, Dict[str, object]] = {}

    def _get_or_create_union_entry(
        residue_number: int,
        *,
        residue_name: Optional[str] = None,
        chain_id: Optional[str] = None,
    ) -> Dict[str, object]:
        entry = union_map.get(residue_number)
        if entry is None:
            entry = {
                "residue_number": residue_number,
                "residue_name": residue_name,
                "chain_id": chain_id,
                "has_hbond": False,
                "has_charge_complement": False,
                "has_pi_contact": False,
                "interaction_category": "none",
                "partner_residues": [],
                "hbond_count": 0,
                "charge_pair_count": 0,
                "pi_contact_count": 0,
                "hbond_partner_residues": [],
                "charge_partner_residues": [],
                "pi_partner_residues": [],
                "is_charged": False,
                "charge_type": "neutral",
                "has_sasa_occlusion": False,
                "has_ca_distance": False,
                "delta_sasa": None,
                "min_ca_distance": None,
                "nearest_ligand_residue": None,
                "distance_threshold": None,
            }
            union_map[residue_number] = entry
        else:
            if residue_name and not entry.get("residue_name"):
                entry["residue_name"] = residue_name
            if chain_id and not entry.get("chain_id"):
                entry["chain_id"] = chain_id
        return entry

    for info in interaction_rows:
        res_num = _as_int(info.get("residue_number"))
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            res_num,
            residue_name=(info.get("residue_name") or "").strip() or None,
            chain_id=(info.get("chain_id") or "").strip() or None,
        )
        entry["has_hbond"] = _as_bool(info.get("has_hbond"))
        entry["has_charge_complement"] = _as_bool(info.get("has_charge_complement"))
        entry["has_pi_contact"] = _as_bool(info.get("has_pi_contact"))
        entry["interaction_category"] = (info.get("interaction_category") or "none").strip() or "none"
        entry["hbond_count"] = _as_int(info.get("hbond_count")) or 0
        entry["charge_pair_count"] = _as_int(info.get("charge_pair_count")) or 0
        entry["pi_contact_count"] = _as_int(info.get("pi_contact_count")) or 0
        entry["is_charged"] = _as_bool(info.get("is_charged"))
        entry["charge_type"] = (info.get("charge_type") or "neutral").strip() or "neutral"

    for info in occlusion_rows:
        res_num = _as_int(info.get("residue_number"))
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            res_num,
            residue_name=(info.get("residue_name") or "").strip() or None,
            chain_id=(info.get("chain_id") or "").strip() or None,
        )
        entry["has_sasa_occlusion"] = bool(entry["has_sasa_occlusion"]) or _as_bool(info.get("is_occluded"))
        delta_sasa = _as_float(info.get("delta_sasa"))
        if delta_sasa is not None:
            entry["delta_sasa"] = delta_sasa

    for info in ca_distance_rows:
        res_num = _as_int(info.get("residue_number"))
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            res_num,
            residue_name=(info.get("residue_name") or "").strip() or None,
            chain_id=(info.get("chain_id") or "").strip() or None,
        )
        entry["has_ca_distance"] = bool(entry["has_ca_distance"]) or _as_bool(info.get("passes_filter"))
        min_ca_distance = _as_float(info.get("min_ca_distance"))
        if min_ca_distance is not None:
            entry["min_ca_distance"] = min_ca_distance
        nearest_ligand = _as_int(info.get("nearest_ligand_residue"))
        if nearest_ligand is not None:
            entry["nearest_ligand_residue"] = nearest_ligand
        distance_threshold = _as_float(info.get("distance_threshold"))
        if distance_threshold is not None:
            entry["distance_threshold"] = distance_threshold

    union_residue_info = sorted(union_map.values(), key=lambda x: int(x["residue_number"]))
    receptor_chain = None
    if union_residue_info:
        receptor_chain = union_residue_info[0].get("chain_id")

    return {
        "dataset_type": "union",
        "receptor_chain": receptor_chain,
        "residue_info": union_residue_info,
    }


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Supplementary Figure 4 from a target's per-atom CSP inputs.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--target", type=str, default="7jq8")
    parser.add_argument("--output", type=Path, default=None, help="Output path for suppl_fig_4.png")
    parser.add_argument("--title", type=str, default=None)
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    target = args.target.strip().lower()
    target_dir = args.outputs_dir / target

    csp_hn_path = target_dir / "csp_table.csv"
    csp_ca_path = target_dir / "csp_table_CA.csv"
    interaction_path = target_dir / "interaction_filter.csv"
    occlusion_path = target_dir / "occlusion_analysis.csv"
    ca_distance_path = target_dir / "ca_distance_filter.csv"

    required_files = [csp_hn_path, csp_ca_path, interaction_path, occlusion_path, ca_distance_path]
    missing = [str(path) for path in required_files if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required input files for target '{target}': {', '.join(missing)}")

    results_hn = _load_csp_results(csp_hn_path)
    results_ca = _load_csp_results(csp_ca_path)
    binding_results = _build_union_binding_results(interaction_path, occlusion_path, ca_distance_path)

    output_path = args.output or (args.figures_dir / "suppl_fig_4.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    title = args.title if args.title is not None else f"{target} per-atom perturbations"

    plot_per_atom_classification_panels(
        results_hn=results_hn,
        results_ca=results_ca,
        binding_results=binding_results,
        out_png=str(output_path),
        title=title,
        pdb_id=target,
        structure_pdb_path=str(target_dir / f"{target}_csp.pdb"),
        receptor_chain=str(binding_results.get("receptor_chain") or ""),
    )
    print(f"Supplementary Figure 4 written to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
