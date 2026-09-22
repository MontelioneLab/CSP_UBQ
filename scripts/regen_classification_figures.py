#!/usr/bin/env python3
"""Rebuild per-target CSP classification PNGs and recompose case-study figures.

Reads existing csp_table / master_alignment CSVs; does not recompute CSPs or
rerun the pipeline. Reuses cached PyMOL case-study panels (no interactive view).

Usage:
    python scripts/regen_classification_figures.py
    python scripts/regen_classification_figures.py --outputs outputs --ids 2FIN_6809
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg", force=True)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.case_study import generate_case_study_figure  # noqa: E402
from scripts.case_study_2 import generate_case_study_2_figure  # noqa: E402
from scripts.case_study_3 import generate_case_study_3_figure  # noqa: E402
from scripts.case_study_z import generate_case_study_z_figure  # noqa: E402
from scripts.csp import CSPResult  # noqa: E402
from scripts.plot_csp_z_vs_ca_distance import plot_per_target  # noqa: E402
from scripts.regen_hsqc_case_study_offsets import (  # noqa: E402
    _b,
    _f,
    load_binding_results,
    load_metadata_index,
)
from scripts.visualize import AA_THREE_TO_ONE, plot_csp_classification_bars, plot_per_atom_classification_panels  # noqa: E402

ONE_TO_THREE = {one: three for three, one in AA_THREE_TO_ONE.items()}


def _first_float(row: Dict[str, str], *keys: str) -> Optional[float]:
    for key in keys:
        val = _f(row.get(key))
        if val is not None:
            return val
    return None


def _first_bool(row: Dict[str, str], *keys: str) -> Optional[bool]:
    for key in keys:
        val = _b(row.get(key))
        if val is not None:
            return val
    return None


def load_csp_table(path: Path) -> List[CSPResult]:
    """Load HN, CA, or HA/CA CSP tables into CSPResult rows for classification plots."""
    results: List[CSPResult] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            apo_resi = _f(row.get("apo_resi"))
            holo_resi = _f(row.get("holo_resi")) or _f(row.get("residue_number"))
            if holo_resi is None:
                continue
            if apo_resi is None:
                apo_resi = holo_resi
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
                    csp_A=_first_float(row, "csp_A", "csp_CA", "csp_HA_CA"),
                    significant=_first_bool(
                        row, "significant", "csp_CA_significant", "csp_HA_CA_significant"
                    ),
                    significant_1sd=_first_bool(
                        row,
                        "significant_1sd",
                        "csp_CA_significant_1sd",
                        "csp_HA_CA_significant_1sd",
                    ),
                    significant_2sd=_first_bool(
                        row,
                        "significant_2sd",
                        "csp_CA_significant_2sd",
                        "csp_HA_CA_significant_2sd",
                    ),
                    H_holo_original=_f(row.get("H_holo_original")) or _f(row.get("H_holo")),
                    N_holo_original=_f(row.get("N_holo_original")) or _f(row.get("N_holo")),
                    H_offset=_f(row.get("H_offset")),
                    N_offset=_f(row.get("N_offset")),
                    z_score=_first_float(row, "csp_z", "csp_CA_z", "csp_HA_CA_z"),
                    CA_apo=_f(row.get("CA_apo")),
                    CA_holo=_f(row.get("CA_holo")),
                    CA_holo_original=_f(row.get("CA_holo_original")) or _f(row.get("CA_holo")),
                    CA_offset=_f(row.get("CA_offset")),
                    dCA=_f(row.get("dCA")),
                    HA_apo=_f(row.get("HA_apo")),
                    HA_holo=_f(row.get("HA_holo")),
                    HA_holo_original=_f(row.get("HA_holo_original")) or _f(row.get("HA_holo")),
                    HA_offset=_f(row.get("HA_offset")),
                    dHA=_f(row.get("dHA")),
                )
            )
    return results


def load_plot_binding_results(target_dir: Path) -> Optional[dict]:
    """Binding lookup from master_alignment, with 3-letter residue names for alignment."""
    binding = load_binding_results(target_dir)
    if not binding:
        return None
    for info in binding.get("residue_info") or []:
        name = (info.get("residue_name") or "").strip()
        if len(name) == 1:
            info["residue_name"] = ONE_TO_THREE.get(name.upper(), name)
    return binding


def ids_from_csp_table(target_dir: Path) -> Dict[str, str]:
    table = target_dir / "csp_table.csv"
    if not table.is_file():
        return {}
    with table.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), None)
    if not row:
        return {}
    return {
        "apo_bmrb": (row.get("apo_bmrb") or "").strip(),
        "holo_bmrb": (row.get("holo_bmrb") or "").strip(),
        "holo_pdb": (row.get("holo_pdb") or "").strip(),
    }


def _glob_one(target_dir: Path, pattern: str) -> bool:
    return any(target_dir.glob(pattern))


def process_target(
    target_dir: Path,
    meta_index: Dict[Tuple[str, str], Dict[str, str]],
) -> List[str]:
    """Rebuild classification plots and recompose existing case-study figures.

    Returns a list of short action labels for logging.
    """
    csp_path = target_dir / "csp_table.csv"
    if not csp_path.is_file():
        raise RuntimeError("missing csp_table.csv")

    results_hn = load_csp_table(csp_path)
    if not results_hn:
        raise RuntimeError("no CSP rows")

    binding = load_plot_binding_results(target_dir)
    if binding is None:
        raise RuntimeError("missing master_alignment.csv binding rows")

    parts = target_dir.name.split("_", 1)
    dirname_pdb = parts[0]
    dirname_apo = parts[1] if len(parts) > 1 else ""
    ids = ids_from_csp_table(target_dir)
    holo_pdb = ids.get("holo_pdb") or dirname_pdb
    apo_bmrb = ids.get("apo_bmrb") or dirname_apo
    meta = meta_index.get((holo_pdb.upper(), apo_bmrb), {})
    holo_bmrb = (ids.get("holo_bmrb") or meta.get("holo_bmrb") or "").strip()
    apo_pdb = (meta.get("apo_pdb") or "").strip() or None

    done: List[str] = []

    plot_csp_classification_bars(
        results_hn,
        binding,
        str(target_dir / "csp_classification_bars_original.png"),
        title=f"{holo_pdb} CSP Classification (original)",
        significance_field="significant",
    )
    done.append("hn_bars")

    ca_table = target_dir / "csp_table_CA.csv"
    if ca_table.is_file():
        results_ca = load_csp_table(ca_table)
        if results_ca:
            plot_csp_classification_bars(
                results_ca,
                binding,
                str(target_dir / "csp_classification_bars_original_CA.png"),
                title=f"{holo_pdb} CSP Classification (CA, original)",
                significance_field="significant",
            )
            done.append("ca_bars")
    else:
        results_ca = None

    ha_ca_table = target_dir / "csp_table_HA_CA.csv"
    if ha_ca_table.is_file():
        results_ha_ca = load_csp_table(ha_ca_table)
        if results_ha_ca:
            plot_csp_classification_bars(
                results_ha_ca,
                binding,
                str(target_dir / "csp_classification_bars_original_HA_CA.png"),
                title=f"{holo_pdb} CSP Classification (HA/CA, original)",
                significance_field="significant",
            )
            done.append("ha_ca_bars")
    else:
        results_ha_ca = None

    per_atom_path = target_dir / "per_atom_classification_panels.png"
    if per_atom_path.is_file():
        structure_pdb = target_dir / f"{holo_pdb}_csp.pdb"
        if not structure_pdb.is_file():
            structure_pdb = target_dir / f"{holo_pdb.lower()}_csp.pdb"
        plot_per_atom_classification_panels(
            results_hn=results_hn,
            results_ca=results_ca,
            results_ha_ca=results_ha_ca,
            binding_results=binding,
            out_png=str(per_atom_path),
            title=f"{holo_pdb} per-atom perturbations",
            pdb_id=holo_pdb,
            structure_pdb_path=str(structure_pdb) if structure_pdb.is_file() else None,
        )
        done.append("per_atom")

    if (target_dir / "master_alignment.csv").is_file():
        scatter = plot_per_target(str(target_dir))
        if scatter:
            done.append("scatter")

    cs_kwargs: Dict[str, Any] = dict(
        target_dir=str(target_dir),
        pdb_id=holo_pdb,
        apo_bmrb=apo_bmrb or None,
        holo_bmrb=holo_bmrb or None,
        apo_pdb=apo_pdb,
        view_key=target_dir.name,
        reuse_existing_panels=True,
        allow_interactive_view=False,
    )

    if _glob_one(target_dir, "*_case_study.png"):
        generate_case_study_figure(**cs_kwargs)
        done.append("case_study")
    if _glob_one(target_dir, "*_case_study_2.png"):
        generate_case_study_2_figure(**cs_kwargs)
        done.append("case_study_2")
    if _glob_one(target_dir, "*_case_study_3.png"):
        if not apo_bmrb or not holo_bmrb:
            raise RuntimeError("case_study_3 needs apo_bmrb and holo_bmrb")
        generate_case_study_3_figure(
            **cs_kwargs,
            do_fetch=False,
        )
        done.append("case_study_3")
    if _glob_one(target_dir, "*_case_study_z.png"):
        generate_case_study_z_figure(**cs_kwargs)
        done.append("case_study_z")

    return done


def discover_targets(outputs_dir: Path, ids: Sequence[str]) -> List[Path]:
    wanted = {x.strip() for x in ids if x.strip()}
    targets: List[Path] = []
    for path in sorted(outputs_dir.iterdir(), key=lambda p: p.name.lower()):
        if not path.is_dir() or path.name.startswith("."):
            continue
        if wanted and path.name not in wanted:
            continue
        if (path / "csp_table.csv").is_file():
            targets.append(path)
    return targets


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=PROJECT_ROOT / "outputs",
        help="Outputs root containing per-target folders (default: outputs)",
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
        help="Optional comma-separated target folder names (e.g. 2FIN_6809)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    outputs_dir = args.outputs if args.outputs.is_absolute() else PROJECT_ROOT / args.outputs
    input_csv = args.input if args.input.is_absolute() else PROJECT_ROOT / args.input
    ids = [x.strip() for x in args.ids.split(",") if x.strip()]
    targets = discover_targets(outputs_dir, ids)
    if not targets:
        print(f"No target directories with csp_table.csv under {outputs_dir}")
        return 1

    meta_index = load_metadata_index(input_csv)
    ok = 0
    fail = 0
    failures: List[str] = []
    print(f"Regenerating classification figures for {len(targets)} targets under {outputs_dir}")
    for i, target_dir in enumerate(targets, start=1):
        print(f"[{i}/{len(targets)}] {target_dir.name}", flush=True)
        try:
            done = process_target(target_dir, meta_index)
            print(f"  OK {target_dir.name}: {', '.join(done)}", flush=True)
            ok += 1
        except Exception as exc:
            fail += 1
            failures.append(f"{target_dir.name}: {exc}")
            print(f"  FAIL {target_dir.name}: {exc}", flush=True)

    print("=" * 60)
    print(f"DONE ok={ok} fail={fail} targets={len(targets)}")
    if failures:
        print("Failures:")
        for line in failures:
            print(f"  {line}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
