#!/usr/bin/env python3
"""
Write Supplementary_Data_1.csv–Supplementary_Data_6.csv from the pipeline
outputs already used by Figures 1–3.

- Supplementary_Data_1: tidy ¹H/¹⁵N offset-grid table for Figure 1a (PDB 7JQ8)
- Supplementary_Data_2–5: per-residue tables for Figure 2 rows A–D
- Supplementary_Data_6: per-system TP/FP/TN/FN counts for the Figure 3 cohort

Example::

    python scripts/create_supplementary_data.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .config import Referencing
    from .create_fig_1 import resolve_fig1_panel_a_csv_from_pipeline
    from .create_fig_2 import FIG2_TARGETS, Fig2Target
    from .create_fig_3 import (
        CA_DISTANCE_COLUMN,
        CSP_COLUMN,
        PREDICTOR_COLUMNS,
        SIGNIFICANT_COLUMN,
        _as_bool,
    )
    from .csp import load_grid_csv
    from .merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from .target_resolution import (
        build_resolution_caches,
        load_target_rows,
        resolve_row,
        resolve_target_rows,
    )
except Exception:
    _sys_path = str(_REPO_ROOT)
    if _sys_path not in sys.path:
        sys.path.insert(0, _sys_path)
    from scripts.config import Referencing  # type: ignore
    from scripts.create_fig_1 import resolve_fig1_panel_a_csv_from_pipeline  # type: ignore
    from scripts.create_fig_2 import FIG2_TARGETS, Fig2Target  # type: ignore
    from scripts.create_fig_3 import (  # type: ignore
        CA_DISTANCE_COLUMN,
        CSP_COLUMN,
        PREDICTOR_COLUMNS,
        SIGNIFICANT_COLUMN,
        _as_bool,
    )
    from scripts.csp import load_grid_csv  # type: ignore
    from scripts.merge_csv import filter_recorded_csp_dataframe, parse_optional_bool  # type: ignore
    from scripts.target_resolution import (  # type: ignore
        build_resolution_caches,
        load_target_rows,
        resolve_row,
        resolve_target_rows,
    )


FIG2_SUPPLEMENTARY_FILES: Sequence[Tuple[str, Fig2Target]] = (
    ("Supplementary_Data_2.csv", FIG2_TARGETS[0]),
    ("Supplementary_Data_3.csv", FIG2_TARGETS[1]),
    ("Supplementary_Data_4.csv", FIG2_TARGETS[2]),
    ("Supplementary_Data_5.csv", FIG2_TARGETS[3]),
)

CASE_STUDY_COLUMNS: Sequence[str] = (
    "Residue",
    "CSP magnitude",
    "CSP Mask",
    "Binding Site Mask",
    "CSP Classification",
)

FIG3_COUNT_COLUMNS: Sequence[str] = (
    "apo_pdb",
    "apo_bmrb",
    "holo_pdb",
    "holo_bmrb",
    "TP",
    "FP",
    "TN",
    "FN",
)

VALID_CLASSIFICATIONS = frozenset({"TP", "FP", "TN", "FN"})


def _as_bool_mask(value: object) -> bool:
    parsed = parse_optional_bool(value)
    return bool(parsed)


def _parse_grid_comment_metadata(grid_csv: Path) -> Dict[str, str]:
    """Read ``# key,value`` comment lines from a pipeline offset-grid CSV."""
    meta: Dict[str, str] = {}
    with grid_csv.open(newline="", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped.startswith("#"):
                continue
            body = stripped[1:].strip()
            if "," not in body:
                continue
            key, value = body.split(",", 1)
            key = key.strip()
            value = value.strip()
            if key:
                meta[key] = value
    return meta


def write_supplementary_data_1(
    grid_csv: Path,
    output_path: Path,
    *,
    holo_pdb: str,
    apo_bmrb: str,
    holo_bmrb: str,
) -> int:
    """Write a tidy H/N offset-grid table from the Figure 1a heatmap CSV."""
    grid_result = load_grid_csv(str(grid_csv))
    if grid_result is None:
        raise ValueError(f"Failed to parse offset-grid CSV: {grid_csv}")

    h_values = list(grid_result["h_values"])
    n_values = list(grid_result["n_values"])
    count_matrix = list(grid_result["count_matrix"])
    comment_meta = _parse_grid_comment_metadata(grid_csv)
    cutoff = comment_meta.get("cutoff", str(Referencing().grid_cutoff))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        try:
            source = grid_csv.resolve().relative_to(_REPO_ROOT.resolve()).as_posix()
        except ValueError:
            source = grid_csv.as_posix()
        handle.write(f"# source,{source}\n")
        handle.write(f"# holo_pdb,{holo_pdb}\n")
        handle.write(f"# apo_bmrb,{apo_bmrb}\n")
        handle.write(f"# holo_bmrb,{holo_bmrb}\n")
        handle.write(f"# cutoff,{cutoff}\n")
        handle.write(f"# best_h_offset,{grid_result['best_h_offset']}\n")
        handle.write(f"# best_n_offset,{grid_result['best_n_offset']}\n")
        handle.write(f"# best_count,{grid_result['best_count']}\n")
        writer = csv.writer(handle)
        writer.writerow(["h_offset_ppm", "n_offset_ppm", "n_aligned_peaks"])
        for i, h_off in enumerate(h_values):
            if i >= len(count_matrix):
                break
            row_counts = count_matrix[i]
            for j, n_off in enumerate(n_values):
                if j >= len(row_counts):
                    break
                writer.writerow([h_off, n_off, int(row_counts[j])])
                n_rows += 1
    return n_rows


def _resolve_fig2_target_dir(
    target: Fig2Target,
    outputs_index: Dict[str, Path],
    bmrb_cache,
) -> Path:
    path = resolve_row(target.row, outputs_index, bmrb_cache)
    if path is None:
        raise FileNotFoundError(
            "No outputs subdirectory for "
            f"holo_pdb={target.holo_pdb} apo_bmrb={target.apo_bmrb} "
            f"holo_bmrb={target.holo_bmrb}."
        )
    alignment = path / "master_alignment.csv"
    if not alignment.is_file():
        raise FileNotFoundError(f"Missing master_alignment.csv in {path}")
    return path


def collect_figure_2_rows(alignment_path: Path) -> List[Dict[str, object]]:
    """Residue rows shown in the Figure 2 CSP classification bar plot."""
    records: List[Dict[str, object]] = []
    with alignment_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cls = (row.get("classification") or "").strip()
            if cls not in VALID_CLASSIFICATIONS:
                continue
            csp_str = (row.get("csp_A") or "").strip()
            if not csp_str:
                continue
            try:
                float(csp_str)
            except ValueError:
                continue
            try:
                residue = int(float(row.get("holo_resi", "")))
            except (TypeError, ValueError):
                continue
            records.append(
                {
                    "Residue": residue,
                    "CSP magnitude": csp_str,
                    "CSP Mask": _as_bool_mask(row.get("significant", "")),
                    "Binding Site Mask": _as_bool_mask(row.get("is_binding_site", "")),
                    "CSP Classification": cls,
                }
            )
    return records


def write_figure_2_table(alignment_path: Path, output_path: Path) -> int:
    records = collect_figure_2_rows(alignment_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CASE_STUDY_COLUMNS))
        writer.writeheader()
        writer.writerows(records)
    return len(records)


def _load_fig3_targets_meta(targets_csv: Path) -> Dict[str, Dict[str, str]]:
    """Map canonical ``HOLO_apoBmrb`` system id to apo/holo identifiers."""
    meta: Dict[str, Dict[str, str]] = {}
    with targets_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = {fn.lower(): fn for fn in (reader.fieldnames or [])}
        for raw in reader:
            holo = (raw.get(fields.get("holo_pdb", "holo_pdb"), "") or "").strip()
            apo_bmrb = (raw.get(fields.get("apo_bmrb", "apo_bmrb"), "") or "").strip()
            if not holo or not apo_bmrb:
                continue
            apo_pdb = ""
            if "apo_pdb" in fields:
                apo_pdb = (raw.get(fields["apo_pdb"], "") or "").strip()
            holo_bmrb = ""
            if "holo_bmrb" in fields:
                holo_bmrb = (raw.get(fields["holo_bmrb"], "") or "").strip()
            meta[f"{holo.upper()}_{apo_bmrb}"] = {
                "apo_pdb": apo_pdb,
                "apo_bmrb": apo_bmrb,
                "holo_pdb": holo.upper(),
                "holo_bmrb": holo_bmrb,
            }
    return meta


def count_figure_3_confusion(alignment_path: Path) -> Optional[Dict[str, int]]:
    """TP/FP/TN/FN counts using the same residue filters as Figure 3 histograms."""
    df = pd.read_csv(alignment_path)
    required_cols = {SIGNIFICANT_COLUMN, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS}
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"[WARN] Skipping {alignment_path}: missing columns {missing}")
        return None

    df = filter_recorded_csp_dataframe(
        df, csp_column=CSP_COLUMN, significant_column=SIGNIFICANT_COLUMN
    )
    if df.empty:
        return None

    is_significant = df[SIGNIFICANT_COLUMN].map(_as_bool)
    is_binding = pd.DataFrame({c: df[c].map(_as_bool) for c in PREDICTOR_COLUMNS}).any(
        axis=1
    )
    distances = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")
    valid = distances.notna()
    if not valid.any():
        return None

    sig = is_significant[valid]
    bind = is_binding[valid]
    return {
        "TP": int((sig & bind).sum()),
        "FP": int((sig & ~bind).sum()),
        "FN": int((~sig & bind).sum()),
        "TN": int((~sig & ~bind).sum()),
    }


def _ids_from_alignment(alignment_path: Path) -> Dict[str, str]:
    with alignment_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            return {
                "apo_bmrb": (row.get("apo_bmrb") or "").strip(),
                "holo_bmrb": (row.get("holo_bmrb") or "").strip(),
                "holo_pdb": (row.get("holo_pdb") or "").strip().upper(),
            }
    return {}


def collect_figure_3_rows(
    outputs_dir: Path,
    targets_csv: Path,
) -> List[Dict[str, object]]:
    rows = load_target_rows(targets_csv)
    target_paths = resolve_target_rows(rows, outputs_dir)
    meta = _load_fig3_targets_meta(targets_csv)
    records: List[Dict[str, object]] = []
    for target_dir in target_paths:
        alignment_path = target_dir / "master_alignment.csv"
        if not alignment_path.is_file():
            continue
        counts = count_figure_3_confusion(alignment_path)
        if counts is None:
            continue
        info = dict(meta.get(target_dir.name, {}))
        if not info:
            info.update(_ids_from_alignment(alignment_path))
            info.setdefault("apo_pdb", "")
            info.setdefault("holo_pdb", target_dir.name.split("_")[0].upper())
        records.append(
            {
                "apo_pdb": info.get("apo_pdb", ""),
                "apo_bmrb": info.get("apo_bmrb", ""),
                "holo_pdb": info.get("holo_pdb", ""),
                "holo_bmrb": info.get("holo_bmrb", ""),
                "TP": counts["TP"],
                "FP": counts["FP"],
                "TN": counts["TN"],
                "FN": counts["FN"],
            }
        )
    records.sort(
        key=lambda rec: (
            str(rec.get("holo_pdb", "")).upper(),
            str(rec.get("apo_bmrb", "")),
        )
    )
    return records


def write_figure_3_table(records: Sequence[Dict[str, object]], output_path: Path) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FIG3_COUNT_COLUMNS))
        writer.writeheader()
        writer.writerows(records)
    return len(records)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write Supplementary_Data_1.csv–Supplementary_Data_6.csv for Figures 1–3."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=_REPO_ROOT / "outputs",
        help="Pipeline outputs root (default: outputs).",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=_REPO_ROOT / "figures",
        help="Directory for Supplementary_Data_*.csv (default: figures).",
    )
    parser.add_argument(
        "--fig1-targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ.csv",
        help="Targets CSV used to resolve Figure 1a (default: data/CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--fig3-targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ_ph0.5_temp5C.csv",
        help="Figure 3 cohort CSV (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    figures_dir = args.figures_dir if args.figures_dir.is_absolute() else _REPO_ROOT / args.figures_dir
    fig1_targets = (
        args.fig1_targets_csv
        if args.fig1_targets_csv.is_absolute()
        else _REPO_ROOT / args.fig1_targets_csv
    )
    fig3_targets = (
        args.fig3_targets_csv
        if args.fig3_targets_csv.is_absolute()
        else _REPO_ROOT / args.fig3_targets_csv
    )

    fig1_target = FIG2_TARGETS[0]
    grid_csv = resolve_fig1_panel_a_csv_from_pipeline(outputs_dir, fig1_targets)
    if grid_csv is None or not grid_csv.is_file():
        print("Error: Figure 1a offset-grid CSV not found.", file=sys.stderr)
        return 1

    sd1_path = figures_dir / "Supplementary_Data_1.csv"
    n_grid = write_supplementary_data_1(
        grid_csv,
        sd1_path,
        holo_pdb=fig1_target.holo_pdb,
        apo_bmrb=fig1_target.apo_bmrb,
        holo_bmrb=fig1_target.holo_bmrb,
    )
    print(f"Wrote {sd1_path} ({n_grid} grid points)")

    outputs_index, bmrb_cache = build_resolution_caches(outputs_dir)
    for filename, target in FIG2_SUPPLEMENTARY_FILES:
        target_dir = _resolve_fig2_target_dir(target, outputs_index, bmrb_cache)
        out_path = figures_dir / filename
        n_rows = write_figure_2_table(target_dir / "master_alignment.csv", out_path)
        print(
            f"Wrote {out_path} ({n_rows} residues; "
            f"{target.panel} {target.holo_pdb})"
        )

    sd6_path = figures_dir / "Supplementary_Data_6.csv"
    fig3_records = collect_figure_3_rows(outputs_dir, fig3_targets)
    if not fig3_records:
        print("Error: no Figure 3 systems produced confusion counts.", file=sys.stderr)
        return 1
    n_systems = write_figure_3_table(fig3_records, sd6_path)
    print(f"Wrote {sd6_path} ({n_systems} systems)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
