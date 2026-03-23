#!/usr/bin/env python3
"""
Confusion-matrix stacked histograms split by apo/holo buffer similarity (pH, temperature).

Uses the same definitions as analyze_targets(_ca): TP/FP/FN/TN vs minimum CA distance.

Rows are classified using apo_holo_exp_conditions.csv (or recomputed from the same rules):
  - "similar": both apo and holo report pH and temperature, |ΔpH| ≤ --ph-max-diff (default 0.5),
    |ΔT| ≤ --temp-max-diff-c (default 5 °C)
  - "complementary": all other rows (missing metrics or outside tolerances)

Stricter pH match (e.g. |ΔpH| ≤ 0.1):  --ph-max-diff 0.1

Each CSP_UBQ row is mapped to its pipeline output folder (including holo_pdb_1, _2, … when
duplicates exist). Folder names are matched case-insensitively against outputs/.

Default mode is CA-inclusive CSPs (csp_table_CA.csv). Use --mode nh for N/H (master_alignment.csv).

Examples:
  python3 scripts/plot_confusion_histograms_by_exp_conditions.py \\
    --output-dir outputs/summary_statistics
  python3 scripts/plot_confusion_histograms_by_exp_conditions.py \\
    --ph-max-diff 0.1 --output-dir outputs/summary_statistics
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.config import Paths  # noqa: E402

DEFAULT_PH_MAX_DIFF = 0.5
DEFAULT_TEMP_MAX_DIFF_C = 5.0


def load_csp_ubq_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def logical_output_dirname_per_row(rows: Sequence[Dict[str, str]]) -> List[Optional[str]]:
    """Same holo_pdb + duplicate-index scheme as scripts/pipeline.py."""
    holo_pdb_to_rows: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        h = (row.get("holo_pdb") or "").strip()
        if h:
            holo_pdb_to_rows.setdefault(h, []).append(row)

    out: List[Optional[str]] = []
    for row in rows:
        h = (row.get("holo_pdb") or "").strip()
        apo = (row.get("apo_bmrb") or "").strip()
        holo_b = (row.get("holo_bmrb") or "").strip()
        if not h:
            out.append(None)
            continue
        matches = holo_pdb_to_rows[h]
        if len(matches) <= 1:
            out.append(h)
            continue
        duplicate_index: Optional[int] = None
        for idx, m in enumerate(matches, start=1):
            if (m.get("apo_bmrb") or "").strip() == apo and (m.get("holo_bmrb") or "").strip() == holo_b:
                duplicate_index = idx
                break
        if duplicate_index is not None:
            out.append(f"{h}_{duplicate_index}")
        else:
            out.append(None)
    return out


def build_outputs_index(outputs_dir: Path) -> Dict[str, Path]:
    """Lowercase directory name -> actual path (first wins)."""
    idx: Dict[str, Path] = {}
    if not outputs_dir.is_dir():
        return idx
    for p in sorted(outputs_dir.iterdir()):
        if p.is_dir() and not p.name.startswith("."):
            k = p.name.lower()
            if k not in idx:
                idx[k] = p
    return idx


def resolve_dir(index: Dict[str, Path], logical: Optional[str]) -> Optional[Path]:
    if not logical:
        return None
    return index.get(logical.lower())


def row_meets_ph_temp_criteria(
    exp_row: Dict[str, str],
    *,
    ph_max_diff: float = DEFAULT_PH_MAX_DIFF,
    temp_max_diff_c: float = DEFAULT_TEMP_MAX_DIFF_C,
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
    return abs(aph - hph) <= ph_max_diff and abs(at - ht) <= temp_max_diff_c


def collect_nh_confusion(target_dirs: Sequence[Path]) -> List[Any]:
    """Lazy-imports analyze_targets (requires matplotlib)."""
    from scripts.analyze_targets import (  # noqa: WPS433
        CA_DISTANCE_COLUMN,
        PREDICTOR_COLUMNS,
        SIGNIFICANT_COLUMN,
        AlignmentParsingError,
        ConfusionRecord,
        load_alignment,
    )

    records: List[Any] = []
    for target_dir in target_dirs:
        alignment_path = target_dir / "master_alignment.csv"
        if not alignment_path.exists():
            continue
        try:
            df = load_alignment(alignment_path)
        except AlignmentParsingError:
            continue
        predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)
        if SIGNIFICANT_COLUMN not in df.columns or CA_DISTANCE_COLUMN not in df.columns:
            continue
        mask = df[CA_DISTANCE_COLUMN].notna()
        if not mask.any():
            continue
        subset = df.loc[mask]
        for idx, row in subset.iterrows():
            records.append(
                ConfusionRecord(
                    distance=float(row[CA_DISTANCE_COLUMN]),
                    is_binding=bool(predicted.loc[idx]),
                    is_significant=bool(row[SIGNIFICANT_COLUMN]),
                )
            )
    return records


def collect_ca_confusion(target_dirs: Sequence[Path]) -> List[Any]:
    """Lazy-imports analyze_targets_ca (requires matplotlib)."""
    from scripts.analyze_targets_ca import (  # noqa: WPS433
        CA_DISTANCE_COLUMN,
        PREDICTOR_COLUMNS,
        SIGNIFICANT_COLUMN,
        AlignmentParsingError,
        ConfusionRecord,
        load_ca_alignment,
    )

    records: List[Any] = []
    for target_dir in target_dirs:
        if not (target_dir / "csp_table_CA.csv").exists():
            continue
        try:
            df = load_ca_alignment(target_dir)
        except AlignmentParsingError:
            continue
        available = [c for c in PREDICTOR_COLUMNS if c in df.columns]
        if not available:
            continue
        if SIGNIFICANT_COLUMN not in df.columns or CA_DISTANCE_COLUMN not in df.columns:
            continue
        predicted = df[available].any(axis=1)
        mask = df[CA_DISTANCE_COLUMN].notna()
        if not mask.any():
            continue
        subset = df.loc[mask]
        for idx, row in subset.iterrows():
            records.append(
                ConfusionRecord(
                    distance=float(row[CA_DISTANCE_COLUMN]),
                    is_binding=bool(predicted.loc[idx]),
                    is_significant=bool(row[SIGNIFICANT_COLUMN]),
                )
            )
    return records


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Confusion-matrix stacked histograms for similar vs complementary pH/temperature subsets."
    )
    p.add_argument(
        "--csp-csv",
        type=Path,
        default=None,
        help="CSP_UBQ-style CSV with apo_bmrb, holo_bmrb, holo_pdb (default: config Paths.input_csv).",
    )
    p.add_argument(
        "--exp-conditions-csv",
        type=Path,
        default=None,
        help="apo_holo_exp_conditions.csv (default: repo root apo_holo_exp_conditions.csv).",
    )
    p.add_argument(
        "--outputs-dir",
        type=Path,
        default=None,
        help="Pipeline outputs root (default: config Paths.outputs_dir).",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for PNG outputs (default: outputs/summary_statistics).",
    )
    p.add_argument(
        "--mode",
        choices=("ca", "nh"),
        default="ca",
        help="CA-inclusive (csp_table_CA) or N/H (master_alignment) confusion data.",
    )
    p.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional CSV listing each CSP row, dirname, subset, and whether output dir resolved.",
    )
    p.add_argument(
        "--ph-max-diff",
        type=float,
        default=DEFAULT_PH_MAX_DIFF,
        metavar="DELTA",
        help=f"Max |apo pH − holo pH| for 'similar' subset (default: {DEFAULT_PH_MAX_DIFF}).",
    )
    p.add_argument(
        "--temp-max-diff-c",
        type=float,
        default=DEFAULT_TEMP_MAX_DIFF_C,
        metavar="DELTA_C",
        help=f"Max |Δtemperature| in °C for 'similar' subset (default: {DEFAULT_TEMP_MAX_DIFF_C}).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    cfg = Paths()
    csp_path = args.csp_csv or Path(cfg.input_csv)
    exp_path = args.exp_conditions_csv or (_REPO_ROOT / "apo_holo_exp_conditions.csv")
    outputs_root = args.outputs_dir or Path(cfg.outputs_dir)
    out_dir = args.output_dir or (Path(cfg.outputs_dir) / "summary_statistics")

    if not csp_path.is_file():
        print(f"Missing CSP CSV: {csp_path}", file=sys.stderr)
        return 1
    if not exp_path.is_file():
        print(f"Missing experimental conditions CSV: {exp_path}", file=sys.stderr)
        return 1
    if not outputs_root.is_dir():
        print(f"Missing outputs directory: {outputs_root}", file=sys.stderr)
        return 1

    csp_rows = load_csp_ubq_rows(csp_path)
    dirnames = logical_output_dirname_per_row(csp_rows)

    with exp_path.open(newline="", encoding="utf-8") as f:
        exp_rows = list(csv.DictReader(f))

    if len(exp_rows) != len(csp_rows):
        print(
            f"[WARN] Row count mismatch: {csp_path.name} ({len(csp_rows)}) vs "
            f"{exp_path.name} ({len(exp_rows)}); pairing by row index anyway.",
            file=sys.stderr,
        )
    n = min(len(csp_rows), len(exp_rows))

    index = build_outputs_index(outputs_root)
    similar_dirs: List[Path] = []
    complement_dirs: List[Path] = []
    manifest_rows: List[Dict[str, str]] = []

    seen_similar: Set[str] = set()
    seen_comp: Set[str] = set()

    ph_tol = float(args.ph_max_diff)
    temp_tol = float(args.temp_max_diff_c)

    for i in range(n):
        logical = dirnames[i] if i < len(dirnames) else None
        resolved = resolve_dir(index, logical)
        similar = row_meets_ph_temp_criteria(
            exp_rows[i], ph_max_diff=ph_tol, temp_max_diff_c=temp_tol
        )
        subset = "similar" if similar else "complementary"
        manifest_rows.append(
            {
                "row_index": str(i),
                "logical_output_dir": logical or "",
                "resolved_path": str(resolved) if resolved else "",
                "subset": subset,
                "ph_max_diff": str(ph_tol),
                "temp_max_diff_c": str(temp_tol),
            }
        )
        if resolved is None:
            continue
        key = str(resolved.resolve())
        if similar:
            if key not in seen_similar:
                seen_similar.add(key)
                similar_dirs.append(resolved)
        else:
            if key not in seen_comp:
                seen_comp.add(key)
                complement_dirs.append(resolved)

    suffix = "ca" if args.mode == "ca" else "nh"
    out_dir = out_dir if out_dir.is_absolute() else _REPO_ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    sim_png = out_dir / (
        f"confusion_matrix_stacked_histogram_similar_ph{ph_tol:g}_temp{temp_tol:g}_{suffix}.png"
    )
    comp_png = out_dir / (
        f"confusion_matrix_stacked_histogram_complement_ph{ph_tol:g}_temp{temp_tol:g}_{suffix}.png"
    )

    if args.mode == "ca":
        from scripts.analyze_targets_ca import (  # noqa: WPS433
            render_confusion_matrix_stacked_histogram as render_confusion_ca,
        )

        sim_recs = collect_ca_confusion(similar_dirs)
        comp_recs = collect_ca_confusion(complement_dirs)
        render_confusion_ca(sim_recs, sim_png, show_title=True)
        render_confusion_ca(comp_recs, comp_png, show_title=True)
    else:
        from scripts.analyze_targets import (  # noqa: WPS433
            render_confusion_matrix_stacked_histogram as render_confusion_nh,
        )

        sim_recs = collect_nh_confusion(similar_dirs)
        comp_recs = collect_nh_confusion(complement_dirs)
        render_confusion_nh(sim_recs, sim_png, show_title=True)
        render_confusion_nh(comp_recs, comp_png, show_title=True)

    if args.manifest and manifest_rows:
        man_path = args.manifest if args.manifest.is_absolute() else _REPO_ROOT / args.manifest
        man_path.parent.mkdir(parents=True, exist_ok=True)
        with man_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
            w.writeheader()
            w.writerows(manifest_rows)

    print(f"Similar subset: {len(similar_dirs)} targets with resolved dirs, {len(sim_recs)} confusion records -> {sim_png}")
    print(f"Complementary subset: {len(complement_dirs)} targets, {len(comp_recs)} records -> {comp_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
