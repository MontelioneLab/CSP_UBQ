#!/usr/bin/env python3
"""Count false-positive CSPs by DSSP secondary structure (helix / sheet / loop).

Defaults match the Figure 3 / paper cohort:
  - targets from ``data/CSP_UBQ_ph0.5_temp5C.csv`` (140 systems)
  - significance: ``significant_max_05_cleaned_mean`` (= max(cleaned mean, 0.05))
  - only residues with a recorded CSP and finite ``min_ca_distance_distance``
  - binding-site union = occlusion ∪ CA-distance filter ∪ interactions
    (same predictors as ``create_fig_3_thresholds.collect_distance_categories``)

Requires the ``ss`` column written by ``backfill_dssp_ss.py``.

Writes:
  - <outputs>/fp_by_ss_summary.csv          aggregate counts
  - <outputs>/fp_by_ss_by_target.csv        per-target counts
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.merge_csv import has_recorded_csp, parse_optional_bool
from scripts.secondary_structure import SS_NAME
from scripts.target_resolution import load_target_rows, resolve_target_rows

SIGNIFICANT_COLUMN = "significant_max_05_cleaned_mean"
SS_COLUMN = "ss"
CSP_COLUMN = "csp_A"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
DEFAULT_TARGETS_CSV = _ROOT / "data" / "CSP_UBQ_ph0.5_temp5C.csv"
SS_ORDER = ("H", "E", "C", "")  # unknown last

# Same binding-site predictors as create_fig_3_thresholds.py
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


def find_master_alignment_files(outputs_dir: Path) -> List[Path]:
    if not outputs_dir.is_dir():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    return sorted(outputs_dir.glob("*/master_alignment.csv"))


def _normalize_ss(raw: str) -> str:
    lab = (raw or "").strip().upper()
    if lab in ("H", "E", "C"):
        return lab
    return ""


def _ss_display(code: str) -> str:
    if code in SS_NAME:
        return SS_NAME[code]
    return "unknown"


def _as_bool(value: object) -> bool:
    parsed = parse_optional_bool(value)
    return bool(parsed)


def _is_binding_site(row: Dict[str, str]) -> bool:
    return any(_as_bool(row.get(col, "")) for col in PREDICTOR_COLUMNS)


def _has_finite_ca_distance(row: Dict[str, str]) -> bool:
    raw = (row.get(CA_DISTANCE_COLUMN) or "").strip()
    if not raw:
        return False
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return False
    return value == value  # reject NaN


def count_fps_in_master(
    path: Path,
    *,
    significant_column: str = SIGNIFICANT_COLUMN,
    require_ca_distance: bool = True,
) -> Tuple[Counter, Dict[str, str], Optional[str]]:
    """
    Return (fp_counts_by_ss, meta, skip_reason).

    ``fp_counts_by_ss`` keys are H/E/C/'' (empty = unknown).
    """
    counts: Counter = Counter()
    meta = {"holo_pdb": "", "apo_bmrb": "", "target": path.parent.name}
    try:
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            if SS_COLUMN not in fields:
                return counts, meta, f"missing {SS_COLUMN!r} column (run backfill_dssp_ss.py)"
            if significant_column not in fields:
                return counts, meta, f"missing significance column {significant_column!r}"
            if require_ca_distance and CA_DISTANCE_COLUMN not in fields:
                return counts, meta, f"missing {CA_DISTANCE_COLUMN!r}"
            for row in reader:
                if not meta["holo_pdb"]:
                    meta["holo_pdb"] = (row.get("holo_pdb") or "").strip()
                if not meta["apo_bmrb"]:
                    meta["apo_bmrb"] = (row.get("apo_bmrb") or "").strip()

                if not has_recorded_csp(
                    row,
                    csp_column=CSP_COLUMN,
                    significant_column=significant_column,
                ):
                    continue
                if require_ca_distance and not _has_finite_ca_distance(row):
                    continue

                significant = parse_optional_bool(row.get(significant_column, ""))
                if significant is None or not significant:
                    continue
                if _is_binding_site(row):
                    continue

                counts[_normalize_ss(row.get(SS_COLUMN, ""))] += 1
    except Exception as exc:
        return counts, meta, str(exc)
    return counts, meta, None


def aggregate_targets(
    masters: Sequence[Path],
    *,
    significant_column: str = SIGNIFICANT_COLUMN,
    require_ca_distance: bool = True,
) -> Tuple[Counter, List[Dict[str, object]], List[Tuple[str, str]]]:
    total: Counter = Counter()
    per_target: List[Dict[str, object]] = []
    skipped: List[Tuple[str, str]] = []

    for path in masters:
        counts, meta, reason = count_fps_in_master(
            path,
            significant_column=significant_column,
            require_ca_distance=require_ca_distance,
        )
        if reason:
            skipped.append((path.parent.name, reason))
            continue
        total.update(counts)
        for ss_code in ("H", "E", "C", ""):
            n = int(counts.get(ss_code, 0))
            if n == 0 and ss_code == "":
                continue
            per_target.append(
                {
                    "holo_pdb": meta["holo_pdb"] or path.parent.name.split("_", 1)[0],
                    "apo_bmrb": meta["apo_bmrb"],
                    "target": meta["target"],
                    "ss": ss_code or "",
                    "ss_name": _ss_display(ss_code),
                    "fp_count": n,
                }
            )
    return total, per_target, skipped


def write_summary_csv(path: Path, total: Counter) -> None:
    grand = sum(total.values())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["ss", "ss_name", "fp_count", "fp_fraction"]
        )
        writer.writeheader()
        for ss_code in SS_ORDER:
            n = int(total.get(ss_code, 0))
            if ss_code == "" and n == 0:
                continue
            frac = (n / grand) if grand else 0.0
            writer.writerow(
                {
                    "ss": ss_code,
                    "ss_name": _ss_display(ss_code),
                    "fp_count": n,
                    "fp_fraction": f"{frac:.6f}",
                }
            )
        writer.writerow(
            {
                "ss": "ALL",
                "ss_name": "all",
                "fp_count": grand,
                "fp_fraction": "1.000000" if grand else "0.000000",
            }
        )


def write_per_target_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["holo_pdb", "apo_bmrb", "target", "ss", "ss_name", "fp_count"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _select_masters(
    outputs_dir: Path,
    *,
    targets_csv: Optional[Path],
    all_targets: bool,
) -> List[Path]:
    if all_targets or targets_csv is None:
        return find_master_alignment_files(outputs_dir)
    rows = load_target_rows(targets_csv)
    if not rows:
        return []
    dirs = resolve_target_rows(rows, outputs_dir)
    masters: List[Path] = []
    for d in dirs:
        ma = d / "master_alignment.csv"
        if ma.is_file():
            masters.append(ma)
    return masters


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_ROOT / "outputs",
        help="Root with per-target master_alignment.csv "
        "(default: outputs/)",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=DEFAULT_TARGETS_CSV,
        help="Targets CSV (default: data/CSP_UBQ_ph0.5_temp5C.csv). "
        "Ignored with --all-targets.",
    )
    parser.add_argument(
        "--all-targets",
        action="store_true",
        help="Use every master_alignment.csv under --outputs (ignore --targets-csv).",
    )
    parser.add_argument(
        "--significant-column",
        default=SIGNIFICANT_COLUMN,
        help=f"Significance mask column (default: {SIGNIFICANT_COLUMN})",
    )
    parser.add_argument(
        "--no-require-ca-distance",
        action="store_true",
        help="Do not require a finite min_ca_distance_distance (deviates from Fig 3).",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Aggregate summary CSV path "
        "(default: <outputs>/fp_by_ss_summary.csv)",
    )
    parser.add_argument(
        "--by-target",
        type=Path,
        default=None,
        help="Per-target CSV path "
        "(default: <outputs>/fp_by_ss_by_target.csv)",
    )
    args = parser.parse_args(argv)

    targets_csv = args.targets_csv
    if targets_csv is not None and not targets_csv.is_absolute():
        targets_csv = _ROOT / targets_csv

    masters = _select_masters(
        args.outputs,
        targets_csv=None if args.all_targets else targets_csv,
        all_targets=args.all_targets,
    )
    cohort = "all targets" if args.all_targets else f"targets-csv={targets_csv.name}"
    print(
        f"Found {len(masters)} master_alignment.csv under {args.outputs} ({cohort})",
        flush=True,
    )
    if not masters:
        return 1

    require_ca = not args.no_require_ca_distance
    total, per_target, skipped = aggregate_targets(
        masters,
        significant_column=args.significant_column,
        require_ca_distance=require_ca,
    )

    summary_path = args.summary or (args.outputs / "fp_by_ss_summary.csv")
    by_target_path = args.by_target or (args.outputs / "fp_by_ss_by_target.csv")
    write_summary_csv(summary_path, total)
    write_per_target_csv(by_target_path, per_target)

    helix = int(total.get("H", 0))
    sheet = int(total.get("E", 0))
    loop = int(total.get("C", 0))
    unknown = int(total.get("", 0))
    grand = helix + sheet + loop + unknown

    print(
        f"FP by SS (significance={args.significant_column}, "
        f"require_ca={require_ca}): "
        f"helix={helix}  sheet={sheet}  loop={loop}"
        + (f"  unknown={unknown}" if unknown else "")
        + f"  total={grand}",
        flush=True,
    )
    print(f"Wrote {summary_path}", flush=True)
    print(f"Wrote {by_target_path}", flush=True)

    if skipped:
        print(f"Skipped {len(skipped)} targets:", flush=True)
        for name, reason in skipped[:15]:
            print(f"  {name}: {reason}", flush=True)
        if len(skipped) > 15:
            print(f"  ... and {len(skipped) - 15} more", flush=True)
        if not per_target and grand == 0:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
