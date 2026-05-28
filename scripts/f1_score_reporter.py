#!/usr/bin/env python3
"""
F1 Score Reporter - Calculate and report F1 scores for user-specified targets
across N-H, N-H-CA, and HA-CA CSP analyses.

Provide either:

- ``--targets-csv``: ``holo_pdb`` column defines the target set (no extra validation).
- ``--targets``: comma-separated holo PDB IDs; each must appear in ``data/CSP_UBQ.csv``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Set

import pandas as pd

DEFAULT_VALIDATION_CSV = Path("data/CSP_UBQ.csv")

try:
    from .analyze_targets import (
        AlignmentParsingError,
        TargetResult,
        load_alignment,
        compute_f1_score,
        PREDICTOR_COLUMNS,
    )
    from .analyze_targets_ca import to_bool
    from .target_resolution import TargetRow, load_target_rows, resolve_target_rows
except ImportError:
    from analyze_targets import (
        AlignmentParsingError,
        TargetResult,
        load_alignment,
        compute_f1_score,
        PREDICTOR_COLUMNS,
    )
    from analyze_targets_ca import to_bool
    from target_resolution import TargetRow, load_target_rows, resolve_target_rows

MODE_CONFIGS: Dict[str, Dict[str, str]] = {
    "nh": {"label": "N-H", "csv": "master_alignment.csv", "significant": "significant"},
    "nh_ca": {"label": "N-H-CA", "csv": "csp_table_CA.csv", "significant": "csp_CA_significant"},
    "ha_ca": {"label": "HA-CA", "csv": "csp_table_HA_CA.csv", "significant": "csp_HA_CA_significant"},
}


def _compute_f1_from_columns(df: pd.DataFrame, actual: pd.Series, predicted: pd.Series) -> TargetResult:
    normalized = df.copy()
    normalized["significant"] = actual
    return compute_f1_score(normalized, predicted)


def load_holo_pdb_set(targets_csv: Path) -> Set[str]:
    """Lowercase holo_pdb IDs from CSV (holo_pdb column)."""
    out: Set[str] = set()
    with targets_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "holo_pdb" not in reader.fieldnames:
            raise ValueError(f"{targets_csv} must contain a 'holo_pdb' column")
        for row in reader:
            holo_pdb = (row.get("holo_pdb") or "").strip().lower()
            if holo_pdb:
                out.add(holo_pdb)
    return out


def _filter_rows_by_holo(rows: Sequence[TargetRow], holo_targets: Set[str]) -> List[TargetRow]:
    wanted = {h.strip().lower() for h in holo_targets if h.strip()}
    return [row for row in rows if row.holo_pdb.strip().lower() in wanted]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate and report F1 scores for a subset of targets across "
            "N-H, N-H-CA, and HA-CA CSP analyses."
        )
    )
    target_grp = parser.add_mutually_exclusive_group(required=True)
    target_grp.add_argument(
        "--targets",
        type=str,
        help="Comma-separated list of holo PDB IDs (e.g., '1cf4,2lox,2mkr')",
    )
    target_grp.add_argument(
        "--targets-csv",
        type=Path,
        help="CSV with a holo_pdb column defining targets (like average_FP_percent.py).",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=sorted(MODE_CONFIGS),
        default=["nh", "nh_ca", "ha_ca"],
        help="Which CSP families to report (default: %(default)s).",
    )
    return parser.parse_args(argv)


def validate_targets(target_list: List[str], validation_csv: Path) -> tuple[Set[str], List[str]]:
    if not validation_csv.exists():
        raise FileNotFoundError(f"Validation CSV not found: {validation_csv}")

    try:
        df = pd.read_csv(validation_csv)
        if "holo_pdb" not in df.columns:
            raise ValueError(f"{validation_csv} must have a 'holo_pdb' column")

        valid_targets = set(df["holo_pdb"].astype(str).str.strip().str.lower())
        normalized_input = [t.strip().lower() for t in target_list]
        invalid_targets = [t for t in normalized_input if t not in valid_targets]
        valid_set = {t for t in normalized_input if t in valid_targets}
        return valid_set, invalid_targets
    except Exception as exc:
        raise ValueError(f"Failed to read {validation_csv}: {exc}") from exc


def _load_mode_dataframe(target_dir: Path, mode_key: str) -> pd.DataFrame:
    cfg = MODE_CONFIGS[mode_key]
    csv_path = target_dir / cfg["csv"]
    if mode_key == "nh":
        return load_alignment(csv_path)

    df = pd.read_csv(csv_path)
    significant_col = cfg["significant"]
    required_columns = (significant_col, *PREDICTOR_COLUMNS)
    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise AlignmentParsingError(
            f"Alignment file {csv_path} is missing required columns: {', '.join(missing_columns)}"
        )
    for column in (significant_col, *PREDICTOR_COLUMNS):
        if column in df.columns:
            df[column] = df[column].apply(to_bool)
    return df


def _resolve_target_paths_from_rows(
    rows: Sequence[TargetRow],
    outputs_dir: Path,
    validation_rows: Sequence[TargetRow],
) -> List[Path]:
    """Resolve requested target rows to canonical outputs dirs, backfilling apo_bmrb when needed."""
    direct_rows: List[TargetRow] = []
    missing_apo_holos: Set[str] = set()
    for row in rows:
        if row.apo_bmrb.strip():
            direct_rows.append(row)
        else:
            missing_apo_holos.add(row.holo_pdb.strip().lower())

    resolved = resolve_target_rows(direct_rows, outputs_dir, log_warnings=False) if direct_rows else []
    if missing_apo_holos:
        fallback_rows = _filter_rows_by_holo(validation_rows, missing_apo_holos)
        resolved.extend(resolve_target_rows(fallback_rows, outputs_dir, log_warnings=False))

    deduped: List[Path] = []
    seen: Set[str] = set()
    for path in resolved:
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        deduped.append(path)
    return deduped


def calculate_f1_for_targets(target_paths: Sequence[Path], mode_key: str) -> List[TargetResult]:
    results: List[TargetResult] = []

    for target_dir in sorted(target_paths, key=lambda p: p.name.lower()):
        alignment_path = target_dir / MODE_CONFIGS[mode_key]["csv"]
        if not alignment_path.exists():
            print(f"[WARN] {alignment_path.name} not found for {target_dir.name}. Skipping.", file=sys.stderr)
            continue

        try:
            df = _load_mode_dataframe(target_dir, mode_key)
            available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
            if not available_predictors:
                raise AlignmentParsingError(f"No predictor columns found in {alignment_path}")
            predicted = df[available_predictors].any(axis=1)
            metrics = _compute_f1_from_columns(df, df[MODE_CONFIGS[mode_key]["significant"]], predicted)
            metrics.target = target_dir.name
            results.append(metrics)
        except AlignmentParsingError as exc:
            print(f"[WARN] Failed to parse {alignment_path}: {exc}. Skipping.", file=sys.stderr)
            continue
        except Exception as exc:
            print(f"[WARN] Unexpected error processing {target_dir.name}: {exc}. Skipping.", file=sys.stderr)
            continue

    return results


def print_results_table(results: List[TargetResult], mode_key: str) -> None:
    if not results:
        print("No results to display.")
        return

    df = pd.DataFrame([result.__dict__ for result in results]).sort_values("f1", ascending=False)

    print("\n" + "=" * 70)
    print(f"F1 Score Report - {MODE_CONFIGS[mode_key]['label']} CSPs")
    print("=" * 70)
    print("\nTarget Results:")
    print("-" * 70)

    display_df = df[["target", "f1", "true_positives", "false_positives", "false_negatives", "total_rows"]].copy()
    display_df.columns = ["Target", "F1 Score", "TP", "FP", "FN", "Total Rows"]
    display_df["F1 Score"] = display_df["F1 Score"].map("{:.3f}".format)

    print(display_df.to_string(index=False))
    print("-" * 70)

    f1_scores = df["f1"].values
    print("\nSummary Statistics:")
    print("-" * 70)
    print(f"  Count:  {len(results)}")
    print(f"  Mean:   {f1_scores.mean():.3f}")
    print(f"  Median: {pd.Series(f1_scores).median():.3f}")
    print(f"  Min:    {f1_scores.min():.3f}")
    print(f"  Max:    {f1_scores.max():.3f}")
    print("=" * 70 + "\n")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parents[1]
    validation_path = DEFAULT_VALIDATION_CSV if DEFAULT_VALIDATION_CSV.is_absolute() else repo_root / DEFAULT_VALIDATION_CSV
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else repo_root / args.outputs_dir
    validation_rows = load_target_rows(validation_path)

    target_paths: List[Path]

    if args.targets_csv is not None:
        targets_csv = args.targets_csv if args.targets_csv.is_absolute() else repo_root / args.targets_csv
        if not targets_csv.is_file():
            print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
            return 1
        try:
            requested_rows = load_target_rows(targets_csv)
            requested_holos = load_holo_pdb_set(targets_csv)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if not requested_holos:
            print("Error: No holo_pdb values in targets CSV.", file=sys.stderr)
            return 1
        target_paths = _resolve_target_paths_from_rows(requested_rows, outputs_dir, validation_rows)
        print(f"Filtering to {len(requested_holos)} holo_pdb IDs from {targets_csv.name}")
    else:
        assert args.targets is not None
        target_list = [t.strip() for t in args.targets.split(",") if t.strip()]
        if not target_list:
            print("Error: No targets provided. Use --targets with comma-separated PDB IDs.", file=sys.stderr)
            return 1
        try:
            valid_targets, invalid_targets = validate_targets(target_list, validation_path)
        except Exception as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if invalid_targets:
            print(
                f"Warning: The following targets were not found in {validation_path}: {', '.join(invalid_targets)}",
                file=sys.stderr,
            )
        if not valid_targets:
            print("Error: No valid targets found. Exiting.", file=sys.stderr)
            return 1
        requested_rows = _filter_rows_by_holo(validation_rows, valid_targets)
        target_paths = _resolve_target_paths_from_rows(requested_rows, outputs_dir, validation_rows)

    if not outputs_dir.exists():
        print(f"Error: Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1
    if not target_paths:
        print("Error: No matching output directories were found for the requested targets.", file=sys.stderr)
        return 1

    any_results = False
    for mode_key in args.modes:
        results = calculate_f1_for_targets(target_paths, mode_key)
        if not results:
            print(
                f"Warning: No results calculated for {MODE_CONFIGS[mode_key]['label']} "
                f"({MODE_CONFIGS[mode_key]['csv']}).",
                file=sys.stderr,
            )
            continue
        any_results = True
        print_results_table(results, mode_key)

    if not any_results:
        print("Error: No results calculated for any requested CSP family.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
