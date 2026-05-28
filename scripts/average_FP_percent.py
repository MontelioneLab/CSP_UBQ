#!/usr/bin/env python3
"""
Compute FP-focused summary statistics for all supported CSP families:
N-H, N-H-CA, and HA-CA.

For each receptor and CSP family:
- FP % = (number of FP residues / total residues) * 100
- FP = significant CSP outside the predicted binding site

Also reports the mean of FP/(TP+FP) over receptors with at least one significant
residue (TP+FP > 0), and the mean per-receptor F1 from the same definition as
``f1_score_reporter.py`` (``analyze_targets.compute_f1_score``).

By default only targets listed in the given ``--targets-csv`` are included. Each row is
resolved to ``outputs/{HOLO_PDB}_{apo_bmrb}/`` via :mod:`scripts.target_resolution`
(canonical basename lookup). Set ``CSP_LEGACY_OUTPUT_DIRS=1`` to also match older
``outputs/<pdb>/`` or ``<pdb>_<n>/`` trees using ``master_alignment.csv`` BMRB pairs.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

try:
    from .analyze_targets import (
        AlignmentParsingError,
        compute_f1_score,
        load_alignment,
        PREDICTOR_COLUMNS,
    )
    from .analyze_targets_ca import (
        AlignmentParsingError as AlignmentParsingErrorCA,
        load_ca_alignment,
        to_bool,
    )
    from .target_resolution import load_target_rows, resolve_target_rows
except ImportError:
    from analyze_targets import (
        AlignmentParsingError,
        compute_f1_score,
        load_alignment,
        PREDICTOR_COLUMNS,
    )
    from analyze_targets_ca import (
        AlignmentParsingError as AlignmentParsingErrorCA,
        load_ca_alignment,
        to_bool,
    )
    from target_resolution import load_target_rows, resolve_target_rows


_REQUIRED_TARGET_COLS = frozenset({"apo_bmrb", "holo_bmrb", "holo_pdb"})

MODE_CONFIGS: Dict[str, Dict[str, str]] = {
    "nh": {
        "label": "N-H",
        "csv": "master_alignment.csv",
        "significant": "significant",
    },
    "nh_ca": {
        "label": "N-H-CA",
        "csv": "csp_table_CA.csv",
        "significant": "csp_CA_significant",
    },
    "ha_ca": {
        "label": "HA-CA",
        "csv": "csp_table_HA_CA.csv",
        "significant": "csp_HA_CA_significant",
    },
}


def _compute_f1_from_columns(df: pd.DataFrame, actual: pd.Series, predicted: pd.Series) -> float:
    """Replicate analyze_targets.compute_f1_score without assuming a hardcoded significance column."""
    normalized = df.copy()
    normalized["significant"] = actual
    return compute_f1_score(normalized, predicted).f1


def _validate_target_columns(fieldnames: Optional[List[str]], csv_path: Path) -> None:
    if not fieldnames:
        raise ValueError(f"{csv_path} has no header row")
    missing = _REQUIRED_TARGET_COLS - {fn.lower() for fn in fieldnames}
    if missing:
        raise ValueError(
            f"{csv_path} must contain columns {sorted(_REQUIRED_TARGET_COLS)}; missing {sorted(missing)}"
        )


def _load_mode_dataframe(target_dir: Path, mode_key: str) -> pd.DataFrame:
    cfg = MODE_CONFIGS[mode_key]
    csv_path = target_dir / cfg["csv"]
    if mode_key == "nh":
        return load_alignment(csv_path)

    # N-H-CA / HA-CA rows live in csp_table_*.csv, but binding-site predictor columns come from
    # per-residue overlays (interaction / distance / occlusion CSVs). Reuse analyze_targets_ca
    # merge logic rather than expecting predictors inside the CSP table CSV.
    if mode_key in ("nh_ca", "ha_ca"):
        try:
            df = load_ca_alignment(target_dir, mode_key)
        except AlignmentParsingErrorCA as exc:
            raise AlignmentParsingError(str(exc)) from exc
        significant_col = cfg["significant"]
        if significant_col not in df.columns:
            raise AlignmentParsingError(f"Merged alignment missing significance column {significant_col}")
        df[significant_col] = df[significant_col].apply(to_bool)
        for column in PREDICTOR_COLUMNS:
            if column in df.columns:
                df[column] = df[column].apply(to_bool)
        return df

    raise AssertionError(f"Unhandled mode_key={mode_key!r} — update _load_mode_dataframe after adding CSP modes.")


def _try_append_receptor_metrics(
    target_dir: Path,
    mode_key: str,
    fp_percents: List[float],
    target_to_pct: List[tuple[str, float, int, int]],
    fp_over_tp_fp: List[float],
    f1_scores: List[float],
) -> bool:
    """Parse the mode-specific CSP table and append aggregate metrics."""
    cfg = MODE_CONFIGS[mode_key]
    alignment_path = target_dir / cfg["csv"]
    if not alignment_path.exists():
        return False

    try:
        df = _load_mode_dataframe(target_dir, mode_key)
    except AlignmentParsingError as exc:
        print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
        return False

    significant_col = cfg["significant"]
    if significant_col not in df.columns:
        return False

    available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
    if not available_predictors:
        return False

    predicted = df[available_predictors].any(axis=1)
    actual = df[significant_col]
    total = int(len(df))
    if total == 0:
        return False

    n_fp = int((actual & ~predicted).sum())
    n_tp = int((actual & predicted).sum())
    pct = 100.0 * n_fp / total

    fp_percents.append(pct)
    target_to_pct.append((target_dir.name, pct, n_fp, total))

    sig_denom = n_tp + n_fp
    if sig_denom > 0:
        fp_over_tp_fp.append(n_fp / sig_denom)

    f1_scores.append(_compute_f1_from_columns(df, actual, predicted))
    return True


def _print_mode_summary(
    mode_key: str,
    fp_percents: List[float],
    target_to_pct: List[tuple[str, float, int, int]],
    fp_over_tp_fp: List[float],
    f1_scores: List[float],
) -> None:
    label = MODE_CONFIGS[mode_key]["label"]
    print(f"\n{label} CSPs")
    print("-" * (len(label) + 5))

    if not fp_percents:
        print("No receptors with classification data found.")
        return

    avg_pct = sum(fp_percents) / len(fp_percents)
    min_pct = min(fp_percents)
    min_target = min(target_to_pct, key=lambda x: x[1])

    print(f"Receptors analyzed: {len(fp_percents)}")
    print(f"Average % of FP residues: {avg_pct:.2f}%")
    print(
        f"Minimum % of FP residues: {min_pct:.2f}% "
        f"(receptor: {min_target[0]}, {min_target[2]} FP / {min_target[3]} total)"
    )
    if fp_over_tp_fp:
        mean_fp_ratio = sum(fp_over_tp_fp) / len(fp_over_tp_fp)
        print(f"Mean FP / (TP+FP): {mean_fp_ratio:.3f} (over {len(fp_over_tp_fp)} receptors with TP+FP > 0)")
    else:
        print("Mean FP / (TP+FP): n/a (no receptors with TP+FP > 0)")
    if f1_scores:
        mean_f1 = sum(f1_scores) / len(f1_scores)
        print(f"Mean F1 score: {mean_f1:.3f} (over {len(f1_scores)} receptors)")
    else:
        print("Mean F1 score: n/a")


def discover_targets(outputs_dir: Path) -> list[Path]:
    """Return sorted list of target directories under outputs_dir."""
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    return sorted(p for p in outputs_dir.iterdir() if p.is_dir() and not p.name.startswith("."))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Average and minimum % of FP residues; mean FP/(TP+FP) and mean F1 "
            "for N-H, N-H-CA, and HA-CA CSP analyses."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: outputs)",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help=(
            "CSV with apo_bmrb, holo_bmrb, holo_pdb (one output folder per row). "
            "Default: data/CSP_UBQ_ph0.5_temp5C.csv. Ignored with --all-targets."
        ),
    )
    parser.add_argument(
        "--all-targets",
        action="store_true",
        help="Use every subdirectory under --outputs-dir (ignore --targets-csv).",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else repo_root / args.outputs_dir

    metrics_by_mode = {
        key: {
            "fp_percents": [],
            "target_to_pct": [],
            "fp_over_tp_fp": [],
            "f1_scores": [],
        }
        for key in MODE_CONFIGS
    }

    if args.all_targets:
        print("Using all subdirectories under outputs (--all-targets).")
        target_paths = discover_targets(outputs_dir)
    else:
        targets_csv = args.targets_csv if args.targets_csv.is_absolute() else repo_root / args.targets_csv
        if not targets_csv.is_file():
            print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
            return 1
        try:
            with targets_csv.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                _validate_target_columns(reader.fieldnames, targets_csv)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        rows = load_target_rows(targets_csv)
        target_paths = resolve_target_rows(rows, outputs_dir)
        print(
            f"Resolving {len(rows)} rows from {targets_csv.name} to {len(target_paths)} output dirs "
            "(BMRB-congruent master_alignment match; see scripts/target_resolution)"
        )

    for target_dir in target_paths:
        for mode_key, mode_metrics in metrics_by_mode.items():
            _try_append_receptor_metrics(
                target_dir,
                mode_key,
                mode_metrics["fp_percents"],
                mode_metrics["target_to_pct"],
                mode_metrics["fp_over_tp_fp"],
                mode_metrics["f1_scores"],
            )

    if not any(metrics["fp_percents"] for metrics in metrics_by_mode.values()):
        print("No receptors with classification data found.")
        return 0

    for mode_key, mode_metrics in metrics_by_mode.items():
        _print_mode_summary(
            mode_key,
            mode_metrics["fp_percents"],
            mode_metrics["target_to_pct"],
            mode_metrics["fp_over_tp_fp"],
            mode_metrics["f1_scores"],
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
