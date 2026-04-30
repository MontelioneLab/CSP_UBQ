#!/usr/bin/env python3
"""
Sweep buffer-difference tolerances and aggregate per-target classification metrics.

For every (pH allowance, temperature allowance) combination on a 2D grid, identify
the subset of CSP_UBQ targets whose paired apo/holo buffer rows satisfy
|ΔpH| <= ph_allowance AND |ΔT| <= temp_allowance, then aggregate the same four
statistics ``average_FP_percent.py`` reports:

  - n              = number of targets in the subset (with classifications)
  - mean_F1        = mean per-target F1 (same definition as analyze_targets.py)
  - pct_allosteric = mean over targets with TP+FP > 0 of FP / (TP+FP)
  - fp_pct         = mean over subset of 100 * FP / (TP+FP+TN+FN)

Per-target metrics are computed once via a single pass over ``outputs/`` and the
resolution logic from ``average_FP_percent.py`` is replayed cell-by-cell so that
duplicate-suffix bookkeeping matches its default behaviour exactly. Writes one
summary CSV plus four heatmap PNGs into ``outputs/buffer_threshold_sweep/``
(configurable).
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.analyze_targets import (  # noqa: E402
    AlignmentParsingError,
    PREDICTOR_COLUMNS,
    compute_f1_score,
    load_alignment,
)
from scripts.average_FP_percent import (  # noqa: E402
    CLASSIFICATION_COLUMN,
    VALID_CLASSIFICATIONS,
    _build_holo_pdb_groups,
    _build_outputs_index,
    _candidate_dirs_for_holo,
    _logical_output_dirname,
)
from scripts.config import Paths  # noqa: E402
from scripts.filter_csp_ubq_by_buffer import row_meets  # noqa: E402


@dataclass(frozen=True)
class TargetMetrics:
    n_fp: int
    n_tp: int
    n_total: int
    f1: float


def _parse_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _exp_diffs(exp_row: Dict[str, str]) -> Optional[Tuple[float, float]]:
    """Return (|ΔpH|, |ΔT|) if all four values are present and numeric, else None."""
    aph = _parse_float(exp_row.get("apo_pH"))
    hph = _parse_float(exp_row.get("holo_pH"))
    at = _parse_float(exp_row.get("apo_temperature_C"))
    ht = _parse_float(exp_row.get("holo_temperature_C"))
    if aph is None or hph is None or at is None or ht is None:
        return None
    return abs(aph - hph), abs(at - ht)


def _per_target_metrics(target_dir: Path) -> Optional[TargetMetrics]:
    """Return (n_fp, n_tp, n_total, f1) for a target dir; None if unusable."""
    alignment_path = target_dir / "master_alignment.csv"
    if not alignment_path.is_file():
        return None
    try:
        df = load_alignment(alignment_path)
    except AlignmentParsingError as exc:
        print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
        return None
    if CLASSIFICATION_COLUMN not in df.columns:
        return None
    cls = df[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper()
    valid_mask = cls.isin(VALID_CLASSIFICATIONS)
    classified = df.loc[valid_mask]
    total = len(classified)
    if total == 0:
        return None
    cls_upper = classified[CLASSIFICATION_COLUMN].astype(str).str.strip().str.upper()
    n_fp = int((cls_upper == "FP").sum())
    n_tp = int((cls_upper == "TP").sum())
    predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)
    f1 = compute_f1_score(df, predicted).f1
    return TargetMetrics(n_fp=n_fp, n_tp=n_tp, n_total=total, f1=f1)


def _first_row_bmrb_pair_from_df(alignment_path: Path) -> Optional[Tuple[str, str]]:
    """First non-empty (apo_bmrb, holo_bmrb) in a master_alignment.csv (cached use)."""
    if not alignment_path.is_file():
        return None
    with alignment_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a = (row.get("apo_bmrb") or "").strip()
            h = (row.get("holo_bmrb") or "").strip()
            if a or h:
                return a, h
    return None


def _build_outputs_caches(
    outputs_dir: Path,
) -> Tuple[
    Dict[str, Path],
    Dict[Path, Optional[TargetMetrics]],
    Dict[Path, Optional[Tuple[str, str]]],
]:
    """Pre-compute (1) outputs_index, (2) per-dir metrics, (3) per-dir first BMRB pair."""
    outputs_index = _build_outputs_index(outputs_dir)
    metrics_cache: Dict[Path, Optional[TargetMetrics]] = {}
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]] = {}
    for path in outputs_index.values():
        metrics_cache[path] = _per_target_metrics(path)
        bmrb_cache[path] = _first_row_bmrb_pair_from_df(path / "master_alignment.csv")
    return outputs_index, metrics_cache, bmrb_cache


def _resolve_with_caches(
    outputs_index: Dict[str, Path],
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]],
    logical: str,
    row: Dict[str, str],
) -> Optional[Path]:
    """Same logic as ``average_FP_percent._resolve_output_dir_path`` using cached BMRB pairs."""
    direct = outputs_index.get(logical.lower())
    if direct is not None:
        return direct
    raw_h = (row.get("holo_pdb") or "").strip()
    apo = (row.get("apo_bmrb") or "").strip()
    holo_b = (row.get("holo_bmrb") or "").strip()
    want = (apo, holo_b)
    candidates = _candidate_dirs_for_holo(outputs_index, raw_h)
    matches = [p for p in candidates if bmrb_cache.get(p) == want]
    if len(matches) == 1:
        return matches[0]
    return None


def load_paired_rows(
    csp_path: Path,
    exp_path: Path,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], List[Optional[Tuple[float, float]]]]:
    """Load CSP and EXP rows paired by row index, plus precomputed (|ΔpH|, |ΔT|)."""
    with csp_path.open(newline="", encoding="utf-8") as f:
        csp_rows = list(csv.DictReader(f))
    with exp_path.open(newline="", encoding="utf-8") as f:
        exp_rows = list(csv.DictReader(f))
    n_pairs = min(len(csp_rows), len(exp_rows))
    if len(csp_rows) != len(exp_rows):
        print(
            f"WARN: row counts differ (csp={len(csp_rows)}, exp={len(exp_rows)}); "
            f"using first {n_pairs} pairs.",
            file=sys.stderr,
        )
    csp_rows = csp_rows[:n_pairs]
    exp_rows = exp_rows[:n_pairs]
    diffs = [_exp_diffs(er) for er in exp_rows]
    return csp_rows, exp_rows, diffs


def aggregate_grid(
    csp_rows: List[Dict[str, str]],
    exp_rows: List[Dict[str, str]],
    diffs: List[Optional[Tuple[float, float]]],
    outputs_index: Dict[str, Path],
    metrics_cache: Dict[Path, Optional[TargetMetrics]],
    bmrb_cache: Dict[Path, Optional[Tuple[str, str]]],
    ph_values: np.ndarray,
    temp_values: np.ndarray,
) -> pd.DataFrame:
    """Replay average_FP_percent's resolution + aggregation per (T, pH) cell."""
    records: List[dict] = []
    for t in temp_values:
        for p in ph_values:
            ph_max = float(p)
            t_max = float(t)
            kept_indices = [
                i for i, d in enumerate(diffs)
                if d is not None and d[0] <= ph_max and d[1] <= t_max
            ]
            subset_rows = [csp_rows[i] for i in kept_indices]
            holo_groups = _build_holo_pdb_groups(subset_rows)

            seen_paths: set[str] = set()
            metrics: List[TargetMetrics] = []
            for row in subset_rows:
                logical = _logical_output_dirname(row, holo_groups)
                if not logical:
                    continue
                path = _resolve_with_caches(outputs_index, bmrb_cache, logical, row)
                if path is None:
                    continue
                key = str(path.resolve())
                if key in seen_paths:
                    continue
                seen_paths.add(key)
                m = metrics_cache.get(path)
                if m is None:
                    continue
                metrics.append(m)

            n = len(metrics)
            if n == 0:
                records.append(
                    {
                        "temp_allowance": t_max,
                        "ph_allowance": ph_max,
                        "n": 0,
                        "mean_F1": np.nan,
                        "pct_allosteric": np.nan,
                        "fp_pct": np.nan,
                    }
                )
                continue
            f1_arr = np.array([m.f1 for m in metrics])
            n_fp_arr = np.array([m.n_fp for m in metrics])
            n_tp_arr = np.array([m.n_tp for m in metrics])
            n_total_arr = np.array([m.n_total for m in metrics])
            sig_denom = n_tp_arr + n_fp_arr
            mean_f1 = float(np.mean(f1_arr))
            fp_pct = float(np.mean(100.0 * n_fp_arr / n_total_arr))
            sig_mask = sig_denom > 0
            if sig_mask.any():
                pct_allosteric = float(np.mean(n_fp_arr[sig_mask] / sig_denom[sig_mask]))
            else:
                pct_allosteric = float("nan")
            records.append(
                {
                    "temp_allowance": t_max,
                    "ph_allowance": ph_max,
                    "n": n,
                    "mean_F1": mean_f1,
                    "pct_allosteric": pct_allosteric,
                    "fp_pct": fp_pct,
                }
            )
    return pd.DataFrame.from_records(records)


def _save_heatmap(
    df: pd.DataFrame,
    value_col: str,
    out_path: Path,
    title: str,
    cbar_label: str,
    fmt: str,
    cmap: str,
    annotate: bool,
) -> None:
    pivot = df.pivot(index="ph_allowance", columns="temp_allowance", values=value_col)
    pivot = pivot.sort_index(ascending=False)
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)

    fig_w = max(8.0, 0.5 * pivot.shape[1] + 3.5)
    fig_h = max(6.0, 0.35 * pivot.shape[0] + 2.5)
    plt.figure(figsize=(fig_w, fig_h))

    if value_col == "n":
        plot_data = pivot.fillna(0).astype(int)
        annot_data = plot_data if annotate else False
    else:
        plot_data = pivot
        annot_data = annotate

    ax = sns.heatmap(
        plot_data,
        annot=annot_data,
        fmt=fmt if annotate else "",
        cmap=cmap,
        cbar_kws={"label": cbar_label},
        linewidths=0.3,
        linecolor="white",
        square=False,
    )
    ax.set_xlabel("Temperature allowance (°C)")
    ax.set_ylabel("pH allowance")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> int:
    cfg = Paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csp",
        type=Path,
        default=Path(cfg.input_csv),
        help="Input CSP_UBQ CSV (default: data/CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--exp",
        type=Path,
        default=Path(cfg.exp_conditions_csv),
        help="apo_holo_exp_conditions CSV (default: data/apo_holo_exp_conditions.csv).",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: outputs).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs/buffer_threshold_sweep"),
        help="Where to write sweep_metrics.csv and the heatmap PNGs.",
    )
    parser.add_argument("--ph-min", type=float, default=0.1)
    parser.add_argument("--ph-max", type=float, default=2.0)
    parser.add_argument("--ph-step", type=float, default=0.1)
    parser.add_argument("--temp-min", type=float, default=2.0)
    parser.add_argument("--temp-max", type=float, default=20.0)
    parser.add_argument("--temp-step", type=float, default=1.0)
    parser.add_argument(
        "--no-annot",
        action="store_true",
        help="Skip per-cell numeric labels in heatmaps.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    csp_path = args.csp if args.csp.is_absolute() else repo_root / args.csp
    exp_path = args.exp if args.exp.is_absolute() else repo_root / args.exp
    outputs_dir = (
        args.outputs_dir if args.outputs_dir.is_absolute() else repo_root / args.outputs_dir
    )
    out_dir = args.out_dir if args.out_dir.is_absolute() else repo_root / args.out_dir

    for label, p in (("--csp", csp_path), ("--exp", exp_path)):
        if not p.is_file():
            print(f"Error: {label} file not found: {p}", file=sys.stderr)
            return 1
    if not outputs_dir.is_dir():
        print(f"Error: --outputs-dir not found: {outputs_dir}", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    ph_values = np.round(
        np.arange(args.ph_min, args.ph_max + args.ph_step / 2.0, args.ph_step), 6
    )
    temp_values = np.round(
        np.arange(args.temp_min, args.temp_max + args.temp_step / 2.0, args.temp_step),
        6,
    )

    print(f"Pre-caching outputs/ metrics...")
    outputs_index, metrics_cache, bmrb_cache = _build_outputs_caches(outputs_dir)
    n_with_metrics = sum(1 for m in metrics_cache.values() if m is not None)
    print(
        f"Indexed {len(outputs_index)} outputs subdirs "
        f"({n_with_metrics} with classifications)."
    )

    print(f"Loading paired rows from {csp_path.name} and {exp_path.name}...")
    csp_rows, exp_rows, diffs = load_paired_rows(csp_path, exp_path)
    n_with_diffs = sum(1 for d in diffs if d is not None)
    print(f"Paired rows: {len(csp_rows)} ({n_with_diffs} with complete buffer metadata).")

    print(
        f"Aggregating {len(ph_values)} pH x {len(temp_values)} T "
        f"= {len(ph_values) * len(temp_values)} cells..."
    )
    df = aggregate_grid(
        csp_rows,
        exp_rows,
        diffs,
        outputs_index,
        metrics_cache,
        bmrb_cache,
        ph_values,
        temp_values,
    )
    df = df.sort_values(["temp_allowance", "ph_allowance"]).reset_index(drop=True)

    csv_path = out_dir / "sweep_metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} ({len(df)} rows).")

    annotate = not args.no_annot
    _save_heatmap(
        df,
        value_col="n",
        out_path=out_dir / "heatmap_n.png",
        title="Number of targets in subset",
        cbar_label="n targets",
        fmt="d",
        cmap="viridis",
        annotate=annotate,
    )
    _save_heatmap(
        df,
        value_col="mean_F1",
        out_path=out_dir / "heatmap_mean_f1.png",
        title="Mean per-target F1",
        cbar_label="mean F1",
        fmt=".3f",
        cmap="viridis",
        annotate=annotate,
    )
    _save_heatmap(
        df,
        value_col="pct_allosteric",
        out_path=out_dir / "heatmap_pct_allosteric.png",
        title="Mean FP / (TP + FP)  (% allosteric CSPs)",
        cbar_label="FP / (TP + FP)",
        fmt=".3f",
        cmap="rocket_r",
        annotate=annotate,
    )
    _save_heatmap(
        df,
        value_col="fp_pct",
        out_path=out_dir / "heatmap_fp_pct.png",
        title="Mean FP / (TP + FP + TN + FN)  (FP %)",
        cbar_label="FP % of classified residues",
        fmt=".1f",
        cmap="rocket_r",
        annotate=annotate,
    )
    print(f"Wrote 4 heatmap PNGs to {out_dir}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
