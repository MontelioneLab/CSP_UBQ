#!/usr/bin/env python3
"""
Per-target single-atom (H, N, CA, HA) 1D shift perturbation analysis.

For each target directory under outputs/ that has a CSP table, this script:
  - Loads unreferenced apo/holo shifts (H/N from csp_table.csv, CA from
    csp_table_CA.csv, HA from csp_table_HA_CA.csv)
  - Runs an independent 1D offset grid for each nucleus and writes
    offset_grid_1d_*.csv
  - Computes per-residue 1D shift metrics:
      CSP_X_1d = |(holo_original + offset) - apo|
      significance when CSP_X_1d >= max(cleaned mean, 0.05 ppm)
  - Writes outputs/{target}/1d_analysis.csv with one row per residue
    containing all single-atom 1D metrics for that residue.

Shared helpers :func:`fraction_rows_with_both_ca_shifts_1d` and
:func:`target_basenames_passing_ca_shift_coverage` define the CA-shift coverage
rule for SI Fig. S10 / S12 / S15 and the standalone 1D F1 boxplot (default
``DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE``: strictly more than half of rows with both
``CA_apo`` and ``CA_holo``).

No scaling coefficients are applied to the 1D CSPs; they are plain absolute
differences in ppm for each atom type.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Set, Tuple

import statistics
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Reuse the canonical outlier-removal significance-threshold routine that the
# 2D/HN-weighted pipeline in scripts/csp.py uses (outlier_z=3 iteratively,
# threshold = mean + significance_z * SD of the cleaned subset, defaults
# yielding threshold = mean of outlier-free subset).
try:
    from .csp import (  # type: ignore
        _floor_primary_hn_cutoff,
        _fmt_grid_num,
        compute_threshold_with_outlier_removal,
        run_offset_grid_search_1d,
    )
    from .config import Referencing as _Referencing  # type: ignore
    from .config import thresholds as _thresholds  # type: ignore
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(os.path.abspath(__file__))))
    from scripts.csp import (  # type: ignore
        _floor_primary_hn_cutoff,
        _fmt_grid_num,
        compute_threshold_with_outlier_removal,
        run_offset_grid_search_1d,
    )
    from scripts.config import Referencing as _Referencing  # type: ignore
    from scripts.config import thresholds as _thresholds  # type: ignore


# Default cutoff for SI figures pairing CA-inclusive analyses with the 1D table: strictly greater than
# this fraction of ``1d_analysis.csv`` rows must have both CA_apo and CA_holo populated.
DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE = 0.5


@dataclass
class Atom1DStats:
    delta: Optional[float]
    csp_1d: Optional[float]
    z_1d: Optional[float]
    significant: Optional[bool] = None


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-residue 1D single-atom shift perturbations (H, N, CA, HA) "
            "from csp_table.csv and csp_table_CA.csv and write 1d_analysis.csv "
            "for each target, then compute per-target F1 scores for each atom type "
            "and render summary heatmaps/statistics."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        help=(
            "Optional CSV file containing target IDs to filter. "
            "Should have a 'holo_pdb' column. Only targets listed in this file will be included."
        ),
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("outputs") / "summary_statistics" / "f1_heatmap_1d_single_atom.png",
        help="Destination for the combined H/N/CA 1D F1 heatmap figure (default: %(default)s).",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("outputs") / "summary_statistics" / "f1_summary_1d_single_atom.csv",
        help="Destination for the summary statistics CSV of 1D F1 scores (default: %(default)s).",
    )
    parser.add_argument(
        "--boxplot-image",
        type=Path,
        help=(
            "Destination for the F1 score boxplot figure. "
            "If not provided, defaults to the same directory as --summary-csv with name 'f1_boxplot_1d_single_atom.png'."
        ),
    )
    return parser.parse_args(list(argv))


def discover_targets(outputs_dir: Path) -> List[Path]:
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")
    return sorted(
        path for path in outputs_dir.iterdir() if path.is_dir() and not path.name.startswith(".")
    )


def fraction_rows_with_both_ca_shifts_1d(target_dir: Path) -> Optional[float]:
    """Return the fraction of rows in ``1d_analysis.csv`` with both ``CA_apo`` and ``CA_holo``.

    Uses the residue row count from ``1d_analysis.csv``. Returns ``None`` if the file is missing,
    unreadable, empty, or lacks the CA columns.
    """
    one_d_path = target_dir / "1d_analysis.csv"
    if not one_d_path.exists():
        return None
    try:
        df = pd.read_csv(one_d_path)
    except Exception:
        return None
    if df.empty or "CA_apo" not in df.columns or "CA_holo" not in df.columns:
        return None
    total = int(len(df))
    if total == 0:
        return None
    both = int((df["CA_apo"].notna() & df["CA_holo"].notna()).sum())
    return both / total


def target_basenames_passing_ca_shift_coverage(
    outputs_dir: Path,
    *,
    min_coverage: float = DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
    allowed_basenames: Optional[Mapping[str, object]] = None,
) -> Tuple[Set[str], Dict[str, float]]:
    """Per-target dirs whose CA row fraction in ``1d_analysis.csv`` is **strictly** ``> min_coverage``.

    Shared eligibility rule for SI Fig. S10 / S12 / S15 and the standalone 1D F1 boxplot so CA-shift
    targets are gated identically from the pipeline 1D table.

    Args:
        outputs_dir: Pipeline ``outputs/`` root.
        min_coverage: Exclusive lower bound (default 0.5 → strictly more than half the rows).
        allowed_basenames: If set, only these subdirectory names are scanned.

    Returns:
        ``passing``: basenames above the cutoff.
        ``coverage_by_basename``: all scanned targets that had a readable finite coverage fraction.
    """
    passing: Set[str] = set()
    coverage_by_basename: Dict[str, float] = {}
    if not outputs_dir.exists():
        return passing, coverage_by_basename

    outs = outputs_dir.resolve()
    min_c = float(min_coverage)
    for path in sorted(outs.iterdir()):
        if not path.is_dir() or path.name.startswith("."):
            continue
        if allowed_basenames is not None and path.name not in allowed_basenames:
            continue
        frac = fraction_rows_with_both_ca_shifts_1d(path)
        if frac is None:
            continue
        coverage_by_basename[path.name] = float(frac)
        if frac > min_c:
            passing.add(path.name)
    return passing, coverage_by_basename


def load_allowed_targets(targets_csv: Optional[Path]) -> Optional[Dict[str, bool]]:
    if not targets_csv:
        return None
    if not targets_csv.exists():
        raise FileNotFoundError(f"Targets CSV not found: {targets_csv}")
    allowed: Dict[str, bool] = {}
    import pandas as pd

    df = pd.read_csv(targets_csv)
    if "holo_pdb" not in df.columns:
        raise ValueError(f"Targets CSV {targets_csv} must have a 'holo_pdb' column")
    for val in df["holo_pdb"].astype(str).str.strip():
        if val:
            allowed[val] = True
    return allowed


# Predictor columns mirror those used in analyze_targets.py
PREDICTOR_COLUMNS: Tuple[str, ...] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


@dataclass
class TargetResult:
    target: str
    atom_type: str
    f1: float
    true_positives: int
    false_positives: int
    false_negatives: int
    total_rows: int


def compute_atom_stats(values: List[float]) -> Tuple[float, float]:
    """Return (mean, sd) for a list of values; sd is population SD."""
    if not values:
        return 0.0, 0.0
    if len(values) == 1:
        return float(values[0]), 0.0
    mean_val = statistics.mean(values)
    sd_val = statistics.pstdev(values)
    return float(mean_val), float(sd_val)


def _merge_1d_with_alignment(target_dir: Path) -> Optional[pd.DataFrame]:
    """
    Load 1d_analysis.csv and master_alignment.csv for a target and merge them on (resi, aa).
    Returns a merged DataFrame or None if required files are missing.
    """
    one_d_path = target_dir / "1d_analysis.csv"
    align_path = target_dir / "master_alignment.csv"
    target_name = target_dir.name
    
    if not one_d_path.exists() or not align_path.exists():
        return None

    try:
        df_1d = pd.read_csv(one_d_path)
        df_align = pd.read_csv(align_path)
    except Exception:
        return None

    # Ensure merge keys exist in 1d_analysis.csv
    if "resi" not in df_1d.columns or "aa" not in df_1d.columns:
        return None

    # Handle master_alignment.csv: it may have 'resi'/'aa' or 'holo_resi'/'holo_aa'
    if "resi" not in df_align.columns or "aa" not in df_align.columns:
        # Try using holo_resi/holo_aa and rename them
        if "holo_resi" in df_align.columns and "holo_aa" in df_align.columns:
            df_align = df_align.copy()
            df_align["resi"] = df_align["holo_resi"]
            df_align["aa"] = df_align["holo_aa"]
        else:
            return None

    def _normalize(df: pd.DataFrame, resi_col: str = "resi", aa_col: str = "aa") -> pd.DataFrame:
        df = df.copy()
        df[resi_col] = pd.to_numeric(df[resi_col], errors="coerce")
        df = df[df[resi_col].notna()]
        df[resi_col] = df[resi_col].astype(int)
        df[aa_col] = df[aa_col].astype(str).str.strip()
        return df

    df_1d = _normalize(df_1d)
    df_align = _normalize(df_align)

    merged = pd.merge(
        df_1d,
        df_align,
        on=["resi", "aa"],
        how="inner",
        suffixes=("_1d", "_align"),
    )
    if merged.empty:
        return None
    
    return merged


def _compute_f1_for_atom(
    df: pd.DataFrame,
    atom_type: str,
    target_name: str = "",
) -> Optional[TargetResult]:
    """
    Compute F1 score for a single atom type using the 2D-convention significance
    column (csp_X_1d_significant) as ground truth and predictor columns as the
    prediction.

    The significance column is computed per target via iterative outlier removal
    (mirroring scripts/csp.py for the HN-weighted pipeline): outliers with
    z > outlier_z are removed iteratively, then the cutoff for significance
    is ``max(cleaned mean, 0.05 ppm)``.

    For backward compatibility, if the significance column is missing we fall
    back to the legacy ``z_X_1d > 0`` rule.
    """
    # After merge, columns may carry a '_1d' or '_align' suffix. Probe several
    # naming patterns and use the first candidate that has any non-null data.
    def _pick_column(base: str) -> Optional[str]:
        for candidate in [f"{base}_1d", f"{base}_align", base]:
            if candidate in df.columns and df[candidate].notna().any():
                return candidate
        return None

    sig_col = _pick_column(f"csp_{atom_type}_1d_significant")
    z_col = _pick_column(f"z_{atom_type}_1d")

    if sig_col is not None:
        valid_mask = df[sig_col].notna()
        if not valid_mask.any():
            return None
        raw = df.loc[valid_mask, sig_col]
        # Column is stored as "True"/"False" strings in 1d_analysis.csv but may
        # also arrive as bool/numeric depending on loader. Normalise robustly.
        if raw.dtype == bool:
            actual = raw.astype(bool)
        else:
            actual = raw.astype(str).str.strip().str.lower().isin({"true", "1", "1.0"})
    elif z_col is not None:
        valid_mask = df[z_col].notna()
        if not valid_mask.any():
            return None
        actual = df.loc[valid_mask, z_col].astype(float) > 0.0
    else:
        return None

    # Prediction: any predictor column true
    # After merge, predictor columns may have _align suffix
    available_predictors = []
    for pred_col in PREDICTOR_COLUMNS:
        # Try base name first, then _align suffix
        if pred_col in df.columns:
            available_predictors.append(pred_col)
        elif f"{pred_col}_align" in df.columns:
            available_predictors.append(f"{pred_col}_align")
    
    if not available_predictors:
        return None
    predicted = df.loc[valid_mask, available_predictors].any(axis=1)

    true_positives = int((actual & predicted).sum())
    false_positives = int((~actual & predicted).sum())
    false_negatives = int((actual & ~predicted).sum())

    denom = 2 * true_positives + false_positives + false_negatives
    f1 = (2 * true_positives / denom) if denom else 0.0

    return TargetResult(
        target="",
        atom_type=atom_type,
        f1=float(f1),
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        total_rows=int(valid_mask.sum()),
    )


def collect_1d_f1_results(
    outputs_dir: Path,
    allowed_targets: Optional[Dict[str, bool]],
) -> Tuple[List[TargetResult], List[TargetResult], List[TargetResult], List[TargetResult]]:
    """Collect per-target F1 scores for H, N, CA, and HA 1D CSPs."""
    results_H: List[TargetResult] = []
    results_N: List[TargetResult] = []
    results_CA: List[TargetResult] = []
    results_HA: List[TargetResult] = []

    for target_dir in discover_targets(outputs_dir):
        target_name = target_dir.name
        if allowed_targets is not None and target_name not in allowed_targets:
            continue

        merged = _merge_1d_with_alignment(target_dir)
        if merged is None:
            continue

        for atom, bucket in (("H", results_H), ("N", results_N), ("CA", results_CA), ("HA", results_HA)):
            res = _compute_f1_for_atom(merged, atom, target_name=target_name)
            if res is not None:
                res.target = target_name
                bucket.append(res)

    return results_H, results_N, results_CA, results_HA


def render_1d_f1_heatmaps(
    results_H: List[TargetResult],
    results_N: List[TargetResult],
    results_CA: List[TargetResult],
    results_HA: List[TargetResult],
    output_image: Path,
) -> None:
    """Render a 1x4 panel of F1 heatmaps for H, N, CA, and HA 1D CSPs."""
    atom_to_results = {
        "H": results_H,
        "N": results_N,
        "CA": results_CA,
        "HA": results_HA,
    }
    titles = {
        "H": "F1 Scores (H 1D CSPs)",
        "N": "F1 Scores (N 1D CSPs)",
        "CA": "F1 Scores (CA 1D CSPs)",
        "HA": "F1 Scores (HA 1D CSPs)",
    }

    fig, axes = plt.subplots(1, 4, figsize=(16, max(3, len(results_H) + len(results_N) + len(results_CA) + len(results_HA)) * 0.15))

    for idx, (atom, ax) in enumerate(zip(["H", "N", "CA", "HA"], axes)):
        res_list = atom_to_results[atom]
        if not res_list:
            ax.text(
                0.5,
                0.5,
                "No data",
                ha="center",
                va="center",
                fontsize=12,
            )
            ax.set_axis_off()
            continue

        df = pd.DataFrame([r.__dict__ for r in res_list])
        df = df.set_index("target").sort_values("f1", ascending=False)
        heatmap_data = df[["f1"]]

        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt=".2f",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            cbar=(idx == 3),
            ax=ax,
        )
        ax.set_xlabel("Metric", fontsize=10)
        ax.set_ylabel("Target", fontsize=10)
        ax.set_title(titles[atom], fontsize=12)
        ax.tick_params(axis="x", labelsize=8)
        ax.tick_params(axis="y", labelsize=8)

    plt.tight_layout()
    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def write_1d_f1_summary_csv(
    results_H: List[TargetResult],
    results_N: List[TargetResult],
    results_CA: List[TargetResult],
    results_HA: List[TargetResult],
    output_csv: Path,
) -> None:
    """Write summary statistics for 1D F1 score distributions per atom type."""

    def _summarize(atom_type: str, items: List[TargetResult]) -> Dict[str, object]:
        if not items:
            return {
                "atom_type": atom_type,
                "mean_f1": np.nan,
                "median_f1": np.nan,
                "std_f1": np.nan,
                "min_f1": np.nan,
                "max_f1": np.nan,
                "q1_f1": np.nan,
                "q3_f1": np.nan,
                "n_targets": 0,
            }
        vals = np.array([r.f1 for r in items], dtype=float)
        return {
            "atom_type": atom_type,
            "mean_f1": float(np.mean(vals)),
            "median_f1": float(np.median(vals)),
            "std_f1": float(np.std(vals, ddof=0)),
            "min_f1": float(np.min(vals)),
            "max_f1": float(np.max(vals)),
            "q1_f1": float(np.percentile(vals, 25)),
            "q3_f1": float(np.percentile(vals, 75)),
            "n_targets": int(len(vals)),
        }

    rows = [
        _summarize("H", results_H),
        _summarize("N", results_N),
        _summarize("CA", results_CA),
        _summarize("HA", results_HA),
    ]
    df = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)


def render_1d_f1_boxplot(
    results_H: List[TargetResult],
    results_N: List[TargetResult],
    results_CA: List[TargetResult],
    results_HA: List[TargetResult],
    output_image: Path,
) -> None:
    """Render a boxplot comparing F1 score distributions for H, N, CA, and HA 1D CSPs."""
    # Prepare data for boxplot
    data_for_plot = []
    labels = []
    
    if results_H:
        data_for_plot.append([r.f1 for r in results_H])
        labels.append("H")
    
    if results_N:
        data_for_plot.append([r.f1 for r in results_N])
        labels.append("N")
    
    if results_CA:
        data_for_plot.append([r.f1 for r in results_CA])
        labels.append("CA")

    if results_HA:
        data_for_plot.append([r.f1 for r in results_HA])
        labels.append("HA")
    
    if not data_for_plot:
        return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    bp = ax.boxplot(
        data_for_plot,
        labels=labels,
        patch_artist=True,
        showmeans=True,
        meanline=True,
    )
    
    # Customize boxplot colors
    colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"]
    for patch, color in zip(bp["boxes"], colors[:len(bp["boxes"])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_ylabel("F1 Score", fontsize=12)
    ax.set_xlabel("Atom Type", fontsize=12)
    ax.set_title("F1 Score Distribution for 1D Single-Atom CSPs", fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.set_ylim(0.0, 1.0)
    
    plt.tight_layout()
    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close()


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with open(path, "r", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _residue_key(row: Mapping[str, str]) -> Tuple[str, str]:
    return ((row.get("holo_resi") or "").strip(), (row.get("holo_aa") or "").strip())


def _index_by_residue(rows: List[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, str]]:
    lookup: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in rows:
        key = _residue_key(row)
        if key[0] and key[1]:
            lookup[key] = row
    return lookup


def _parse_shift(row: Optional[Mapping[str, str]], field: str) -> Optional[float]:
    if row is None:
        return None
    raw = (row.get(field) or "").strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _unreferenced_pair(
    primary: Mapping[str, str],
    fallback: Optional[Mapping[str, str]],
    apo_field: str,
    original_field: str,
) -> Optional[Tuple[float, float]]:
    apo = _parse_shift(primary, apo_field)
    original = _parse_shift(primary, original_field)
    if apo is None or original is None:
        apo = _parse_shift(fallback, apo_field) if apo is None else apo
        original = _parse_shift(fallback, original_field) if original is None else original
    if apo is None or original is None:
        return None
    return apo, original


def _save_1d_grid_csv(
    path: Path,
    *,
    atom: str,
    min_offset: float,
    max_offset: float,
    step: float,
    cutoff: float,
    grid_result: Mapping[str, object],
) -> None:
    offsets = list(grid_result["offset_values"])  # type: ignore[index]
    counts = list(grid_result["counts"])  # type: ignore[index]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        handle.write(f"# atom,{atom}\n")
        handle.write(f"# min,{min_offset}\n")
        handle.write(f"# max,{max_offset}\n")
        handle.write(f"# step,{step}\n")
        handle.write(f"# cutoff,{cutoff}\n")
        handle.write("offset_values," + ",".join(str(value) for value in offsets) + "\n")
        handle.write("counts," + ",".join(str(count) for count in counts) + "\n")
        handle.write(f"best_offset,{grid_result['best_offset']}\n")
        handle.write(f"best_count,{grid_result['best_count']}\n")


def compute_1d_metrics_for_target(target_dir: Path) -> None:
    """
    For a single target directory, compute 1D single-atom metrics and write
    ``1d_analysis.csv``.

    Each nucleus is referenced by its own 1D grid search on unreferenced
    apo/holo pairs (H/N from ``csp_table.csv``, CA from ``csp_table_CA.csv``,
    HA from ``csp_table_HA_CA.csv``), using the ``Referencing`` bounds for that
    nucleus. The search maximizes the count of ``|(holo + offset) - apo|``
    below ``grid_cutoff``. CSP is the absolute referenced difference. A residue
    is significant when that CSP is at least ``max(cleaned mean, 0.05 ppm)``.
    """
    csp_table_path = target_dir / "csp_table.csv"
    csp_table_ca_path = target_dir / "csp_table_CA.csv"
    csp_table_ha_ca_path = target_dir / "csp_table_HA_CA.csv"

    hn_rows = _read_csv_rows(csp_table_path)
    ca_rows = _read_csv_rows(csp_table_ca_path)
    ha_rows = _read_csv_rows(csp_table_ha_ca_path)
    if hn_rows:
        rows = [dict(row) for row in hn_rows]
    elif ca_rows:
        rows = [dict(row) for row in ca_rows]
    elif ha_rows:
        rows = [dict(row) for row in ha_rows]
    else:
        return

    ca_lookup = _index_by_residue(ca_rows)
    ha_lookup = _index_by_residue(ha_rows)
    ref = _Referencing()
    atom_sources = {
        "H": ("H_apo", "H_holo_original", ref.grid_h_min, ref.grid_h_max, ref.grid_h_step, ca_lookup),
        "N": ("N_apo", "N_holo_original", ref.grid_n_min, ref.grid_n_max, ref.grid_n_step, ca_lookup),
        "CA": ("CA_apo", "CA_holo_original", ref.grid_ca_min, ref.grid_ca_max, ref.grid_ca_step, ca_lookup),
        "HA": ("HA_apo", "HA_holo_original", ref.grid_ha_min, ref.grid_ha_max, ref.grid_ha_step, ha_lookup),
    }
    for atom, (apo_field, original_field, min_offset, max_offset, step, lookup) in atom_sources.items():
        pairs: List[Optional[Tuple[float, float]]] = []
        points: List[Tuple[float, float]] = []
        for row in rows:
            fallback = lookup.get(_residue_key(row))
            pair = _unreferenced_pair(row, fallback, apo_field, original_field)
            pairs.append(pair)
            if pair is not None:
                points.append(pair)
        row_offset = ""
        if points:
            grid_result = run_offset_grid_search_1d(
                points,
                min_offset=min_offset,
                max_offset=max_offset,
                step=step,
                cutoff=float(ref.grid_cutoff),
            )
            best_offset = float(grid_result["best_offset"])
            row_offset = f"{best_offset:.4f}"
            slug = (
                f"{atom}_{_fmt_grid_num(min_offset)}_{_fmt_grid_num(max_offset)}_"
                f"{_fmt_grid_num(step)}__C_{_fmt_grid_num(ref.grid_cutoff)}"
            )
            _save_1d_grid_csv(
                target_dir / f"offset_grid_1d_{slug}.csv",
                atom=atom,
                min_offset=min_offset,
                max_offset=max_offset,
                step=step,
                cutoff=float(ref.grid_cutoff),
                grid_result=grid_result,
            )
        else:
            best_offset = 0.0
        for row, pair in zip(rows, pairs):
            if pair is None:
                row[apo_field] = row.get(apo_field, "")
                row[original_field] = ""
                row[f"{atom}_offset"] = ""
                row[f"{atom}_holo"] = ""
                continue
            apo_shift, original_shift = pair
            row[apo_field] = f"{apo_shift:.4f}"
            row[original_field] = f"{original_shift:.4f}"
            row[f"{atom}_offset"] = row_offset
            row[f"{atom}_holo"] = f"{original_shift + best_offset:.4f}"

    if not rows:
        return

    # Collect per-atom CSP_1d values to compute per-target mean/sd
    h_csp_vals: List[float] = []
    n_csp_vals: List[float] = []
    ca_csp_vals: List[float] = []
    ha_csp_vals: List[float] = []

    per_row_stats: List[Dict[str, Atom1DStats]] = []

    for row in rows:
        stats_for_row: Dict[str, Atom1DStats] = {}

        # Helper to parse float safely
        def _get_float(field: str) -> Optional[float]:
            val = row.get(field, "").strip()
            if not val:
                return None
            try:
                return float(val)
            except ValueError:
                return None

        # H
        H_apo = _get_float("H_apo")
        H_holo = _get_float("H_holo")
        if H_apo is not None and H_holo is not None:
            dH = H_holo - H_apo
            csp_H_1d = abs(dH)
            stats_for_row["H"] = Atom1DStats(delta=dH, csp_1d=csp_H_1d, z_1d=None)
            h_csp_vals.append(csp_H_1d)
        else:
            stats_for_row["H"] = Atom1DStats(delta=None, csp_1d=None, z_1d=None)

        # N
        N_apo = _get_float("N_apo")
        N_holo = _get_float("N_holo")
        if N_apo is not None and N_holo is not None:
            dN = N_holo - N_apo
            csp_N_1d = abs(dN)
            stats_for_row["N"] = Atom1DStats(delta=dN, csp_1d=csp_N_1d, z_1d=None)
            n_csp_vals.append(csp_N_1d)
        else:
            stats_for_row["N"] = Atom1DStats(delta=None, csp_1d=None, z_1d=None)

        CA_apo = _get_float("CA_apo")
        CA_holo = _get_float("CA_holo")
        if CA_apo is not None and CA_holo is not None:
            dCA = CA_holo - CA_apo
            csp_CA_1d = abs(dCA)
            stats_for_row["CA"] = Atom1DStats(delta=dCA, csp_1d=csp_CA_1d, z_1d=None)
            ca_csp_vals.append(csp_CA_1d)
        else:
            stats_for_row["CA"] = Atom1DStats(delta=None, csp_1d=None, z_1d=None)

        HA_apo = _get_float("HA_apo")
        HA_holo = _get_float("HA_holo")
        if HA_apo is not None and HA_holo is not None:
            dHA = HA_holo - HA_apo
            csp_HA_1d = abs(dHA)
            stats_for_row["HA"] = Atom1DStats(delta=dHA, csp_1d=csp_HA_1d, z_1d=None)
            ha_csp_vals.append(csp_HA_1d)
        else:
            stats_for_row["HA"] = Atom1DStats(delta=None, csp_1d=None, z_1d=None)

        per_row_stats.append(stats_for_row)

    # Z-scores use the full-sample mean and SD. Significance uses the
    # outlier-cleaned mean, floored at 0.05 ppm.
    h_mean, h_sd = compute_atom_stats(h_csp_vals)
    n_mean, n_sd = compute_atom_stats(n_csp_vals)
    ca_mean, ca_sd = compute_atom_stats(ca_csp_vals)
    ha_mean, ha_sd = compute_atom_stats(ha_csp_vals)

    h_threshold_info = compute_threshold_with_outlier_removal(
        h_csp_vals,
        _thresholds.outlier_z_score,
        _thresholds.significance_z_score,
        _thresholds.max_outlier_iterations,
        _thresholds.max_outlier_fraction,
    )
    n_threshold_info = compute_threshold_with_outlier_removal(
        n_csp_vals,
        _thresholds.outlier_z_score,
        _thresholds.significance_z_score,
        _thresholds.max_outlier_iterations,
        _thresholds.max_outlier_fraction,
    )
    ca_threshold_info = compute_threshold_with_outlier_removal(
        ca_csp_vals,
        _thresholds.outlier_z_score,
        _thresholds.significance_z_score,
        _thresholds.max_outlier_iterations,
        _thresholds.max_outlier_fraction,
    )
    ha_threshold_info = compute_threshold_with_outlier_removal(
        ha_csp_vals,
        _thresholds.outlier_z_score,
        _thresholds.significance_z_score,
        _thresholds.max_outlier_iterations,
        _thresholds.max_outlier_fraction,
    )

    h_cutoff = _floor_primary_hn_cutoff(float(h_threshold_info.threshold)) if h_csp_vals else None
    n_cutoff = _floor_primary_hn_cutoff(float(n_threshold_info.threshold)) if n_csp_vals else None
    ca_cutoff = _floor_primary_hn_cutoff(float(ca_threshold_info.threshold)) if ca_csp_vals else None
    ha_cutoff = _floor_primary_hn_cutoff(float(ha_threshold_info.threshold)) if ha_csp_vals else None

    for stats_for_row in per_row_stats:
        # H
        if stats_for_row["H"].csp_1d is not None and h_sd > 0.0:
            stats_for_row["H"].z_1d = (stats_for_row["H"].csp_1d - h_mean) / h_sd
        elif stats_for_row["H"].csp_1d is not None:
            stats_for_row["H"].z_1d = 0.0
        if stats_for_row["H"].csp_1d is not None and h_cutoff is not None:
            stats_for_row["H"].significant = bool(stats_for_row["H"].csp_1d >= h_cutoff)

        # N
        if stats_for_row["N"].csp_1d is not None and n_sd > 0.0:
            stats_for_row["N"].z_1d = (stats_for_row["N"].csp_1d - n_mean) / n_sd
        elif stats_for_row["N"].csp_1d is not None:
            stats_for_row["N"].z_1d = 0.0
        if stats_for_row["N"].csp_1d is not None and n_cutoff is not None:
            stats_for_row["N"].significant = bool(stats_for_row["N"].csp_1d >= n_cutoff)

        # CA
        if stats_for_row["CA"].csp_1d is not None and ca_sd > 0.0:
            stats_for_row["CA"].z_1d = (stats_for_row["CA"].csp_1d - ca_mean) / ca_sd
        elif stats_for_row["CA"].csp_1d is not None:
            stats_for_row["CA"].z_1d = 0.0
        if stats_for_row["CA"].csp_1d is not None and ca_cutoff is not None:
            stats_for_row["CA"].significant = bool(stats_for_row["CA"].csp_1d >= ca_cutoff)

        # HA
        if stats_for_row["HA"].csp_1d is not None and ha_sd > 0.0:
            stats_for_row["HA"].z_1d = (stats_for_row["HA"].csp_1d - ha_mean) / ha_sd
        elif stats_for_row["HA"].csp_1d is not None:
            stats_for_row["HA"].z_1d = 0.0
        if stats_for_row["HA"].csp_1d is not None and ha_cutoff is not None:
            stats_for_row["HA"].significant = bool(stats_for_row["HA"].csp_1d >= ha_cutoff)

    # Write per-target 1d_analysis.csv
    output_path = target_dir / "1d_analysis.csv"
    fieldnames = [
        "apo_bmrb",
        "holo_bmrb",
        "holo_pdb",
        "chain",
        "apo_resi",
        "apo_aa",
        "holo_resi",
        "holo_aa",
        # Generic sequence columns used by merge_csv alignment
        "resi",
        "aa",
        # H metrics
        "H_apo",
        "H_holo_original",
        "H_offset",
        "H_holo",
        "dH_1d",
        "CSP_H_1d",
        "z_H_1d",
        "csp_H_1d_significant",
        # N metrics
        "N_apo",
        "N_holo_original",
        "N_offset",
        "N_holo",
        "dN_1d",
        "CSP_N_1d",
        "z_N_1d",
        "csp_N_1d_significant",
        # CA metrics
        "CA_apo",
        "CA_holo_original",
        "CA_offset",
        "CA_holo",
        "dCA_1d",
        "CSP_CA_1d",
        "z_CA_1d",
        "csp_CA_1d_significant",
        # HA metrics
        "HA_apo",
        "HA_holo_original",
        "HA_offset",
        "HA_holo",
        "dHA_1d",
        "CSP_HA_1d",
        "z_HA_1d",
        "csp_HA_1d_significant",
    ]

    with open(output_path, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        for row, stats_for_row in zip(rows, per_row_stats):
            out_row: Dict[str, object] = {}
            # Identification
            out_row["apo_bmrb"] = row.get("apo_bmrb", "")
            out_row["holo_bmrb"] = row.get("holo_bmrb", "")
            out_row["holo_pdb"] = row.get("holo_pdb", "")
            out_row["chain"] = row.get("chain", "")
            out_row["apo_resi"] = row.get("apo_resi", "")
            out_row["apo_aa"] = row.get("apo_aa", "")
            out_row["holo_resi"] = row.get("holo_resi", "")
            out_row["holo_aa"] = row.get("holo_aa", "")
            # Alignment helper columns (use holo sequence)
            out_row["resi"] = row.get("holo_resi", "")
            out_row["aa"] = row.get("holo_aa", "")

            # H
            out_row["H_apo"] = row.get("H_apo", "")
            out_row["H_holo_original"] = row.get("H_holo_original", "")
            out_row["H_offset"] = row.get("H_offset", "")
            out_row["H_holo"] = row.get("H_holo", "")
            h_stats = stats_for_row["H"]
            out_row["dH_1d"] = f"{h_stats.delta:.4f}" if h_stats.delta is not None else ""
            out_row["CSP_H_1d"] = f"{h_stats.csp_1d:.4f}" if h_stats.csp_1d is not None else ""
            out_row["z_H_1d"] = f"{h_stats.z_1d:.4f}" if h_stats.z_1d is not None else ""
            out_row["csp_H_1d_significant"] = (
                "True" if h_stats.significant is True else ("False" if h_stats.significant is False else "")
            )

            # N
            out_row["N_apo"] = row.get("N_apo", "")
            out_row["N_holo_original"] = row.get("N_holo_original", "")
            out_row["N_offset"] = row.get("N_offset", "")
            out_row["N_holo"] = row.get("N_holo", "")
            n_stats = stats_for_row["N"]
            out_row["dN_1d"] = f"{n_stats.delta:.4f}" if n_stats.delta is not None else ""
            out_row["CSP_N_1d"] = f"{n_stats.csp_1d:.4f}" if n_stats.csp_1d is not None else ""
            out_row["z_N_1d"] = f"{n_stats.z_1d:.4f}" if n_stats.z_1d is not None else ""
            out_row["csp_N_1d_significant"] = (
                "True" if n_stats.significant is True else ("False" if n_stats.significant is False else "")
            )

            # CA (may be absent)
            out_row["CA_apo"] = row.get("CA_apo", "")
            out_row["CA_holo_original"] = row.get("CA_holo_original", "")
            out_row["CA_offset"] = row.get("CA_offset", "")
            out_row["CA_holo"] = row.get("CA_holo", "")
            ca_stats = stats_for_row["CA"]
            out_row["dCA_1d"] = f"{ca_stats.delta:.4f}" if ca_stats.delta is not None else ""
            out_row["CSP_CA_1d"] = f"{ca_stats.csp_1d:.4f}" if ca_stats.csp_1d is not None else ""
            out_row["z_CA_1d"] = f"{ca_stats.z_1d:.4f}" if ca_stats.z_1d is not None else ""
            out_row["csp_CA_1d_significant"] = (
                "True" if ca_stats.significant is True else ("False" if ca_stats.significant is False else "")
            )

            # HA
            out_row["HA_apo"] = row.get("HA_apo", "")
            out_row["HA_holo_original"] = row.get("HA_holo_original", "")
            out_row["HA_offset"] = row.get("HA_offset", "")
            out_row["HA_holo"] = row.get("HA_holo", "")
            ha_stats = stats_for_row["HA"]
            out_row["dHA_1d"] = f"{ha_stats.delta:.4f}" if ha_stats.delta is not None else ""
            out_row["CSP_HA_1d"] = f"{ha_stats.csp_1d:.4f}" if ha_stats.csp_1d is not None else ""
            out_row["z_HA_1d"] = f"{ha_stats.z_1d:.4f}" if ha_stats.z_1d is not None else ""
            out_row["csp_HA_1d_significant"] = (
                "True" if ha_stats.significant is True else ("False" if ha_stats.significant is False else "")
            )

            writer.writerow(out_row)


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir.resolve()
    output_image = args.output_image.resolve()
    summary_csv = args.summary_csv.resolve()
    allowed_targets = load_allowed_targets(args.targets_csv.resolve()) if args.targets_csv else None

    # Determine boxplot output path
    if args.boxplot_image:
        boxplot_image = args.boxplot_image.resolve()
    else:
        # Default to same directory as summary_csv with descriptive name
        boxplot_image = summary_csv.parent / "f1_boxplot_1d_single_atom.png"

    for target_dir in discover_targets(outputs_dir):
        target_name = target_dir.name
        if allowed_targets is not None and target_name not in allowed_targets:
            continue
        try:
            compute_1d_metrics_for_target(target_dir)
        except Exception as exc:
            print(f"[WARN] Failed to compute 1D metrics for {target_name}: {exc}")

    # Collect F1 results across targets
    results_H, results_N, results_CA, results_HA = collect_1d_f1_results(outputs_dir, allowed_targets)
    if not results_H and not results_N and not results_CA and not results_HA:
        return 0

    render_1d_f1_heatmaps(results_H, results_N, results_CA, results_HA, output_image)
    write_1d_f1_summary_csv(results_H, results_N, results_CA, results_HA, summary_csv)
    render_1d_f1_boxplot(results_H, results_N, results_CA, results_HA, boxplot_image)

    return 0


if __name__ == "__main__":
    import sys as _sys

    raise SystemExit(main(_sys.argv[1:]))
