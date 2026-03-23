#!/usr/bin/env python3
"""
Per-target single-atom (H, N, CA) 1D shift perturbation analysis.

For each target directory under outputs/ that has a csp_table.csv, this script:
  - Loads csp_table.csv (and csp_table_CA.csv if present)
  - Computes per-residue 1D shift metrics for each available atom type:
      ΔH = H_holo - H_apo, ΔN = N_holo - N_apo, ΔCA = CA_holo - CA_apo
      CSP_H_1d = |ΔH|, CSP_N_1d = |ΔN|, CSP_CA_1d = |ΔCA|
      z-scores per atom type, using target-specific mean and SD of CSP_X_1d
  - Writes outputs/{holo_pdb}/1d_analysis.csv with one row per residue
    containing all single-atom 1D metrics for that residue.

No scaling coefficients are applied to the 1D CSPs; they are plain absolute
differences in ppm for each atom type.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import statistics
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


@dataclass
class Atom1DStats:
    delta: Optional[float]
    csp_1d: Optional[float]
    z_1d: Optional[float]


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-residue 1D single-atom shift perturbations (H, N, CA) "
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
    Compute F1 score for a single atom type using z_X_1d > 0 as ground truth
    and predictor columns as the prediction.
    """
    # After merge, z-score columns may have suffixes. For CA in particular,
    # the reliable z-scores often come from master_alignment.csv rather than 1d_analysis.csv.
    base_z_col = f"z_{atom_type}_1d"

    # Prefer columns that actually contain non-null data. Try several naming patterns:
    #   - '{base}_align'  (from master_alignment.csv after merge)
    #   - '{base}_1d'     (from 1d_analysis.csv after merge)
    #   - '{base}'        (either side if no suffix was applied)
    z_col: Optional[str] = None
    for candidate in [f"{base_z_col}_align", f"{base_z_col}_1d", base_z_col]:
        exists = candidate in df.columns
        has_data = exists and df[candidate].notna().any() if exists else False
        if exists and has_data:
            z_col = candidate
            break

    if z_col is None:
        return None

    # Ground truth: z-score > 0
    valid_mask = df[z_col].notna()
    if not valid_mask.any():
        return None

    z_vals = df.loc[valid_mask, z_col].astype(float)
    actual = z_vals > 0.0

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
) -> Tuple[List[TargetResult], List[TargetResult], List[TargetResult]]:
    """Collect per-target F1 scores for H, N, and CA 1D CSPs."""
    results_H: List[TargetResult] = []
    results_N: List[TargetResult] = []
    results_CA: List[TargetResult] = []

    for target_dir in discover_targets(outputs_dir):
        target_name = target_dir.name
        if allowed_targets is not None and target_name not in allowed_targets:
            continue

        merged = _merge_1d_with_alignment(target_dir)
        if merged is None:
            continue

        for atom, bucket in (("H", results_H), ("N", results_N), ("CA", results_CA)):
            res = _compute_f1_for_atom(merged, atom, target_name=target_name)
            if res is not None:
                res.target = target_name
                bucket.append(res)

    return results_H, results_N, results_CA


def render_1d_f1_heatmaps(
    results_H: List[TargetResult],
    results_N: List[TargetResult],
    results_CA: List[TargetResult],
    output_image: Path,
) -> None:
    """Render a 1x3 panel of F1 heatmaps for H, N, and CA 1D CSPs."""
    atom_to_results = {
        "H": results_H,
        "N": results_N,
        "CA": results_CA,
    }
    titles = {
        "H": "F1 Scores (H 1D CSPs)",
        "N": "F1 Scores (N 1D CSPs)",
        "CA": "F1 Scores (CA 1D CSPs)",
    }

    fig, axes = plt.subplots(1, 3, figsize=(12, max(3, len(results_H) + len(results_N) + len(results_CA)) * 0.15))

    for idx, (atom, ax) in enumerate(zip(["H", "N", "CA"], axes)):
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
            cbar=(idx == 2),
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
    ]
    df = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)


def render_1d_f1_boxplot(
    results_H: List[TargetResult],
    results_N: List[TargetResult],
    results_CA: List[TargetResult],
    output_image: Path,
) -> None:
    """Render a boxplot comparing F1 score distributions for H, N, and CA 1D CSPs."""
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
    colors = ["#66c2a5", "#fc8d62", "#8da0cb"]
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


def compute_1d_metrics_for_target(target_dir: Path) -> None:
    """
    For a single target directory, compute 1D single-atom metrics and write 1d_analysis.csv.
    """
    csp_table_path = target_dir / "csp_table.csv"
    if not csp_table_path.exists():
        return

    # Optional CA-inclusive table; used to extend CA fields if needed.
    csp_table_ca_path = target_dir / "csp_table_CA.csv"

    rows: List[Dict[str, str]] = []
    with open(csp_table_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(dict(row))

    # Build CA lookup from csp_table_CA.csv keyed by (holo_resi, holo_aa)
    ca_lookup: Dict[Tuple[str, str], Dict[str, str]] = {}
    if csp_table_ca_path.exists():
        with open(csp_table_ca_path, "r", newline="") as f_ca:
            ca_reader = csv.DictReader(f_ca)
            for ca_row in ca_reader:
                key = (
                    (ca_row.get("holo_resi") or "").strip(),
                    (ca_row.get("holo_aa") or "").strip(),
                )
                if key[0] and key[1]:
                    ca_lookup[key] = dict(ca_row)

    if not rows:
        return

    # Collect per-atom CSP_1d values to compute per-target mean/sd
    h_csp_vals: List[float] = []
    n_csp_vals: List[float] = []
    ca_csp_vals: List[float] = []

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

        # CA: may come from csp_table_CA.csv if not in csp_table.csv
        CA_apo = _get_float("CA_apo")
        CA_holo = _get_float("CA_holo")

        if CA_apo is None or CA_holo is None:
            # Try to pull CA values from csp_table_CA.csv
            key = (
                (row.get("holo_resi") or "").strip(),
                (row.get("holo_aa") or "").strip(),
            )
            ca_row = ca_lookup.get(key)
            if ca_row is not None:
                ca_apo_raw = (ca_row.get("CA_apo") or "").strip()
                ca_holo_raw = (ca_row.get("CA_holo") or "").strip()
                try:
                    if CA_apo is None and ca_apo_raw:
                        CA_apo = float(ca_apo_raw)
                        row["CA_apo"] = ca_apo_raw
                except ValueError:
                    CA_apo = None
                try:
                    if CA_holo is None and ca_holo_raw:
                        CA_holo = float(ca_holo_raw)
                        row["CA_holo"] = ca_holo_raw
                except ValueError:
                    CA_holo = None

                # Also propagate CA_offset if present
                ca_offset_raw = (ca_row.get("CA_offset") or "").strip()
                if ca_offset_raw and not row.get("CA_offset"):
                    row["CA_offset"] = ca_offset_raw
        if CA_apo is not None and CA_holo is not None:
            dCA = CA_holo - CA_apo
            csp_CA_1d = abs(dCA)
            stats_for_row["CA"] = Atom1DStats(delta=dCA, csp_1d=csp_CA_1d, z_1d=None)
            ca_csp_vals.append(csp_CA_1d)
        else:
            stats_for_row["CA"] = Atom1DStats(delta=None, csp_1d=None, z_1d=None)

        per_row_stats.append(stats_for_row)

    # Compute per-target mean/sd and z-scores
    h_mean, h_sd = compute_atom_stats(h_csp_vals)
    n_mean, n_sd = compute_atom_stats(n_csp_vals)
    ca_mean, ca_sd = compute_atom_stats(ca_csp_vals)

    idx_h = 0
    idx_n = 0
    idx_ca = 0

    for i, stats_for_row in enumerate(per_row_stats):
        # H
        if stats_for_row["H"].csp_1d is not None and h_sd > 0.0:
            z = (stats_for_row["H"].csp_1d - h_mean) / h_sd
            stats_for_row["H"].z_1d = z
        elif stats_for_row["H"].csp_1d is not None:
            stats_for_row["H"].z_1d = 0.0

        # N
        if stats_for_row["N"].csp_1d is not None and n_sd > 0.0:
            z = (stats_for_row["N"].csp_1d - n_mean) / n_sd
            stats_for_row["N"].z_1d = z
        elif stats_for_row["N"].csp_1d is not None:
            stats_for_row["N"].z_1d = 0.0

        # CA
        if stats_for_row["CA"].csp_1d is not None and ca_sd > 0.0:
            z = (stats_for_row["CA"].csp_1d - ca_mean) / ca_sd
            stats_for_row["CA"].z_1d = z
        elif stats_for_row["CA"].csp_1d is not None:
            stats_for_row["CA"].z_1d = 0.0

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
        "H_holo",
        "H_offset",
        "dH_1d",
        "CSP_H_1d",
        "z_H_1d",
        # N metrics
        "N_apo",
        "N_holo",
        "N_offset",
        "dN_1d",
        "CSP_N_1d",
        "z_N_1d",
        # CA metrics
        "CA_apo",
        "CA_holo",
        "CA_offset",
        "dCA_1d",
        "CSP_CA_1d",
        "z_CA_1d",
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
            out_row["H_holo"] = row.get("H_holo", "")
            out_row["H_offset"] = row.get("H_offset", "")
            h_stats = stats_for_row["H"]
            out_row["dH_1d"] = f"{h_stats.delta:.4f}" if h_stats.delta is not None else ""
            out_row["CSP_H_1d"] = f"{h_stats.csp_1d:.4f}" if h_stats.csp_1d is not None else ""
            out_row["z_H_1d"] = f"{h_stats.z_1d:.4f}" if h_stats.z_1d is not None else ""

            # N
            out_row["N_apo"] = row.get("N_apo", "")
            out_row["N_holo"] = row.get("N_holo", "")
            out_row["N_offset"] = row.get("N_offset", "")
            n_stats = stats_for_row["N"]
            out_row["dN_1d"] = f"{n_stats.delta:.4f}" if n_stats.delta is not None else ""
            out_row["CSP_N_1d"] = f"{n_stats.csp_1d:.4f}" if n_stats.csp_1d is not None else ""
            out_row["z_N_1d"] = f"{n_stats.z_1d:.4f}" if n_stats.z_1d is not None else ""

            # CA (may be absent)
            out_row["CA_apo"] = row.get("CA_apo", "")
            out_row["CA_holo"] = row.get("CA_holo", "")
            out_row["CA_offset"] = row.get("CA_offset", "")
            ca_stats = stats_for_row["CA"]
            out_row["dCA_1d"] = f"{ca_stats.delta:.4f}" if ca_stats.delta is not None else ""
            out_row["CSP_CA_1d"] = f"{ca_stats.csp_1d:.4f}" if ca_stats.csp_1d is not None else ""
            out_row["z_CA_1d"] = f"{ca_stats.z_1d:.4f}" if ca_stats.z_1d is not None else ""

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
    results_H, results_N, results_CA = collect_1d_f1_results(outputs_dir, allowed_targets)
    if not results_H and not results_N and not results_CA:
        return 0

    render_1d_f1_heatmaps(results_H, results_N, results_CA, output_image)
    write_1d_f1_summary_csv(results_H, results_N, results_CA, summary_csv)
    render_1d_f1_boxplot(results_H, results_N, results_CA, boxplot_image)

    return 0


if __name__ == "__main__":
    import sys as _sys

    raise SystemExit(main(_sys.argv[1:]))


