#!/usr/bin/env python3
"""
Analyze CA-inclusive CSP outputs and visualize per-target F1 scores.

This script is similar to analyze_targets.py but:
1. Uses csp_table_CA.csv as the base (instead of master_alignment.csv)
2. Uses csp_CA_significant instead of significant
3. Only processes targets listed in targets_with_ca_shifts.csv
4. Merges CA CSP data with other analysis files (occlusion, interaction, distance)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from math import ceil
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Support running as a script or module
try:
    from .align import align_global
    from .config import classification_colors
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.align import align_global
    from scripts.config import classification_colors


SIGNIFICANT_COLUMN = "csp_CA_significant"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


class AlignmentParsingError(RuntimeError):
    """Raised when a CSV file cannot be parsed."""


@dataclass
class TargetResult:
    target: str
    f1: float
    true_positives: int
    false_positives: int
    false_negatives: int
    total_rows: int


@dataclass
class DistanceRecord:
    distance: float
    is_predicted_positive: bool


@dataclass
class ConfusionRecord:
    """Record for confusion matrix: ground truth = binding site, prediction = CA CSP significance."""
    distance: float
    is_binding: bool
    is_significant: bool


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-target F1 scores from CA-inclusive CSP data (csp_table_CA.csv) "
            "and render them as a heatmap. Only processes targets listed in targets_with_ca_shifts.csv."
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
        required=True,
        help="CSV file containing target IDs to process. Must have a 'holo_pdb' column.",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("outputs") / "f1_heatmap_ca.png",
        help="Destination for the generated heatmap image (default: %(default)s).",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write the per-target F1 summary table as CSV.",
    )
    parser.add_argument(
        "--histogram-image",
        type=Path,
        default=Path("outputs") / "significant_ca_distance_hist_ca.png",
        help=(
            "Destination for the histogram of significant residues by CA distance "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--stacked-histogram-image",
        type=Path,
        default=Path("outputs") / "significant_ca_distance_stacked_hist_ca.png",
        help=(
            "Destination for the stacked histogram splitting predicted positives "
            "and negatives among significant residues (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--confusion-matrix-stacked-histogram-image",
        type=Path,
        default=Path("outputs") / "confusion_matrix_stacked_histogram_ca.png",
        help="Destination for confusion matrix stacked histogram (CA-inclusive CSPs).",
    )
    parser.add_argument(
        "--scatterplot-image",
        type=Path,
        default=Path("outputs") / "f1_comparison_scatterplot.png",
        help=(
            "Destination for the scatterplot comparing F1 scores from CA-inclusive "
            "vs N/H CSPs (default: %(default)s)."
        ),
    )
    return parser.parse_args(list(argv))


def load_targets_from_csv(csv_path: Path) -> Set[str]:
    """
    Load allowed target IDs from a CSV file.
    
    Args:
        csv_path: Path to CSV file with 'holo_pdb' column
        
    Returns:
        Set of target IDs (holo_pdb values)
    """
    try:
        df = pd.read_csv(csv_path)
        if 'holo_pdb' not in df.columns:
            raise ValueError(f"CSV file {csv_path} must have a 'holo_pdb' column")
        # Convert to set and strip whitespace
        targets = set(df['holo_pdb'].astype(str).str.strip())
        return targets
    except Exception as exc:
        raise ValueError(f"Failed to load targets from {csv_path}: {exc}") from exc


def load_ca_csp_reference(csp_table_path: Path) -> Tuple[str, List[int], Dict[int, Dict[str, Any]]]:
    """
    Load CA CSP table and extract the reference sequence.
    
    Args:
        csp_table_path: Path to csp_table_CA.csv
        
    Returns:
        Tuple of (reference_sequence, sequential_positions, csp_data)
    """
    if not csp_table_path.exists():
        raise FileNotFoundError(f"CA CSP table not found: {csp_table_path}")
    
    reference_sequence = []
    sequential_positions = []
    csp_data = {}
    
    # Amino acid mapping from 3-letter to 1-letter codes
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    with open(csp_table_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                holo_resi = int(row['holo_resi'])
                holo_aa = row['holo_aa'].strip()
                
                # Only include rows with valid amino acid data
                if holo_aa and holo_aa != 'P' and holo_aa != '':  # Skip prolines and empty
                    # Convert 3-letter to 1-letter if needed
                    if len(holo_aa) == 3:
                        aa_letter = aa_3to1.get(holo_aa, 'X')
                    else:
                        aa_letter = holo_aa
                    
                    reference_sequence.append(aa_letter)
                    sequential_positions.append(holo_resi)
                    csp_data[holo_resi] = dict(row)
            except (ValueError, KeyError) as e:
                continue
    
    ref_seq_str = ''.join(reference_sequence)
    return ref_seq_str, sequential_positions, csp_data


def load_and_align_csv(csv_path: Path, ref_sequence: str, ref_positions: List[int], 
                      csv_type: str = "unknown") -> Dict[int, Dict[str, Any]]:
    """
    Load other CSV files and align to reference sequence.
    """
    if not csv_path.exists():
        return {}
    
    # Extract sequence from CSV file
    csv_sequence = []
    csv_positions = []
    csv_data = {}
    
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # Try different column name patterns
                residue_number = None
                residue_name = None
                
                if 'residue_number' in row and 'residue_name' in row:
                    residue_number = int(row['residue_number'])
                    residue_name = row['residue_name'].strip()
                elif 'resi' in row and 'aa' in row:
                    residue_number = int(row['resi'])
                    residue_name = row['aa'].strip()
                else:
                    continue
                
                if residue_name and residue_name != 'PRO':  # Skip prolines
                    # Convert 3-letter to 1-letter if needed
                    if len(residue_name) == 3:
                        aa_letter = aa_3to1.get(residue_name, 'X')
                    else:
                        aa_letter = residue_name
                    
                    csv_sequence.append(aa_letter)
                    csv_positions.append(residue_number)
                    csv_data[residue_number] = dict(row)
                    
            except (ValueError, KeyError):
                continue
    
    csv_seq_str = ''.join(csv_sequence)
    
    if not csv_seq_str:
        return {}
    
    # Perform alignment
    aligned_ref, aligned_csv, mapping, score = align_global(ref_sequence, csv_seq_str)
    
    # Map aligned positions back to sequential positions
    aligned_data = {}
    for ref_pos, csv_pos in mapping:
        if ref_pos <= len(ref_positions) and csv_pos <= len(csv_positions):
            sequential_position = ref_positions[ref_pos - 1]  # Convert to 1-based
            csv_position = csv_positions[csv_pos - 1]  # Convert to 1-based
            
            if csv_position in csv_data:
                aligned_data[sequential_position] = csv_data[csv_position]
    
    return aligned_data


def to_bool(value) -> bool:
    """Coerce assorted representations into booleans."""
    if pd.isna(value):
        return False
    if isinstance(value, (bool, int)):
        return bool(value)
    if isinstance(value, float):
        return bool(round(value))
    if isinstance(value, str):
        text = value.strip().lower()
        if not text:
            return False
        if text in {"true", "t", "yes", "y", "1"}:
            return True
        if text in {"false", "f", "no", "n", "0"}:
            return False
    raise AlignmentParsingError(f"Unable to interpret value {value!r} as boolean.")


def load_ca_alignment(target_dir: Path) -> pd.DataFrame:
    """
    Load and merge CA CSP data with other analysis files.
    
    Args:
        target_dir: Directory containing CSV files for a target
        
    Returns:
        DataFrame with merged CA CSP and analysis data
    """
    # Load CA CSP table as reference
    csp_table_path = target_dir / "csp_table_CA.csv"
    if not csp_table_path.exists():
        raise AlignmentParsingError(f"CA CSP table not found: {csp_table_path}")
    
    ref_sequence, ref_positions, csp_data = load_ca_csp_reference(csp_table_path)
    
    if not ref_sequence:
        raise AlignmentParsingError("No valid CA CSP reference sequence found")
    
    # Load and align other CSV files
    csv_types = {
        'occlusion_analysis.csv': 'occlusion',
        'interaction_filter.csv': 'interaction', 
        'ca_distance_filter.csv': 'distance'
    }
    
    aligned_data = {}
    for csv_file, csv_type in csv_types.items():
        csv_path = target_dir / csv_file
        if csv_path.exists():
            aligned_csv_data = load_and_align_csv(csv_path, ref_sequence, ref_positions, csv_type)
            aligned_data[csv_type] = aligned_csv_data
    
    # Merge all data into a DataFrame
    rows = []
    for seq_pos in ref_positions:
        row = {}
        
        # Add CA CSP data
        if seq_pos in csp_data:
            row.update(csp_data[seq_pos])
        
        # Add aligned data from other CSV files
        for csv_type, data in aligned_data.items():
            if seq_pos in data:
                csv_row = data[seq_pos]
                for col_name, value in csv_row.items():
                    if col_name not in ['residue_number', 'residue_name', 'resi', 'aa']:
                        prefixed_name = f"{col_name}_{csv_type}"
                        row[prefixed_name] = value
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Convert boolean columns
    if SIGNIFICANT_COLUMN in df.columns:
        df[SIGNIFICANT_COLUMN] = df[SIGNIFICANT_COLUMN].apply(to_bool)
    
    for col in PREDICTOR_COLUMNS:
        if col in df.columns:
            df[col] = df[col].apply(to_bool)
    
    # Convert CA distance column
    if CA_DISTANCE_COLUMN in df.columns:
        df[CA_DISTANCE_COLUMN] = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")
    
    return df


def compute_f1_score(df: pd.DataFrame, predicted: pd.Series | None = None) -> TargetResult:
    actual = df[SIGNIFICANT_COLUMN]
    if predicted is None:
        # Only use predictor columns that exist in the dataframe
        available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
        if not available_predictors:
            predicted = pd.Series([False] * len(df), index=df.index)
        else:
            predicted = df[available_predictors].any(axis=1)

    true_positives = int((actual & predicted).sum())
    false_positives = int((~actual & predicted).sum())
    false_negatives = int((actual & ~predicted).sum())

    denominator = 2 * true_positives + false_positives + false_negatives
    f1 = (2 * true_positives / denominator) if denominator else 0.0

    return TargetResult(
        target="",
        f1=f1,
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        total_rows=int(len(df)),
    )


def discover_targets(outputs_dir: Path) -> List[Path]:
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")

    return sorted(
        path for path in outputs_dir.iterdir() if path.is_dir() and not path.name.startswith(".")
    )


def load_nh_alignment(alignment_path: Path) -> pd.DataFrame:
    """
    Load N/H CSP alignment data from master_alignment.csv.
    
    Args:
        alignment_path: Path to master_alignment.csv
        
    Returns:
        DataFrame with N/H CSP and analysis data
    """
    try:
        df = pd.read_csv(alignment_path)
    except Exception as exc:
        raise AlignmentParsingError(
            f"Failed to read alignment CSV: {alignment_path}"
        ) from exc

    # Use 'significant' column for N/H CSPs
    nh_significant_col = "significant"
    required_columns = (nh_significant_col, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS)
    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]
    if missing_columns:
        raise AlignmentParsingError(
            f"Alignment file {alignment_path} is missing required columns: "
            f"{', '.join(missing_columns)}"
        )

    for column in (nh_significant_col, *PREDICTOR_COLUMNS):
        if column in df.columns:
            df[column] = df[column].apply(to_bool)

    if CA_DISTANCE_COLUMN in df.columns:
        df[CA_DISTANCE_COLUMN] = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")

    return df


def compute_nh_f1_score(df: pd.DataFrame, predicted: pd.Series | None = None) -> TargetResult:
    """Compute F1 score using N/H CSP significance column."""
    nh_significant_col = "significant"
    if nh_significant_col not in df.columns:
        raise AlignmentParsingError(f"Missing {nh_significant_col} column")
    
    actual = df[nh_significant_col]
    if predicted is None:
        # Only use predictor columns that exist in the dataframe
        available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
        if not available_predictors:
            predicted = pd.Series([False] * len(df), index=df.index)
        else:
            predicted = df[available_predictors].any(axis=1)

    true_positives = int((actual & predicted).sum())
    false_positives = int((~actual & predicted).sum())
    false_negatives = int((actual & ~predicted).sum())

    denominator = 2 * true_positives + false_positives + false_negatives
    f1 = (2 * true_positives / denominator) if denominator else 0.0

    return TargetResult(
        target="",
        f1=f1,
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        total_rows=int(len(df)),
    )


def collect_nh_results(
    outputs_dir: Path, 
    allowed_targets: Set[str]
) -> List[TargetResult]:
    """Collect F1 scores from N/H CSP data for the same targets."""
    results: List[TargetResult] = []

    for target_dir in discover_targets(outputs_dir):
        # Filter by allowed targets
        if target_dir.name not in allowed_targets:
            continue
        
        alignment_path = target_dir / "master_alignment.csv"
        if not alignment_path.exists():
            continue

        try:
            df = load_nh_alignment(alignment_path)
        except AlignmentParsingError as exc:
            continue
        
        # Check required columns
        if "significant" not in df.columns:
            continue

        # Get available predictor columns
        available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
        if not available_predictors:
            continue
        
        predicted = df[available_predictors].any(axis=1)
        metrics = compute_nh_f1_score(df, predicted)
        metrics.target = target_dir.name
        results.append(metrics)

    return results


def collect_results(
    outputs_dir: Path, 
    allowed_targets: Set[str]
) -> Tuple[List[TargetResult], List[DistanceRecord], List[ConfusionRecord]]:
    results: List[TargetResult] = []
    distances: List[DistanceRecord] = []
    confusion_records: List[ConfusionRecord] = []

    for target_dir in discover_targets(outputs_dir):
        # Filter by allowed targets
        if target_dir.name not in allowed_targets:
            continue
        
        csp_table_path = target_dir / "csp_table_CA.csv"
        if not csp_table_path.exists():
            continue

        try:
            df = load_ca_alignment(target_dir)
        except AlignmentParsingError as exc:
            print(f"[WARN] Skipping {target_dir}: {exc}", file=sys.stderr)
            continue
        
        # Check required columns
        if SIGNIFICANT_COLUMN not in df.columns:
            print(f"[WARN] Skipping {target_dir}: missing {SIGNIFICANT_COLUMN} column", file=sys.stderr)
            continue

        # Get available predictor columns
        available_predictors = [col for col in PREDICTOR_COLUMNS if col in df.columns]
        if not available_predictors:
            print(f"[WARN] Skipping {target_dir}: no predictor columns found", file=sys.stderr)
            continue
        
        predicted = df[available_predictors].any(axis=1)
        metrics = compute_f1_score(df, predicted)
        metrics.target = target_dir.name
        results.append(metrics)

        # Collect distance records for significant residues
        if CA_DISTANCE_COLUMN in df.columns:
            significant_mask = df[SIGNIFICANT_COLUMN] & df[CA_DISTANCE_COLUMN].notna()
            if significant_mask.any():
                selected_distances = df.loc[significant_mask, CA_DISTANCE_COLUMN]
                selected_predictions = predicted.loc[significant_mask]
                for distance_value, predicted_value in zip(
                    selected_distances.tolist(), selected_predictions.tolist()
                ):
                    distances.append(
                        DistanceRecord(
                            distance=float(distance_value),
                            is_predicted_positive=bool(predicted_value),
                        )
                    )

            ca_distance_mask = df[CA_DISTANCE_COLUMN].notna()
            if ca_distance_mask.any():
                subset = df.loc[ca_distance_mask]
                for idx, row in subset.iterrows():
                    confusion_records.append(
                        ConfusionRecord(
                            distance=float(row[CA_DISTANCE_COLUMN]),
                            is_binding=bool(predicted.loc[idx]),
                            is_significant=bool(row[SIGNIFICANT_COLUMN]),
                        )
                    )

    return results, distances, confusion_records


def render_heatmap(results: Sequence[TargetResult], output_image: Path) -> None:
    dataframe = pd.DataFrame([result.__dict__ for result in results]).set_index("target")
    dataframe = dataframe.sort_values("f1", ascending=False)
    heatmap_data = dataframe[["f1"]]

    plt.figure(figsize=(6, max(3, len(heatmap_data) * 0.3)))
    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "F1 Score"},
    )
    # set larger font sizes for labels and title
    ax.set_xlabel("Metric", fontsize=16)
    ax.set_ylabel("Target", fontsize=16)
    ax.set_title("Per-target F1 Scores (CA-inclusive CSPs)", fontsize=18)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label("F1 Score", fontsize=16)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_histogram(distance_records: Sequence[DistanceRecord], output_image: Path) -> None:
    distances = [record.distance for record in distance_records]
    if not distances:
        print("No significant residues with CA distances available; skipping histogram.", file=sys.stderr)
        return

    max_distance = max(distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=bins, edgecolor="black", color="#4c72b0")
    plt.xlabel("Minimum CA Distance (Å)")
    plt.ylabel("Number of Significant Residues")
    plt.title("Distribution of Significant Residues by Minimum CA Distance (CA-inclusive CSPs)")
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_stacked_histogram(
    distance_records: Sequence[DistanceRecord],
    output_image: Path,
    positive_count: int,
    negative_count: int,
    tp_color: Optional[str] = None,
    fp_color: Optional[str] = None,
    show_title: bool = True,
    ylabel: Optional[str] = None,
    bold_axes: bool = False,
    axis_label_fontsize: Optional[int] = None,
    legend_fontsize: Optional[int] = None,
) -> None:
    if not distance_records:
        print(
            "No significant residues with CA distances available; skipping stacked histogram.",
            file=sys.stderr,
        )
        return

    positive_distances = [record.distance for record in distance_records if record.is_predicted_positive]
    negative_distances = [record.distance for record in distance_records if not record.is_predicted_positive]

    if not positive_distances and not negative_distances:
        print(
            "No significant residues with CA distances available; skipping stacked histogram.",
            file=sys.stderr,
        )
        return

    tp_color = tp_color if tp_color is not None else classification_colors.TP
    fp_color = fp_color if fp_color is not None else classification_colors.FP

    max_distance = max(positive_distances + negative_distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.hist(
        [positive_distances, negative_distances],
        bins=bins,
        stacked=True,
        color=[tp_color, fp_color],
        edgecolor="black",
        label=[
            f"(TP) Sig. CSP in Binding Site ({positive_count})",
            f"(FP) Sig. CSP -- Allosteric ({negative_count})",
        ],
    )
    y_label = ylabel if ylabel is not None else "Number of Significant Residues"
    xlabel_kw = {"fontweight": "bold" if bold_axes else "normal"}
    ylabel_kw = {"fontweight": "bold" if bold_axes else "normal"}
    if axis_label_fontsize is not None:
        xlabel_kw["fontsize"] = axis_label_fontsize
        ylabel_kw["fontsize"] = axis_label_fontsize
    ax.set_xlabel("Minimum CA Distance (Å)", **xlabel_kw)
    ax.set_ylabel(y_label, **ylabel_kw)
    if bold_axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")
    if show_title:
        plt.title("Predicted Outcomes for Significant Residues by Minimum CA Distance (CA-inclusive CSPs)")
    legend_kw = {}
    if legend_fontsize is not None:
        legend_kw["fontsize"] = legend_fontsize
    plt.legend(**legend_kw)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_f1_comparison_scatterplot(
    ca_results: Sequence[TargetResult],
    nh_results: Sequence[TargetResult],
    output_image: Path,
) -> None:
    """
    Create a scatterplot comparing F1 scores from CA-inclusive vs N/H CSPs.
    
    Args:
        ca_results: F1 scores from CA-inclusive CSPs
        nh_results: F1 scores from N/H CSPs
        output_image: Path to save the plot
    """
    # Create dictionaries for easy lookup
    ca_f1_dict = {result.target: result.f1 for result in ca_results}
    nh_f1_dict = {result.target: result.f1 for result in nh_results}
    
    # Find common targets
    common_targets = set(ca_f1_dict.keys()) & set(nh_f1_dict.keys())
    
    if not common_targets:
        print("No common targets found between CA and N/H CSP data; skipping scatterplot.", file=sys.stderr)
        return
    
    # Extract F1 scores for common targets
    ca_f1_scores = [ca_f1_dict[target] for target in sorted(common_targets)]
    nh_f1_scores = [nh_f1_dict[target] for target in sorted(common_targets)]
    target_labels = sorted(common_targets)
    
    # Calculate average delta: Average (N/H F1 - CA F1)
    deltas = [nh - ca for nh, ca in zip(nh_f1_scores, ca_f1_scores)]
    average_delta = sum(deltas) / len(deltas) if deltas else 0.0
    
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    
    # Create scatter plot
    ax.scatter(ca_f1_scores, nh_f1_scores, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    # Add y=x line (dashed red)
    ax.plot([0, 1], [0, 1], 'r--', linewidth=2, label='y=x', alpha=0.7)
    
    # Set axis limits and labels
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xlabel("F1 Score (CA-inclusive CSPs)", fontsize=14, fontweight='bold')
    ax.set_ylabel("F1 Score (N/H CSPs)", fontsize=14, fontweight='bold')
    ax.set_title("Comparison of F1 Scores: CA-inclusive vs N/H CSPs", fontsize=16, fontweight='bold')
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_aspect('equal', adjustable='box')
    
    # Add legend
    ax.legend(loc='lower right', fontsize=12)
    
    # Add text annotation with number of points and average delta
    annotation_text = f"n = {len(common_targets)} targets\n"
    annotation_text += f"Avg Δ = {average_delta:.4f}\n"
    annotation_text += "(N/H F1 - CA F1)"
    ax.text(
        0.05,
        0.95,
        annotation_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8),
    )
    
    # Print average delta to console
    print(f"Average delta (N/H F1 - CA F1): {average_delta:.4f}")
    
    plt.tight_layout()
    
    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_confusion_matrix_stacked_histogram(
    confusion_records: Sequence[ConfusionRecord],
    output_image: Path,
    tn_color: Optional[str] = None,
    fp_color: Optional[str] = None,
    fn_color: Optional[str] = None,
    tp_color: Optional[str] = None,
    show_title: bool = True,
    bold_axes: bool = False,
    axis_label_fontsize: Optional[int] = None,
    legend_fontsize: Optional[int] = None,
) -> None:
    """Render stacked histogram of TP, FP, FN, TN by minimum CA distance for CA-inclusive CSPs."""
    if not confusion_records:
        print(
            "No residues with CA distances available; skipping confusion matrix stacked histogram.",
            file=sys.stderr,
        )
        return

    tp_distances = [
        r.distance for r in confusion_records
        if r.is_significant and r.is_binding
    ]
    fp_distances = [
        r.distance for r in confusion_records
        if r.is_significant and not r.is_binding
    ]
    fn_distances = [
        r.distance for r in confusion_records
        if not r.is_significant and r.is_binding
    ]
    tn_distances = [
        r.distance for r in confusion_records
        if not r.is_significant and not r.is_binding
    ]

    all_distances = tp_distances + fp_distances + fn_distances + tn_distances
    if not all_distances:
        print(
            "No residues with CA distances available; skipping confusion matrix stacked histogram.",
            file=sys.stderr,
        )
        return

    tn_color = tn_color if tn_color is not None else classification_colors.TN
    fp_color = fp_color if fp_color is not None else classification_colors.FP
    fn_color = fn_color if fn_color is not None else classification_colors.FN
    tp_color = tp_color if tp_color is not None else classification_colors.TP

    max_distance = max(all_distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    colors = [tn_color, fp_color, fn_color, tp_color]  # TN, FP, FN, TP
    labels = [
        f"(TN) Small CSP -- allosteric ({len(tn_distances)})",
        f"(FP) Sig. CSP -- allosteric ({len(fp_distances)})",
        f"(FN) Small CSP in Binding Site ({len(fn_distances)})",
        f"(TP) Sig. CSP in Binding Site ({len(tp_distances)})",
    ]

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.hist(
        [tn_distances, fp_distances, fn_distances, tp_distances],
        bins=bins,
        stacked=True,
        color=colors,
        edgecolor="black",
        label=labels,
    )
    xlabel_kw = {"fontweight": "bold" if bold_axes else "normal"}
    ylabel_kw = {"fontweight": "bold" if bold_axes else "normal"}
    if axis_label_fontsize is not None:
        xlabel_kw["fontsize"] = axis_label_fontsize
        ylabel_kw["fontsize"] = axis_label_fontsize
    ax.set_xlabel("Minimum CA Distance (Å)", **xlabel_kw)
    ax.set_ylabel("Number of Residues", **ylabel_kw)
    if bold_axes:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")
    if show_title:
        plt.title("Confusion Matrix by Minimum CA Distance (CA-inclusive CSPs)")
    legend_kw = {"title": "Confusion Matrix"}
    if legend_fontsize is not None:
        legend_kw["fontsize"] = legend_fontsize
        legend_kw["title_fontsize"] = legend_fontsize
    plt.legend(**legend_kw)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir.resolve()
    output_image = args.output_image.resolve()
    summary_csv = args.summary_csv.resolve() if args.summary_csv else None
    histogram_image = args.histogram_image.resolve()
    stacked_histogram_image = args.stacked_histogram_image.resolve()
    confusion_matrix_stacked_histogram_image = (
        args.confusion_matrix_stacked_histogram_image.resolve()
    )
    scatterplot_image = args.scatterplot_image.resolve()

    # Load allowed targets (required)
    allowed_targets = load_targets_from_csv(args.targets_csv.resolve())
    print(f"Processing {len(allowed_targets)} targets from {args.targets_csv}")

    # Collect CA CSP results
    results, distances, confusion_records = collect_results(outputs_dir, allowed_targets)
    
    # Collect N/H CSP results for comparison
    nh_results = collect_nh_results(outputs_dir, allowed_targets)
    if not results:
        print(
            f"No valid CA CSP data found for targets in {args.targets_csv}",
            file=sys.stderr,
        )
        return 1

    if summary_csv:
        summary_df = pd.DataFrame([result.__dict__ for result in results]).sort_values("target")
        summary_csv.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_csv, index=False)

    positive_count = sum(1 for record in distances if record.is_predicted_positive)
    negative_count = sum(1 for record in distances if not record.is_predicted_positive)

    print(f"Significant CA CSPs predicted significant: {positive_count}")
    print(f"Significant CA CSPs missed: {negative_count}")

    render_heatmap(results, output_image)
    render_histogram(distances, histogram_image)
    render_stacked_histogram(distances, stacked_histogram_image, positive_count, negative_count)
    render_confusion_matrix_stacked_histogram(
        confusion_records, confusion_matrix_stacked_histogram_image
    )
    render_f1_comparison_scatterplot(results, nh_results, scatterplot_image)
    print(f"Heatmap saved to {output_image}")
    print(f"Histogram saved to {histogram_image}")
    print(f"Stacked histogram saved to {stacked_histogram_image}")
    print(f"Confusion matrix stacked histogram saved to {confusion_matrix_stacked_histogram_image}")
    print(f"F1 comparison scatterplot saved to {scatterplot_image}")
    if summary_csv:
        print(f"Summary table written to {summary_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
