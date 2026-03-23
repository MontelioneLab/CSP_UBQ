#!/usr/bin/env python3
"""
Analyze CA-inclusive CSP outputs for targets with matching last authors.

This script is similar to analyze_targets_ca.py but:
1. Identifies targets from CSP_UBQ.csv where apo and holo citation info have the same last author
2. Uses csp_table_CA.csv as the base (instead of master_alignment.csv)
3. Uses csp_CA_significant instead of significant
4. Merges CA CSP data with other analysis files (occlusion, interaction, distance)
"""

from __future__ import annotations

import argparse
import csv
import os
import re
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
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.align import align_global


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


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-target F1 scores from CA-inclusive CSP data (csp_table_CA.csv) "
            "for targets where apo and holo citation info have the same last author. "
            "Targets are identified from CSP_UBQ.csv."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--csprank-csv",
        type=Path,
        default=Path("CSP_UBQ.csv"),
        help="Path to CSP_UBQ.csv file (default: %(default)s).",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("outputs") / "f1_heatmap_same_author.png",
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
        default=Path("outputs") / "significant_ca_distance_hist_same_author.png",
        help=(
            "Destination for the histogram of significant residues by CA distance "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--stacked-histogram-image",
        type=Path,
        default=Path("outputs") / "significant_ca_distance_stacked_hist_same_author.png",
        help=(
            "Destination for the stacked histogram splitting predicted positives "
            "and negatives among significant residues (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--scatterplot-image",
        type=Path,
        default=Path("outputs") / "f1_comparison_scatterplot_same_author.png",
        help=(
            "Destination for the scatterplot comparing F1 scores from CA-inclusive "
            "vs N/H CSPs (default: %(default)s)."
        ),
    )
    return parser.parse_args(list(argv))


def extract_last_author(citation_text: str) -> Optional[str]:
    """Extract the last author from citation text."""
    if not citation_text:
        return None
    
    # Find the Authors line
    authors_match = re.search(r'Authors:\s*(.+?)(?:\n|$)', citation_text, re.MULTILINE)
    if not authors_match:
        return None
    
    authors_line = authors_match.group(1).strip()
    
    # Split by comma and get the last author
    authors = [a.strip() for a in authors_line.split(',')]
    if authors:
        return authors[-1]
    return None


def load_targets_with_same_last_author(csprank_path: Path) -> Set[str]:
    """
    Load target IDs from CSP_UBQ.csv where apo and holo citation info have the same last author.
    
    Args:
        csprank_path: Path to CSP_UBQ.csv file
        
    Returns:
        Set of target IDs (holo_pdb values) with matching last authors
    """
    if not csprank_path.exists():
        raise FileNotFoundError(f"CSP_UBQ.csv not found: {csprank_path}")
    
    targets: Set[str] = set()
    
    try:
        with open(csprank_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                holo_pdb = row.get('holo_pdb', '').strip()
                apo_citation = row.get('apo_citation_info', '')
                holo_citation = row.get('holo_citation_info', '')
                
                # Skip if no holo_pdb
                if not holo_pdb:
                    continue
                
                # Extract last authors
                apo_last_author = extract_last_author(apo_citation)
                holo_last_author = extract_last_author(holo_citation)
                
                # Check if both exist and match
                if apo_last_author and holo_last_author:
                    if apo_last_author == holo_last_author:
                        targets.add(holo_pdb)
    except Exception as exc:
        raise ValueError(f"Failed to load targets from {csprank_path}: {exc}") from exc
    
    return targets


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
) -> Tuple[List[TargetResult], List[DistanceRecord]]:
    results: List[TargetResult] = []
    distances: List[DistanceRecord] = []

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

    return results, distances


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
    ax.set_title("Per-target F1 Scores (CA-inclusive CSPs, Same Last Author)", fontsize=18)
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
    plt.title("Distribution of Significant Residues by Minimum CA Distance (CA-inclusive CSPs, Same Last Author)")
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_stacked_histogram(
    distance_records: Sequence[DistanceRecord],
    output_image: Path,
    positive_count: int,
    negative_count: int,
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

    max_distance = max(positive_distances + negative_distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.hist(
        [positive_distances, negative_distances],
        bins=bins,
        stacked=True,
        color=["#2ca02c", "#d62728"],
        edgecolor="black",
        label=["Predicted Positive", "Predicted Negative"],
    )
    plt.xlabel("Minimum CA Distance (Å)")
    plt.ylabel("Number of Significant Residues")
    plt.title("Predicted Outcomes for Significant Residues by Minimum CA Distance (CA-inclusive CSPs, Same Last Author)")
    plt.legend()
    ax.text(
        0.98,
        0.95,
        (
            f"Predicted significant: {positive_count}\n"
            f"Missed significant: {negative_count}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85),
    )
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
    ax.set_title("Comparison of F1 Scores: CA-inclusive vs N/H CSPs (Same Last Author)", fontsize=16, fontweight='bold')
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_aspect('equal', adjustable='box')
    
    # Add legend
    ax.legend(loc='lower right', fontsize=12)
    
    # Add text annotation with number of points
    ax.text(
        0.05,
        0.95,
        f"n = {len(common_targets)} targets",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8),
    )
    
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
    scatterplot_image = args.scatterplot_image.resolve()

    # Load targets with matching last authors from CSP_UBQ.csv
    allowed_targets = load_targets_with_same_last_author(args.csprank_csv.resolve())
    print(f"Found {len(allowed_targets)} targets with matching last authors from {args.csprank_csv}")
    if allowed_targets:
        print(f"Targets: {', '.join(sorted(allowed_targets))}")

    # Collect CA CSP results
    results, distances = collect_results(outputs_dir, allowed_targets)
    
    # Collect N/H CSP results for comparison
    nh_results = collect_nh_results(outputs_dir, allowed_targets)
    
    if not results:
        print(
            f"No valid CA CSP data found for targets with matching last authors",
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
    render_f1_comparison_scatterplot(results, nh_results, scatterplot_image)
    print(f"Heatmap saved to {output_image}")
    print(f"Histogram saved to {histogram_image}")
    print(f"Stacked histogram saved to {stacked_histogram_image}")
    print(f"F1 comparison scatterplot saved to {scatterplot_image}")
    if summary_csv:
        print(f"Summary table written to {summary_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

