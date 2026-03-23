#!/usr/bin/env python3
"""
F1 Score Reporter - Calculate and report F1 scores for user-specified targets.

This script extends the core functionality of analyze_targets.py to allow users
to provide a comma-separated list of holo PDB IDs and receive F1 score statistics
for each target along with summary statistics.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Set

import pandas as pd

# Import core functions from analyze_targets.py
try:
    from .analyze_targets import (
        AlignmentParsingError,
        TargetResult,
        load_alignment,
        compute_f1_score,
        PREDICTOR_COLUMNS,
    )
except ImportError:
    from analyze_targets import (
        AlignmentParsingError,
        TargetResult,
        load_alignment,
        compute_f1_score,
        PREDICTOR_COLUMNS,
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate and report F1 scores for a user-specified subset of targets. "
            "Targets must exist in CSP_UBQ.csv."
        )
    )
    parser.add_argument(
        "--targets",
        type=str,
        required=True,
        help="Comma-separated list of holo PDB IDs (e.g., '1cf4,2lox,2mkr')",
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
        default=Path("data/CSP_UBQ.csv"),
        help="Path to CSP_UBQ.csv file (default: %(default)s).",
    )
    return parser.parse_args(argv)


def validate_targets(
    target_list: List[str], csprank_csv: Path
) -> tuple[Set[str], List[str]]:
    """
    Validate that all provided targets exist in CSP_UBQ.csv.
    
    Args:
        target_list: List of target PDB IDs to validate
        csprank_csv: Path to CSP_UBQ.csv file
        
    Returns:
        Tuple of (valid_targets_set, invalid_targets_list)
    """
    if not csprank_csv.exists():
        raise FileNotFoundError(f"CSP_UBQ.csv not found: {csprank_csv}")
    
    try:
        df = pd.read_csv(csprank_csv)
        if "holo_pdb" not in df.columns:
            raise ValueError(f"CSP_UBQ.csv must have a 'holo_pdb' column")
        
        # Get all valid targets from CSP_UBQ.csv
        valid_targets = set(df["holo_pdb"].astype(str).str.strip().str.lower())
        
        # Normalize input targets (lowercase, strip whitespace)
        normalized_input = [t.strip().lower() for t in target_list]
        
        # Find invalid targets
        invalid_targets = [t for t in normalized_input if t not in valid_targets]
        
        # Return valid targets (as original case from input) and invalid list
        valid_set = {t for t in normalized_input if t in valid_targets}
        
        return valid_set, invalid_targets
    except Exception as exc:
        raise ValueError(f"Failed to read CSP_UBQ.csv: {exc}") from exc


def calculate_f1_for_targets(
    targets: Set[str], outputs_dir: Path
) -> List[TargetResult]:
    """
    Calculate F1 scores for each target.
    
    Args:
        targets: Set of validated target PDB IDs
        outputs_dir: Root directory containing per-target subdirectories
        
    Returns:
        List of TargetResult objects
    """
    results: List[TargetResult] = []
    
    for target in sorted(targets):
        target_dir = outputs_dir / target
        alignment_path = target_dir / "master_alignment.csv"
        
        if not target_dir.exists():
            print(
                f"[WARN] Target directory not found: {target_dir}. Skipping.",
                file=sys.stderr,
            )
            continue
        
        if not alignment_path.exists():
            print(
                f"[WARN] master_alignment.csv not found for {target}. Skipping.",
                file=sys.stderr,
            )
            continue
        
        try:
            df = load_alignment(alignment_path)
            predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)
            metrics = compute_f1_score(df, predicted)
            metrics.target = target
            results.append(metrics)
        except AlignmentParsingError as exc:
            print(
                f"[WARN] Failed to parse {alignment_path}: {exc}. Skipping.",
                file=sys.stderr,
            )
            continue
        except Exception as exc:
            print(
                f"[WARN] Unexpected error processing {target}: {exc}. Skipping.",
                file=sys.stderr,
            )
            continue
    
    return results


def print_results_table(results: List[TargetResult]) -> None:
    """
    Print formatted table of F1 scores and summary statistics.
    
    Args:
        results: List of TargetResult objects
    """
    if not results:
        print("No results to display.")
        return
    
    # Create DataFrame for individual results
    df = pd.DataFrame([result.__dict__ for result in results])
    df = df.sort_values("f1", ascending=False)
    
    # Format the table
    print("\n" + "=" * 70)
    print("F1 Score Report")
    print("=" * 70)
    print("\nTarget Results:")
    print("-" * 70)
    
    # Format individual results table
    display_df = df[["target", "f1", "true_positives", "false_positives", "false_negatives", "total_rows"]].copy()
    display_df.columns = ["Target", "F1 Score", "TP", "FP", "FN", "Total Rows"]
    
    # Format F1 score to 3 decimal places
    display_df["F1 Score"] = display_df["F1 Score"].map("{:.3f}".format)
    
    # Print table using pandas formatting
    print(display_df.to_string(index=False))
    print("-" * 70)
    
    # Calculate and print summary statistics
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
    """Main entry point."""
    args = parse_args(argv)
    
    # Parse comma-separated targets
    target_list = [t.strip() for t in args.targets.split(",") if t.strip()]
    
    if not target_list:
        print("Error: No targets provided. Use --targets with comma-separated PDB IDs.", file=sys.stderr)
        return 1
    
    # Validate targets
    try:
        valid_targets, invalid_targets = validate_targets(
            target_list, args.csprank_csv.resolve()
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    
    # Report invalid targets
    if invalid_targets:
        print(
            f"Warning: The following targets were not found in CSP_UBQ.csv: {', '.join(invalid_targets)}",
            file=sys.stderr,
        )
    
    if not valid_targets:
        print("Error: No valid targets found. Exiting.", file=sys.stderr)
        return 1
    
    # Calculate F1 scores
    outputs_dir = args.outputs_dir.resolve()
    if not outputs_dir.exists():
        print(f"Error: Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1
    
    results = calculate_f1_for_targets(valid_targets, outputs_dir)
    
    if not results:
        print(
            "Error: No results calculated. Check that master_alignment.csv files exist for the specified targets.",
            file=sys.stderr,
        )
        return 1
    
    # Print results
    print_results_table(results)
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

