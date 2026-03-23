"""
Confusion Matrix Analysis for CSP vs Occlusion and Interactions

This script analyzes the relationship between significant CSPs and positive residues
(defined as the union of occluded NH groups, hydrogen-bonded residues, and 
charge-complementary residues) by creating confusion matrices and calculating 
performance metrics. It uses the master_alignment.csv file for each target which
contains all necessary data in a single file.

Usage:
    python scripts/confusion_matrix_analysis.py --outputs outputs/
"""

from __future__ import annotations

import argparse
import csv
import os
import math
import statistics
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

try:
    from .config import paths
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import paths


@dataclass
class ConfusionMatrix:
    """Confusion matrix components"""
    tp: int  # True Positives: significant CSP AND occluded NH
    fp: int  # False Positives: occluded NH but NOT significant CSP
    tn: int  # True Negatives: NOT significant CSP AND NOT occluded NH
    fn: int  # False Negatives: significant CSP but NOT occluded NH
    
    @property
    def total(self) -> int:
        return self.tp + self.fp + self.tn + self.fn
    
    @property
    def precision(self) -> float:
        """Precision = TP / (TP + FP)"""
        if self.tp + self.fp == 0:
            return 0.0
        return self.tp / (self.tp + self.fp)
    
    @property
    def recall(self) -> float:
        """Recall = TP / (TP + FN)"""
        if self.tp + self.fn == 0:
            return 0.0
        return self.tp / (self.tp + self.fn)
    
    @property
    def f1_score(self) -> float:
        """F1 = 2 * (Precision * Recall) / (Precision + Recall)"""
        p = self.precision
        r = self.recall
        if p + r == 0:
            return 0.0
        return 2 * (p * r) / (p + r)
    
    @property
    def mcc(self) -> float:
        """Matthews Correlation Coefficient"""
        numerator = self.tp * self.tn - self.fp * self.fn
        denominator = math.sqrt((self.tp + self.fp) * (self.tp + self.fn) * 
                              (self.tn + self.fp) * (self.tn + self.fn))
        if denominator == 0:
            return 0.0
        return numerator / denominator


@dataclass
class SystemResults:
    """Results for a single system"""
    system_id: str
    confusion_matrix: ConfusionMatrix
    total_residues: int
    significant_csp_count: int
    occluded_count: int
    proline_count: int
    interacting_count: int
    csp_threshold: str = 'significant'
    positive_strategy: str = 'occluded_only'


def read_master_alignment(csv_path: str) -> List[Dict[str, str]]:
    """Read master alignment CSV and return list of rows"""
    rows = []
    try:
        with open(csv_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
    except FileNotFoundError:
        print(f"Warning: Could not find {csv_path}")
    except Exception as e:
        print(f"Warning: Error reading {csv_path}: {e}")
    return rows




def compute_confusion_matrix(master_alignment_rows: List[Dict[str, str]], significant_column: str = 'significant', positive_strategy: str = 'occluded_only') -> Tuple[ConfusionMatrix, int, int, int, int, int]:
    """
    Compute confusion matrix from master alignment CSV data.
    
    Args:
        master_alignment_rows: Master alignment CSV rows containing all necessary data
        significant_column: Column name for CSP significance ('significant', 'significant_1sd', 'significant_2sd')
        positive_strategy: Strategy for defining positive residues ('occluded_only', 'occluded_or_ca', 'occluded_or_interaction', 'occluded_or_ca_or_interaction')
        
    Returns:
        (confusion_matrix, total_residues, significant_csp_count, occluded_count, proline_count, interacting_count)
    """
    tp = fp = tn = fn = 0
    total_residues = 0
    significant_csp_count = 0
    occluded_count = 0
    interacting_count = 0
    proline_count = 0
    
    for row in master_alignment_rows:
        # Skip rows with missing essential data
        significant_str = row.get(significant_column, '').strip()
        holo_aa = row.get('holo_aa', '').strip()
        
        if not significant_str or not holo_aa:
            continue
            
        try:
            significant = int(significant_str) == 1
        except ValueError:
            continue
        
        # Count prolines (residues without NH groups)
        if holo_aa == 'P':
            proline_count += 1
            continue
        
        # Get occlusion status
        is_occluded_str = row.get('is_occluded_occlusion', '').strip()
        is_occluded = is_occluded_str.lower() == 'true' if is_occluded_str else False
        
        # Get interaction status
        has_hbond_str = row.get('has_hbond_interaction', '').strip()
        has_charge_complement_str = row.get('has_charge_complement_interaction', '').strip()
        has_pi_contact_str = row.get('has_pi_contact_interaction', '').strip()
        has_hbond = has_hbond_str.lower() == 'true' if has_hbond_str else False
        has_charge_complement = has_charge_complement_str.lower() == 'true' if has_charge_complement_str else False
        has_pi_contact = has_pi_contact_str.lower() == 'true' if has_pi_contact_str else False
        is_interacting = has_hbond or has_charge_complement or has_pi_contact
        
        # Get CA distance status
        passes_ca_filter_str = row.get('passes_filter_distance', '').strip()
        passes_ca_filter = passes_ca_filter_str.lower() == 'true' if passes_ca_filter_str else False
        
        # Determine positive status based on strategy
        if positive_strategy == 'occluded_only':
            is_positive = is_occluded
        elif positive_strategy == 'occluded_or_ca':
            is_positive = is_occluded or passes_ca_filter
        elif positive_strategy == 'occluded_or_interaction':
            is_positive = is_occluded or is_interacting
        elif positive_strategy == 'occluded_or_ca_or_interaction':
            is_positive = is_occluded or passes_ca_filter or is_interacting
        else:
            # Default to occluded only
            is_positive = is_occluded
        
        total_residues += 1
        
        if significant:
            significant_csp_count += 1
        if is_occluded:
            occluded_count += 1
        if is_interacting:
            interacting_count += 1
            
        # Confusion matrix classification using selected positive strategy
        if significant and is_positive:
            tp += 1
        elif significant and not is_positive:
            fp += 1
        elif not significant and not is_positive:
            tn += 1
        elif not significant and is_positive:
            fn += 1
    
    print(f"[CONFUSION] TP={tp}, FP={fp}, TN={tn}, FN={fn}")
    print(f"[CONFUSION] Total residues: {total_residues}, Significant CSPs: {significant_csp_count}, Occluded: {occluded_count}, Interacting: {interacting_count}")
    
    return ConfusionMatrix(tp=tp, fp=fp, tn=tn, fn=fn), total_residues, significant_csp_count, occluded_count, proline_count, interacting_count


def analyze_system(system_dir: str, significant_column: str = 'significant', positive_strategy: str = 'occluded_only') -> Optional[SystemResults]:
    """Analyze a single system directory"""
    system_id = os.path.basename(system_dir)
    master_alignment_path = os.path.join(system_dir, "master_alignment.csv")
    
    if not os.path.exists(master_alignment_path):
        print(f"Skipping {system_id}: no master_alignment.csv found")
        return None
    
    master_alignment_rows = read_master_alignment(master_alignment_path)
    
    if not master_alignment_rows:
        print(f"Skipping {system_id}: empty or invalid master_alignment.csv")
        return None
    
    print(f"[ANALYZE] Processing {system_id}: {len(master_alignment_rows)} master alignment rows")
    
    confusion_matrix, total_residues, significant_csp_count, occluded_count, proline_count, interacting_count = compute_confusion_matrix(master_alignment_rows, significant_column, positive_strategy)
    
    return SystemResults(
        system_id=system_id,
        confusion_matrix=confusion_matrix,
        total_residues=total_residues,
        significant_csp_count=significant_csp_count,
        occluded_count=occluded_count,
        proline_count=proline_count,
        interacting_count=interacting_count,
        csp_threshold=significant_column,
        positive_strategy=positive_strategy
    )


def write_per_system_results(results: List[SystemResults], output_path: str) -> None:
    """Write per-system results to CSV"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'system_id', 'tp', 'fp', 'tn', 'fn', 'total_residues',
            'significant_csp_count', 'occluded_count', 'interacting_count', 'proline_count',
            'precision', 'recall', 'f1_score', 'mcc'
        ])
        
        for result in results:
            cm = result.confusion_matrix
            writer.writerow([
                result.system_id,
                cm.tp, cm.fp, cm.tn, cm.fn,
                result.total_residues,
                result.significant_csp_count,
                result.occluded_count,
                result.interacting_count,
                result.proline_count,
                f"{cm.precision:.4f}",
                f"{cm.recall:.4f}",
                f"{cm.f1_score:.4f}",
                f"{cm.mcc:.4f}"
            ])


def compute_aggregated_statistics(results: List[SystemResults]) -> Dict[str, float]:
    """Compute aggregated statistics across all systems"""
    if not results:
        return {}
    
    # Aggregate confusion matrix
    total_tp = sum(r.confusion_matrix.tp for r in results)
    total_fp = sum(r.confusion_matrix.fp for r in results)
    total_tn = sum(r.confusion_matrix.tn for r in results)
    total_fn = sum(r.confusion_matrix.fn for r in results)
    
    aggregated_cm = ConfusionMatrix(tp=total_tp, fp=total_fp, tn=total_tn, fn=total_fn)
    
    # Per-system metrics
    precisions = [r.confusion_matrix.precision for r in results]
    recalls = [r.confusion_matrix.recall for r in results]
    f1_scores = [r.confusion_matrix.f1_score for r in results]
    mccs = [r.confusion_matrix.mcc for r in results]
    
    return {
        'total_systems': len(results),
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_tn': total_tn,
        'total_fn': total_fn,
        'total_residues': sum(r.total_residues for r in results),
        'total_significant_csp': sum(r.significant_csp_count for r in results),
        'total_occluded': sum(r.occluded_count for r in results),
        'total_interacting': sum(r.interacting_count for r in results),
        'total_prolines': sum(r.proline_count for r in results),
        
        # Aggregated metrics
        'aggregated_precision': aggregated_cm.precision,
        'aggregated_recall': aggregated_cm.recall,
        'aggregated_f1_score': aggregated_cm.f1_score,
        'aggregated_mcc': aggregated_cm.mcc,
        
        # Per-system statistics
        'mean_precision': statistics.mean(precisions),
        'median_precision': statistics.median(precisions),
        'std_precision': statistics.stdev(precisions) if len(precisions) > 1 else 0.0,
        
        'mean_recall': statistics.mean(recalls),
        'median_recall': statistics.median(recalls),
        'std_recall': statistics.stdev(recalls) if len(recalls) > 1 else 0.0,
        
        'mean_f1_score': statistics.mean(f1_scores),
        'median_f1_score': statistics.median(f1_scores),
        'std_f1_score': statistics.stdev(f1_scores) if len(f1_scores) > 1 else 0.0,
        
        'mean_mcc': statistics.mean(mccs),
        'median_mcc': statistics.median(mccs),
        'std_mcc': statistics.stdev(mccs) if len(mccs) > 1 else 0.0,
    }


def run_grid_search(system_dirs: List[str], verbose: bool = False) -> Tuple[List[Dict[str, any]], List[List[SystemResults]]]:
    """
    Run grid search over all parameter combinations.
    
    Args:
        system_dirs: List of system directory paths
        verbose: Enable verbose output
        
    Returns:
        Tuple of (grid_results, all_system_results) where:
        - grid_results: List of dictionaries containing aggregated results for each parameter combination
        - all_system_results: List of lists, where each inner list contains SystemResults for one parameter combination
    """
    # Define parameter combinations
    csp_thresholds = ['significant', 'significant_1sd', 'significant_2sd']
    positive_strategies = ['occluded_only', 'occluded_or_ca', 'occluded_or_interaction', 'occluded_or_ca_or_interaction']
    
    grid_results = []
    all_system_results = []
    
    for csp_threshold in csp_thresholds:
        for positive_strategy in positive_strategies:
            if verbose:
                print(f"\n[GRID] Testing CSP threshold: {csp_threshold}, Positive strategy: {positive_strategy}")
            
            # Analyze all systems with this parameter combination
            results = []
            for system_dir in sorted(system_dirs):
                if verbose:
                    print(f"Analyzing {os.path.basename(system_dir)}...")
                
                result = analyze_system(system_dir, csp_threshold, positive_strategy)
                if result:
                    results.append(result)
                    if verbose:
                        cm = result.confusion_matrix
                        print(f"  TP={cm.tp}, FP={cm.fp}, TN={cm.tn}, FN={cm.fn}")
                        print(f"  Precision={cm.precision:.3f}, Recall={cm.recall:.3f}, F1={cm.f1_score:.3f}, MCC={cm.mcc:.3f}")
            
            if not results:
                print(f"Warning: No valid systems found for {csp_threshold} + {positive_strategy}")
                continue
            
            # Store per-system results for this combination
            all_system_results.append(results)
            
            # Compute aggregated statistics for this combination
            stats = compute_aggregated_statistics(results)
            
            # Add parameter information to stats
            stats['csp_threshold'] = csp_threshold
            stats['positive_strategy'] = positive_strategy
            
            grid_results.append(stats)
            
            if verbose:
                print(f"Combination {csp_threshold} + {positive_strategy}: {len(results)} systems")
                print(f"  Aggregated Precision: {stats['aggregated_precision']:.3f}")
                print(f"  Aggregated Recall: {stats['aggregated_recall']:.3f}")
                print(f"  Aggregated F1: {stats['aggregated_f1_score']:.3f}")
                print(f"  Aggregated MCC: {stats['aggregated_mcc']:.3f}")
    
    return grid_results, all_system_results


def write_per_system_grid_results(all_system_results: List[List[SystemResults]], output_path: str) -> None:
    """Write per-system grid search results to CSV"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow([
            'system_id', 'csp_threshold', 'positive_strategy', 'tp', 'fp', 'tn', 'fn',
            'total_residues', 'significant_csp_count', 'occluded_count', 'interacting_count', 'proline_count',
            'precision', 'recall', 'f1_score', 'mcc'
        ])
        
        # Write data rows
        for system_results in all_system_results:
            for result in system_results:
                cm = result.confusion_matrix
                writer.writerow([
                    result.system_id,
                    result.csp_threshold,
                    result.positive_strategy,
                    cm.tp, cm.fp, cm.tn, cm.fn,
                    result.total_residues,
                    result.significant_csp_count,
                    result.occluded_count,
                    result.interacting_count,
                    result.proline_count,
                    f"{cm.precision:.4f}",
                    f"{cm.recall:.4f}",
                    f"{cm.f1_score:.4f}",
                    f"{cm.mcc:.4f}"
                ])


def write_grid_search_results(grid_results: List[Dict[str, any]], output_path: str) -> None:
    """Write grid search results to CSV"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow([
            'csp_threshold', 'positive_strategy', 'total_systems', 'total_tp', 'total_fp', 'total_tn', 'total_fn',
            'total_residues', 'total_significant_csp', 'total_occluded', 'total_interacting', 'total_prolines',
            'aggregated_precision', 'aggregated_recall', 'aggregated_f1_score', 'aggregated_mcc',
            'mean_precision', 'mean_recall', 'mean_f1_score', 'mean_mcc',
            'median_precision', 'median_recall', 'median_f1_score', 'median_mcc',
            'std_precision', 'std_recall', 'std_f1_score', 'std_mcc'
        ])
        
        # Write data rows
        for stats in grid_results:
            writer.writerow([
                stats['csp_threshold'],
                stats['positive_strategy'],
                f"{stats['total_systems']:.0f}",
                f"{stats['total_tp']:.0f}",
                f"{stats['total_fp']:.0f}",
                f"{stats['total_tn']:.0f}",
                f"{stats['total_fn']:.0f}",
                f"{stats['total_residues']:.0f}",
                f"{stats['total_significant_csp']:.0f}",
                f"{stats['total_occluded']:.0f}",
                f"{stats['total_interacting']:.0f}",
                f"{stats['total_prolines']:.0f}",
                f"{stats['aggregated_precision']:.4f}",
                f"{stats['aggregated_recall']:.4f}",
                f"{stats['aggregated_f1_score']:.4f}",
                f"{stats['aggregated_mcc']:.4f}",
                f"{stats['mean_precision']:.4f}",
                f"{stats['mean_recall']:.4f}",
                f"{stats['mean_f1_score']:.4f}",
                f"{stats['mean_mcc']:.4f}",
                f"{stats['median_precision']:.4f}",
                f"{stats['median_recall']:.4f}",
                f"{stats['median_f1_score']:.4f}",
                f"{stats['median_mcc']:.4f}",
                f"{stats['std_precision']:.4f}",
                f"{stats['std_recall']:.4f}",
                f"{stats['std_f1_score']:.4f}",
                f"{stats['std_mcc']:.4f}"
            ])


def write_summary_results(stats: Dict[str, float], output_path: str) -> None:
    """Write aggregated summary results to CSV"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow(['metric', 'value'])
        
        # Write aggregated confusion matrix
        writer.writerow(['total_systems', f"{stats['total_systems']:.0f}"])
        writer.writerow(['total_tp', f"{stats['total_tp']:.0f}"])
        writer.writerow(['total_fp', f"{stats['total_fp']:.0f}"])
        writer.writerow(['total_tn', f"{stats['total_tn']:.0f}"])
        writer.writerow(['total_fn', f"{stats['total_fn']:.0f}"])
        writer.writerow(['total_residues', f"{stats['total_residues']:.0f}"])
        writer.writerow(['total_significant_csp', f"{stats['total_significant_csp']:.0f}"])
        writer.writerow(['total_occluded', f"{stats['total_occluded']:.0f}"])
        writer.writerow(['total_interacting', f"{stats['total_interacting']:.0f}"])
        writer.writerow(['total_prolines', f"{stats['total_prolines']:.0f}"])
        
        # Write aggregated metrics
        writer.writerow(['aggregated_precision', f"{stats['aggregated_precision']:.4f}"])
        writer.writerow(['aggregated_recall', f"{stats['aggregated_recall']:.4f}"])
        writer.writerow(['aggregated_f1_score', f"{stats['aggregated_f1_score']:.4f}"])
        writer.writerow(['aggregated_mcc', f"{stats['aggregated_mcc']:.4f}"])
        
        # Write per-system statistics
        writer.writerow(['mean_precision', f"{stats['mean_precision']:.4f}"])
        writer.writerow(['median_precision', f"{stats['median_precision']:.4f}"])
        writer.writerow(['std_precision', f"{stats['std_precision']:.4f}"])
        
        writer.writerow(['mean_recall', f"{stats['mean_recall']:.4f}"])
        writer.writerow(['median_recall', f"{stats['median_recall']:.4f}"])
        writer.writerow(['std_recall', f"{stats['std_recall']:.4f}"])
        
        writer.writerow(['mean_f1_score', f"{stats['mean_f1_score']:.4f}"])
        writer.writerow(['median_f1_score', f"{stats['median_f1_score']:.4f}"])
        writer.writerow(['std_f1_score', f"{stats['std_f1_score']:.4f}"])
        
        writer.writerow(['mean_mcc', f"{stats['mean_mcc']:.4f}"])
        writer.writerow(['median_mcc', f"{stats['median_mcc']:.4f}"])
        writer.writerow(['std_mcc', f"{stats['std_mcc']:.4f}"])


def generate_confusion_matrix_per_system(outputs_dir: str, verbose: bool = False) -> bool:
    """
    Generate confusion_matrix_per_system.csv from all system directories.
    Uses default parameters (significant, occluded_only).
    Returns True if at least one system was analyzed and the file was written.
    """
    if not os.path.exists(outputs_dir):
        return False
    system_dirs = [
        os.path.join(outputs_dir, item)
        for item in os.listdir(outputs_dir)
        if os.path.isdir(os.path.join(outputs_dir, item))
    ]
    if not system_dirs:
        return False
    results = []
    for system_dir in sorted(system_dirs):
        if verbose:
            print(f"[CONFUSION] Analyzing {os.path.basename(system_dir)}...")
        result = analyze_system(system_dir)
        if result:
            results.append(result)
    if not results:
        return False
    per_system_path = os.path.join(outputs_dir, "confusion_matrix_per_system.csv")
    write_per_system_results(results, per_system_path)
    if verbose:
        print(f"[CONFUSION] Wrote {per_system_path} ({len(results)} systems)")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Confusion Matrix Analysis for CSP vs Occlusion")
    parser.add_argument("--outputs", default=paths.outputs_dir,
                       help="Path to outputs directory containing system subdirectories")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--grid-search", action="store_true",
                       help="Run grid search over parameter combinations")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.outputs):
        print(f"Error: Outputs directory {args.outputs} does not exist")
        return
    
    print(f"Analyzing systems in {args.outputs}")
    
    # Find all system directories
    system_dirs = []
    for item in os.listdir(args.outputs):
        item_path = os.path.join(args.outputs, item)
        if os.path.isdir(item_path):
            system_dirs.append(item_path)
    
    if not system_dirs:
        print("No system directories found")
        return
    
    print(f"Found {len(system_dirs)} system directories")
    
    if args.grid_search:
        # Run grid search over all parameter combinations
        print("Running grid search over parameter combinations...")
        grid_results, all_system_results = run_grid_search(system_dirs, args.verbose)
        
        if not grid_results:
            print("No valid results from grid search")
            return
        
        # Write grid search results
        grid_search_path = os.path.join(args.outputs, "confusion_matrix_grid_search.csv")
        write_grid_search_results(grid_results, grid_search_path)
        print(f"Grid search results written to {grid_search_path}")
        
        # Write per-system grid results
        per_system_grid_path = os.path.join(args.outputs, "confusion_matrix_per_system_grid.csv")
        write_per_system_grid_results(all_system_results, per_system_grid_path)
        print(f"Per-system grid results written to {per_system_grid_path}")
        
        # Print summary of best performing combinations
        print("\nGrid Search Summary:")
        print("Top 3 combinations by F1 score:")
        sorted_results = sorted(grid_results, key=lambda x: x['aggregated_f1_score'], reverse=True)
        for i, result in enumerate(sorted_results[:3]):
            print(f"{i+1}. {result['csp_threshold']} + {result['positive_strategy']}: F1={result['aggregated_f1_score']:.3f}, Precision={result['aggregated_precision']:.3f}, Recall={result['aggregated_recall']:.3f}")
    
    else:
        # Original single analysis (for backward compatibility)
        results = []
        for system_dir in sorted(system_dirs):
            if args.verbose:
                print(f"Analyzing {os.path.basename(system_dir)}...")
            
            result = analyze_system(system_dir)
            if result:
                results.append(result)
                if args.verbose:
                    cm = result.confusion_matrix
                    print(f"  TP={cm.tp}, FP={cm.fp}, TN={cm.tn}, FN={cm.fn}")
                    print(f"  Precision={cm.precision:.3f}, Recall={cm.recall:.3f}, F1={cm.f1_score:.3f}, MCC={cm.mcc:.3f}")
        
        if not results:
            print("No valid systems found for analysis")
            return
        
        print(f"Successfully analyzed {len(results)} systems")
        
        # Write per-system results
        per_system_path = os.path.join(args.outputs, "confusion_matrix_per_system.csv")
        write_per_system_results(results, per_system_path)
        print(f"Per-system results written to {per_system_path}")
        
        # Compute and write aggregated statistics
        stats = compute_aggregated_statistics(results)
        summary_path = os.path.join(args.outputs, "confusion_matrix_summary.csv")
        write_summary_results(stats, summary_path)
        print(f"Summary statistics written to {summary_path}")
        
        # Print summary
        print("\nSummary:")
        print(f"Total systems: {stats['total_systems']}")
        print(f"Total residues: {stats['total_residues']}")
        print(f"Total significant CSPs: {stats['total_significant_csp']}")
        print(f"Total occluded: {stats['total_occluded']}")
        print(f"Total interacting: {stats['total_interacting']}")
        print(f"Aggregated Precision: {stats['aggregated_precision']:.3f}")
        print(f"Aggregated Recall: {stats['aggregated_recall']:.3f}")
        print(f"Aggregated F1: {stats['aggregated_f1_score']:.3f}")
        print(f"Aggregated MCC: {stats['aggregated_mcc']:.3f}")


if __name__ == "__main__":
    main()
