"""
CSP Z-Score vs Distance Heatmap Visualization

This script creates a 2D heatmap showing the distribution of data points
across csp_z (Y-axis, binned with range 1) and min_ca_distance_distance 
(X-axis, auto-determined bins) from all master_alignment.csv files.

Usage:
    python scripts/visualize_csp_z_distance.py [--outputs outputs/] [--output outputs/csp_z_distance_heatmap.png]
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional

try:
    from .config import paths
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import paths


def find_master_alignment_files(outputs_dir: str) -> List[str]:
    """
    Find all master_alignment.csv files in the outputs directory.
    
    Args:
        outputs_dir: Path to the outputs directory
        
    Returns:
        List of paths to master_alignment.csv files
    """
    pattern = os.path.join(outputs_dir, "**", "master_alignment.csv")
    csv_files = glob.glob(pattern, recursive=True)
    return sorted(csv_files)


def collect_data_from_files(csv_files: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read all master_alignment.csv files and extract csp_z and min_ca_distance_distance columns.
    
    Args:
        csv_files: List of paths to master_alignment.csv files
        
    Returns:
        Tuple of (csp_z_values, distance_values) as numpy arrays
    """
    csp_z_values = []
    distance_values = []
    
    for csv_file in csv_files:
        try:
            with open(csv_file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Extract csp_z
                    csp_z_str = row.get('csp_z', '').strip()
                    # Extract min_ca_distance_distance
                    distance_str = row.get('min_ca_distance_distance', '').strip()
                    
                    # Convert to float, skip if invalid
                    try:
                        if csp_z_str and distance_str:
                            csp_z = float(csp_z_str)
                            distance = float(distance_str)
                            # Only include valid finite values
                            if np.isfinite(csp_z) and np.isfinite(distance):
                                csp_z_values.append(csp_z)
                                distance_values.append(distance)
                    except (ValueError, TypeError):
                        continue
        except Exception as e:
            print(f"Warning: Error reading {csv_file}: {e}")
            continue
    
    return np.array(csp_z_values), np.array(distance_values)


def create_bins(csp_z_values: np.ndarray, distance_values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create bin edges for csp_z (range 1) and min_ca_distance_distance (auto-determined).
    Values above the upper bounds are capped: csp_z > 8 → 8, distance > 20 → 20.
    
    Args:
        csp_z_values: Array of csp_z values (already clipped)
        distance_values: Array of min_ca_distance_distance values (already clipped)
        
    Returns:
        Tuple of (csp_z_bins, distance_bins) as numpy arrays of bin edges
    """
    # Y-axis: csp_z with range of 1, capped at 8
    CSP_Z_MAX = 8.0
    csp_z_min = np.floor(np.min(csp_z_values))
    # Always create bins up to the cap (8.0)
    csp_z_max = CSP_Z_MAX
    # Create integer bins with range 1, ensuring we include the max bin
    csp_z_bins = np.arange(csp_z_min, csp_z_max + 1, 1.0)
    # Ensure the last bin edge is exactly at the cap
    if csp_z_bins[-1] != CSP_Z_MAX:
        csp_z_bins = np.append(csp_z_bins, CSP_Z_MAX)
    
    # X-axis: min_ca_distance_distance with auto-determined bins, capped at 20
    DISTANCE_MAX = 20.0
    distance_min = np.min(distance_values)
    # Always create bins up to the cap (20.0)
    distance_max = DISTANCE_MAX
    distance_range = distance_max - distance_min
    
    # Determine sensible bin size based on data range
    # Use 1-2 Å bins depending on data spread
    if distance_range <= 10:
        bin_size = 0.5  # 0.5 Å bins for small ranges
    elif distance_range <= 25:
        bin_size = 1.0  # 1 Å bins for medium ranges
    elif distance_range <= 50:
        bin_size = 2.0  # 2 Å bins for larger ranges
    else:
        bin_size = 5.0  # 5 Å bins for very large ranges
    
    # Create bins up to the cap
    num_bins = int(np.ceil(distance_range / bin_size))
    distance_bins = np.linspace(distance_min, distance_max, num_bins + 1)
    # Ensure the last bin edge is exactly at the cap
    if abs(distance_bins[-1] - DISTANCE_MAX) > 1e-10:
        distance_bins[-1] = DISTANCE_MAX
    
    return csp_z_bins, distance_bins


def create_heatmap(csp_z_values: np.ndarray, distance_values: np.ndarray, 
                   output_path: str) -> None:
    """
    Create and save a 2D heatmap showing the distribution of data points.
    Values above upper bounds are capped: csp_z > 8 → 8, distance > 20 → 20.
    
    Args:
        csp_z_values: Array of csp_z values
        distance_values: Array of min_ca_distance_distance values
        output_path: Path to save the heatmap PNG
    """
    if len(csp_z_values) == 0 or len(distance_values) == 0:
        print("Error: No valid data points found")
        return
    
    # Cap values at upper bounds
    CSP_Z_MAX = 8.0
    DISTANCE_MAX = 20.0
    csp_z_clipped = np.clip(csp_z_values, None, CSP_Z_MAX)
    distance_clipped = np.clip(distance_values, None, DISTANCE_MAX)
    
    # Count how many values were capped
    csp_z_capped_count = np.sum(csp_z_values > CSP_Z_MAX)
    distance_capped_count = np.sum(distance_values > DISTANCE_MAX)
    
    if csp_z_capped_count > 0:
        print(f"Note: {csp_z_capped_count} csp_z values > {CSP_Z_MAX} were capped to {CSP_Z_MAX}")
    if distance_capped_count > 0:
        print(f"Note: {distance_capped_count} distance values > {DISTANCE_MAX} were capped to {DISTANCE_MAX}")
    
    # Create bins using clipped values
    csp_z_bins, distance_bins = create_bins(csp_z_clipped, distance_clipped)
    
    # Create 2D histogram using clipped values
    heatmap_counts, csp_z_edges, distance_edges = np.histogram2d(
        csp_z_clipped, distance_clipped, bins=[csp_z_bins, distance_bins]
    )
    
    # Transpose so that csp_z is on Y-axis (rows) and distance is on X-axis (columns)
    heatmap_counts = heatmap_counts.T
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create heatmap using imshow
    im = ax.imshow(
        heatmap_counts,
        origin='lower',
        aspect='auto',
        cmap='viridis',
        interpolation='nearest',
        extent=[distance_edges[0], distance_edges[-1], csp_z_edges[0], csp_z_edges[-1]]
    )
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Count', rotation=270, labelpad=20)
    
    # Set axis labels
    ax.set_xlabel('min_ca_distance_distance (Å)', fontsize=12)
    ax.set_ylabel('csp_z', fontsize=12)
    ax.set_title('CSP Z-Score vs Minimum CA Distance Heatmap', fontsize=14, fontweight='bold')
    
    # Format tick labels
    # For distance (X-axis), show every few bins to avoid crowding
    num_x_ticks = min(10, len(distance_edges) - 1)
    x_tick_indices = np.linspace(0, len(distance_edges) - 1, num_x_ticks, dtype=int)
    ax.set_xticks([distance_edges[i] for i in x_tick_indices])
    ax.set_xticklabels([f'{distance_edges[i]:.1f}' for i in x_tick_indices], rotation=45, ha='right')
    
    # For csp_z (Y-axis), show every few bins
    num_y_ticks = min(15, len(csp_z_edges) - 1)
    y_tick_indices = np.linspace(0, len(csp_z_edges) - 1, num_y_ticks, dtype=int)
    ax.set_yticks([csp_z_edges[i] for i in y_tick_indices])
    ax.set_yticklabels([f'{csp_z_edges[i]:.1f}' for i in y_tick_indices])
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved heatmap to {output_path}")
    print(f"Total data points: {len(csp_z_values)}")
    print(f"csp_z range (original): {np.min(csp_z_values):.2f} to {np.max(csp_z_values):.2f}")
    print(f"csp_z range (displayed): {np.min(csp_z_clipped):.2f} to {np.max(csp_z_clipped):.2f} (capped at {CSP_Z_MAX})")
    print(f"Distance range (original): {np.min(distance_values):.2f} to {np.max(distance_values):.2f} Å")
    print(f"Distance range (displayed): {np.min(distance_clipped):.2f} to {np.max(distance_clipped):.2f} Å (capped at {DISTANCE_MAX})")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create 2D heatmap of CSP Z-score vs minimum CA distance"
    )
    parser.add_argument(
        "--outputs",
        default=paths.outputs_dir if hasattr(paths, 'outputs_dir') else "outputs",
        help="Path to outputs directory containing system subdirectories (default: outputs/)"
    )
    parser.add_argument(
        "--output",
        default=os.path.join(paths.outputs_dir if hasattr(paths, 'outputs_dir') else "outputs", "csp_z_distance_heatmap.png"),
        help="Output path for the heatmap PNG (default: outputs/csp_z_distance_heatmap.png)"
    )
    
    args = parser.parse_args()
    
    # Find all master_alignment.csv files
    print(f"Searching for master_alignment.csv files in {args.outputs}...")
    csv_files = find_master_alignment_files(args.outputs)
    
    if not csv_files:
        print(f"Error: No master_alignment.csv files found in {args.outputs}")
        sys.exit(1)
    
    print(f"Found {len(csv_files)} master_alignment.csv files")
    
    # Collect data from all files
    print("Reading data from CSV files...")
    csp_z_values, distance_values = collect_data_from_files(csv_files)
    
    if len(csp_z_values) == 0:
        print("Error: No valid data points found in any CSV files")
        sys.exit(1)
    
    # Create and save heatmap
    print("Creating heatmap...")
    create_heatmap(csp_z_values, distance_values, args.output)
    
    print("Done!")


if __name__ == "__main__":
    main()

