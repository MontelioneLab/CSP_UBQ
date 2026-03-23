#!/usr/bin/env python3
"""
Script to analyze H_offset and N_offset values from CSP table CSV files.
Creates histograms showing the distribution of offset values across all targets.
"""

import os
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import argparse


def find_csp_table_files(outputs_dir: str) -> List[str]:
    """
    Find all csp_table.csv files in the outputs directory.
    
    Args:
        outputs_dir: Path to the outputs directory
        
    Returns:
        List of paths to csp_table.csv files
    """
    pattern = os.path.join(outputs_dir, "**", "csp_table.csv")
    csv_files = glob.glob(pattern, recursive=True)
    return sorted(csv_files)


def find_grid_offset_files(outputs_dir: str) -> List[str]:
    """
    Find all grid search CSV files with fixed parameters in the outputs directory.

    Args:
        outputs_dir: Path to the outputs directory

    Returns:
        List of paths to grid search CSV files
    """
    pattern = os.path.join(
        outputs_dir,
        "**",
        "offset_grid_H_-0.12_0.12_0.01__N_-1.2_1.2_0.05__C_0.05.csv",
    )
    csv_files = glob.glob(pattern, recursive=True)
    return sorted(csv_files)


def _extract_best_offsets_from_grid_file(
    csv_file: str,
) -> Tuple[Optional[float], Optional[float], Optional[List[float]], Optional[List[float]]]:
    """
    Extract best H/N offsets and grid definitions from a grid CSV file.

    Args:
        csv_file: Path to the grid CSV file

    Returns:
        Tuple of (best_h_offset, best_n_offset, h_values, n_values)
    """
    best_h_offset: Optional[float] = None
    best_n_offset: Optional[float] = None
    h_values: Optional[List[float]] = None
    n_values: Optional[List[float]] = None

    try:
        with open(csv_file, "r", newline="") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith("h_values"):
                    parts = line.split(",")[1:]
                    if parts:
                        h_values = [float(x) for x in parts]
                elif line.startswith("n_values"):
                    parts = line.split(",")[1:]
                    if parts:
                        n_values = [float(x) for x in parts]
                elif line.startswith("best_h_offset"):
                    parts = line.split(",")
                    if len(parts) > 1:
                        try:
                            best_h_offset = float(parts[1])
                        except ValueError:
                            best_h_offset = None
                elif line.startswith("best_n_offset"):
                    parts = line.split(",")
                    if len(parts) > 1:
                        try:
                            best_n_offset = float(parts[1])
                        except ValueError:
                            best_n_offset = None
    except Exception as e:
        print(f"Error reading grid file {csv_file}: {e}")

    return best_h_offset, best_n_offset, h_values, n_values


def collect_best_grid_offsets(
    outputs_dir: str,
) -> Tuple[List[float], List[float], Dict[str, float], Dict[str, float], Optional[List[float]], Optional[List[float]]]:
    """
    Collect best H/N offsets from grid search CSV files.

    Args:
        outputs_dir: Path to the outputs directory

    Returns:
        Tuple containing:
            - List of best H offsets
            - List of best N offsets
            - Mapping of target name -> best H offset
            - Mapping of target name -> best N offset
            - Representative list of H grid values
            - Representative list of N grid values
    """
    csv_files = find_grid_offset_files(outputs_dir)

    all_best_h_offsets: List[float] = []
    all_best_n_offsets: List[float] = []
    best_h_offsets_by_target: Dict[str, float] = {}
    best_n_offsets_by_target: Dict[str, float] = {}
    representative_h_values: Optional[List[float]] = None
    representative_n_values: Optional[List[float]] = None

    print(f"Found {len(csv_files)} grid CSV files")

    for csv_file in csv_files:
        target_name = Path(csv_file).parent.name
        best_h, best_n, h_values, n_values = _extract_best_offsets_from_grid_file(csv_file)

        if best_h is None and best_n is None:
            print(f"Target {target_name}: No best offsets found in grid CSV")
            continue

        print(f"Target {target_name}: best_h_offset={best_h}, best_n_offset={best_n}")

        if best_h is not None:
            all_best_h_offsets.append(best_h)
            best_h_offsets_by_target[target_name] = best_h

        if best_n is not None:
            all_best_n_offsets.append(best_n)
            best_n_offsets_by_target[target_name] = best_n

        # Capture grid definitions from the first file that provides them
        if representative_h_values is None and h_values:
            representative_h_values = h_values
        if representative_n_values is None and n_values:
            representative_n_values = n_values

    return (
        all_best_h_offsets,
        all_best_n_offsets,
        best_h_offsets_by_target,
        best_n_offsets_by_target,
        representative_h_values,
        representative_n_values,
    )


def extract_offset_values(csv_file: str) -> Tuple[float, float, str]:
    """
    Extract one representative H_offset and N_offset value from a CSP table CSV file.
    
    Args:
        csv_file: Path to the CSP table CSV file
        
    Returns:
        Tuple of (H_offset_value, N_offset_value, target_name)
    """
    h_offset = None
    n_offset = None
    
    # Extract target name from path (e.g., outputs/1cf4/csp_table.csv -> 1cf4)
    target_name = Path(csv_file).parent.name
    
    try:
        with open(csv_file, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Extract H_offset and N_offset values
                try:
                    h_offset_str = row.get('H_offset', '').strip()
                    n_offset_str = row.get('N_offset', '').strip()
                    
                    # Take the first valid values found
                    if h_offset is None and h_offset_str and h_offset_str != '':
                        h_offset = float(h_offset_str)
                    if n_offset is None and n_offset_str and n_offset_str != '':
                        n_offset = float(n_offset_str)
                    
                    # If we have both values, we can break early since they're constant per target
                    if h_offset is not None and n_offset is not None:
                        break
                        
                except (ValueError, TypeError):
                    # Skip rows with invalid offset values
                    continue
                    
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return None, None, target_name
    
    return h_offset, n_offset, target_name


def collect_all_offsets(outputs_dir: str) -> Tuple[List[float], List[float], Dict[str, float], Dict[str, float]]:
    """
    Collect one representative H_offset and N_offset value from each CSP table file.
    
    Args:
        outputs_dir: Path to the outputs directory
        
    Returns:
        Tuple of (all_H_offsets, all_N_offsets, H_offsets_by_target, N_offsets_by_target)
    """
    csv_files = find_csp_table_files(outputs_dir)
    
    all_h_offsets = []
    all_n_offsets = []
    h_offsets_by_target = {}
    n_offsets_by_target = {}
    
    print(f"Found {len(csv_files)} CSP table files")
    
    for csv_file in csv_files:
        h_offset, n_offset, target_name = extract_offset_values(csv_file)
        
        if h_offset is not None or n_offset is not None:
            print(f"Target {target_name}: H_offset={h_offset}, N_offset={n_offset}")
            
            if h_offset is not None:
                all_h_offsets.append(h_offset)
                h_offsets_by_target[target_name] = h_offset
            
            if n_offset is not None:
                all_n_offsets.append(n_offset)
                n_offsets_by_target[target_name] = n_offset
        else:
            print(f"Target {target_name}: No valid offset values found")
    
    return all_h_offsets, all_n_offsets, h_offsets_by_target, n_offsets_by_target


def _count_axis_ticks(max_count: float, min_ticks: int = 4, max_ticks: int = 5) -> List[float]:
    """
    Compute tick values for a count/frequency axis so there are between
    min_ticks and max_ticks (inclusive), using nice increments.

    Args:
        max_count: Maximum count value to span (e.g. max of histogram)
        min_ticks: Minimum number of ticks (default 4)
        max_ticks: Maximum number of ticks (default 5)

    Returns:
        List of tick values including 0, with 4–5 elements when possible.
    """
    v_max = max(1.0, float(np.ceil(max_count)))
    magnitude = 10.0 ** np.floor(np.log10(v_max)) if v_max > 0 else 1.0
    # Nice step candidates: 0.5, 1, 2, 5, 10, 20, ... (to support 4–5 ticks for small counts)
    nice_candidates = [c * magnitude for c in [0.5, 1, 2, 5, 10]]
    # For 5 ticks: step >= v_max/4. For 4 ticks: step >= v_max/3.
    step_5 = v_max / 4.0
    step_4 = v_max / 3.0
    step = None
    for c in nice_candidates:
        n = 1 + int(np.ceil((v_max + 1e-9) / c))
        if min_ticks <= n <= max_ticks and c >= step_4:
            step = c
            break
    if step is None:
        # Prefer step that gives 5 ticks, else 4; otherwise closest in 4–5
        for c in nice_candidates:
            if c >= step_5:
                step = c
                break
        if step is None:
            step = max(step_4, magnitude)
    # Build ticks: 0, step, 2*step, ... up to >= v_max
    ticks = [0.0]
    t = step
    while t <= v_max + step * 0.01 and len(ticks) <= max_ticks:
        ticks.append(t)
        t += step
    # If we have 6, drop the last if it's only slightly above v_max
    if len(ticks) > max_ticks and ticks[-1] > v_max + 1e-9:
        ticks = [x for x in ticks if x <= v_max + step * 0.01]
    return ticks


def _calculate_bin_edges(values: List[float]) -> np.ndarray:
    """
    Calculate bin edges for a sequence of grid center values.

    Args:
        values: Ordered list of grid centers

    Returns:
        Numpy array of bin edges
    """
    if not values:
        raise ValueError("Cannot calculate bin edges from an empty values list")

    if len(values) == 1:
        step = 0.01
    else:
        steps = np.diff(values)
        step = float(np.median(steps))

    start = values[0] - step / 2.0
    end = values[-1] + step / 2.0

    # Use linspace to avoid floating point drift
    edges = np.linspace(start, end, len(values) + 1)
    return edges


def create_grid_heatmap(
    best_h_offsets: List[float],
    best_n_offsets: List[float],
    output_dir: str,
    h_grid_values: Optional[List[float]] = None,
    n_grid_values: Optional[List[float]] = None,
) -> Optional[str]:
    """
    Create a heatmap showing counts of best offsets across the grid.

    Args:
        best_h_offsets: List of best H offsets
        best_n_offsets: List of best N offsets
        output_dir: Directory to save the heatmap
        h_grid_values: Optional list of H grid center values
        n_grid_values: Optional list of N grid center values

    Returns:
        Path to the saved heatmap image, or None if there was nothing to plot
    """
    if not best_h_offsets or not best_n_offsets:
        return None

    os.makedirs(output_dir, exist_ok=True)

    # Fall back to min/max derived ranges if grid definitions are missing
    if not h_grid_values:
        h_min = min(best_h_offsets)
        h_max = max(best_h_offsets)
        h_grid_values = list(np.linspace(h_min, h_max, max(3, len(set(best_h_offsets)))))

    if not n_grid_values:
        n_min = min(best_n_offsets)
        n_max = max(best_n_offsets)
        n_grid_values = list(np.linspace(n_min, n_max, max(3, len(set(best_n_offsets)))))

    h_edges = _calculate_bin_edges(h_grid_values)
    n_edges = _calculate_bin_edges(n_grid_values)

    # Histogram the data
    heatmap_counts, h_edges_hist, n_edges_hist = np.histogram2d(
        best_h_offsets, best_n_offsets, bins=[h_edges, n_edges]
    )

    fig = plt.figure(figsize=(12, 8))  # Slightly wider figsize to account for colorbar column
    gs = GridSpec(
        2,
        3,
        width_ratios=[1.2, 4, 0.3],  # Added narrow column for colorbar
        height_ratios=[4, 1.2],
        hspace=0.0,
        wspace=0.0,
    )

    ax_heatmap = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[0, 0], sharey=ax_heatmap)
    ax_bottom = fig.add_subplot(gs[1, 1], sharex=ax_heatmap)
    ax_empty = fig.add_subplot(gs[1, 0])
    ax_empty.axis("off")
    ax_cbar = fig.add_subplot(gs[0, 2])  # Dedicated axes for colorbar

    title_fontsize = 18
    label_fontsize = 14
    tick_fontsize = 12

    cmap = plt.get_cmap("viridis")
    mesh = ax_heatmap.imshow(
        heatmap_counts,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        extent=[n_edges_hist[0], n_edges_hist[-1], h_edges_hist[0], h_edges_hist[-1]],
    )
    ax_heatmap.set_xlabel("")
    # ax_heatmap.set_ylabel("H offset (ppm)")
    ax_heatmap.set_title(
        "Best Offsets from Grid Search Across Targets",
        pad=20,
        fontsize=title_fontsize,
    )
    ax_heatmap.grid(False)
    ax_heatmap.tick_params(axis="x", labelbottom=False)
    # Use fixed ranges for N offset (-1.225 to 1.225) and H offset (-0.125 to 0.125)
    ax_heatmap.set_xlim(-1.225, 1.225)
    ax_heatmap.set_ylim(-0.125, 0.125)
    # ax_heatmap.yaxis.set_label_position("right")
    # ax_heatmap.yaxis.tick_right()
    ax_heatmap.tick_params(
        axis="y", which="both", labelleft=False, labelright=False, labelsize=tick_fontsize
    )
    ax_heatmap.tick_params(axis="x", labelsize=tick_fontsize)

    cbar = fig.colorbar(mesh, cax=ax_cbar)  # Use dedicated cax instead of ax=ax_heatmap
    cbar.set_label("Target count", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)
    max_heatmap_count = int(np.max(heatmap_counts)) if heatmap_counts.size else 0
    cbar_ticks = list(range(0, max_heatmap_count + 1))
    if cbar_ticks:
        cbar.set_ticks(cbar_ticks)

    # Fixed tick ranges for histograms: N offset -1.225 to 1.225, H offset -0.125 to 0.125,
    # both with equal steps and 0.00 included.
    tick_positions_n = np.linspace(-1.225, 1.225, 11)  # 11 ticks, equal steps, includes 0.00
    tick_positions_h = np.linspace(-0.125, 0.125, 11)  # 11 ticks, equal steps, includes 0.00

    ax_heatmap.set_xticks(tick_positions_n, minor=False)
    ax_heatmap.set_yticks(tick_positions_h, minor=False)

    # Bottom histogram for N offsets
    # Bottom histogram for N offsets. Use the FULL set of bin edges
    # (n_edges) so each bar lines up one-to-one with a heatmap column.
    n_hist_counts, _, _ = ax_bottom.hist(
        best_n_offsets,
        bins=n_edges,
        color=cmap(0.65),
        edgecolor="black",
        alpha=0.7,
    )
    ax_bottom.set_xlabel("N offset (ppm)", fontsize=label_fontsize)
    ax_bottom.set_ylabel("Count", fontsize=label_fontsize)
    ax_bottom.grid(False)
    ax_bottom.invert_yaxis()
    ax_bottom.set_xticks(tick_positions_n)
    ax_bottom.set_xlim(ax_heatmap.get_xlim())
    ax_bottom.tick_params(
        axis="x",
        top=False,
        bottom=True,
        labelbottom=True,
        labelsize=tick_fontsize,
    )
    ax_bottom.tick_params(axis="y", labelsize=tick_fontsize)

    # Left histogram for H offsets (horizontal)
    # Left histogram for H offsets (horizontal), likewise using the full
    # set of bin edges (h_edges) so each bar corresponds to a single
    # heatmap row.
    h_hist_counts, _, _ = ax_left.hist(
        best_h_offsets,
        bins=h_edges,
        orientation="horizontal",
        color=cmap(0.35),
        edgecolor="black",
        alpha=0.7,
    )
    ax_left.set_ylabel("H offset (ppm)", fontsize=label_fontsize)
    ax_left.set_xlabel("Count", fontsize=label_fontsize)
    ax_left.grid(False)
    ax_left.invert_xaxis()
    ax_left.tick_params(
        axis="y",
        labelleft=True,
        labelright=False,
        labelsize=tick_fontsize,
    )
    ax_left.tick_params(axis="x", labelsize=tick_fontsize)
    # Make sure the left histogram uses the same tick positions
    # as the heatmap so the H offset bins line up.
    ax_left.set_yticks(tick_positions_h, minor=False)
    ax_left.set_ylim(ax_heatmap.get_ylim())

    # Count axis: 4–5 ticks with dynamic increments; exclude 0 from margin histograms
    if n_hist_counts.size:
        ticks_n = _count_axis_ticks(float(np.max(n_hist_counts)))
        if ticks_n:
            ticks_n = [t for t in ticks_n if t != 0]
            if ticks_n:
                ax_bottom.set_yticks(ticks_n)
    if h_hist_counts.size:
        ticks_h = _count_axis_ticks(float(np.max(h_hist_counts)))
        if ticks_h:
            ticks_h = [t for t in ticks_h if t != 0]
            if ticks_h:
                ax_left.set_xticks(ticks_h)

    # Sync tick labels visibility and tick placement
    plt.setp(ax_bottom.get_xticklabels(), rotation=45, ha="right")
    plt.setp(ax_bottom.get_yticklabels(), rotation=0)
    fig.tight_layout()
    output_path = os.path.join(output_dir, "offset_grid_heatmap.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return output_path


def create_histograms(all_h_offsets: List[float], all_n_offsets: List[float], 
                     h_offsets_by_target: Dict[str, float], 
                     n_offsets_by_target: Dict[str, float], 
                     output_dir: str) -> None:
    """
    Create histograms for H_offset and N_offset distributions.
    
    Args:
        all_h_offsets: All H_offset values across all targets
        all_n_offsets: All N_offset values across all targets
        h_offsets_by_target: H_offset values grouped by target
        n_offsets_by_target: N_offset values grouped by target
        output_dir: Directory to save the plots
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up the plotting style
    plt.style.use('default')
    fig_width = 12
    fig_height = 8
    
    # 1. Overall histograms for all targets combined
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_height))
    
    # H_offset histogram
    if all_h_offsets:
        ax1.hist(all_h_offsets, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax1.set_xlabel('H_offset (ppm)')
        ax1.set_ylabel('Frequency')
        ax1.set_title(f'H_offset Distribution (All Targets)\nn={len(all_h_offsets)} values')
        ax1.grid(True, alpha=0.3)
        
        # Add statistics
        mean_h = np.mean(all_h_offsets)
        std_h = np.std(all_h_offsets)
        ax1.axvline(mean_h, color='red', linestyle='--', label=f'Mean: {mean_h:.3f}')
        ax1.axvline(mean_h + std_h, color='orange', linestyle='--', alpha=0.7, label=f'+1σ: {mean_h + std_h:.3f}')
        ax1.axvline(mean_h - std_h, color='orange', linestyle='--', alpha=0.7, label=f'-1σ: {mean_h - std_h:.3f}')
        ax1.legend()
    else:
        ax1.text(0.5, 0.5, 'No H_offset data', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('H_offset Distribution (No Data)')
    
    # N_offset histogram
    if all_n_offsets:
        ax2.hist(all_n_offsets, bins=50, alpha=0.7, color='green', edgecolor='black')
        ax2.set_xlabel('N_offset (ppm)')
        ax2.set_ylabel('Frequency')
        ax2.set_title(f'N_offset Distribution (All Targets)\nn={len(all_n_offsets)} values')
        ax2.grid(True, alpha=0.3)
        
        # Add statistics
        mean_n = np.mean(all_n_offsets)
        std_n = np.std(all_n_offsets)
        ax2.axvline(mean_n, color='red', linestyle='--', label=f'Mean: {mean_n:.3f}')
        ax2.axvline(mean_n + std_n, color='orange', linestyle='--', alpha=0.7, label=f'+1σ: {mean_n + std_n:.3f}')
        ax2.axvline(mean_n - std_n, color='orange', linestyle='--', alpha=0.7, label=f'-1σ: {mean_n - std_n:.3f}')
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, 'No N_offset data', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('N_offset Distribution (No Data)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'offset_distributions_all_targets.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Scatter plots by target (since we have one value per target)
    if h_offsets_by_target or n_offsets_by_target:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width, fig_height))
        
        # H_offset scatter plot by target
        if h_offsets_by_target:
            targets_with_h = list(h_offsets_by_target.keys())
            h_values = list(h_offsets_by_target.values())
            
            ax1.scatter(range(len(targets_with_h)), h_values, alpha=0.7, s=50, color='blue')
            ax1.set_ylabel('H_offset (ppm)')
            ax1.set_title('H_offset Values by Target')
            ax1.set_xticks(range(len(targets_with_h)))
            ax1.set_xticklabels(targets_with_h, rotation=45)
            ax1.grid(True, alpha=0.3)
        else:
            ax1.text(0.5, 0.5, 'No H_offset data', ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title('H_offset Values by Target (No Data)')
        
        # N_offset scatter plot by target
        if n_offsets_by_target:
            targets_with_n = list(n_offsets_by_target.keys())
            n_values = list(n_offsets_by_target.values())
            
            ax2.scatter(range(len(targets_with_n)), n_values, alpha=0.7, s=50, color='green')
            ax2.set_ylabel('N_offset (ppm)')
            ax2.set_title('N_offset Values by Target')
            ax2.set_xticks(range(len(targets_with_n)))
            ax2.set_xticklabels(targets_with_n, rotation=45)
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'No N_offset data', ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('N_offset Values by Target (No Data)')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'offset_values_by_target.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Scatter plot of H_offset vs N_offset
    if all_h_offsets and all_n_offsets:
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
        
        # Create scatter plot
        ax.scatter(all_h_offsets, all_n_offsets, alpha=0.6, s=20)
        ax.set_xlabel('H_offset (ppm)')
        ax.set_ylabel('N_offset (ppm)')
        ax.set_title(f'H_offset vs N_offset Correlation\nn={len(all_h_offsets)} points')
        ax.grid(True, alpha=0.3)
        
        # Add correlation coefficient
        if len(all_h_offsets) > 1 and len(all_n_offsets) > 1:
            correlation = np.corrcoef(all_h_offsets, all_n_offsets)[0, 1]
            ax.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                   transform=ax.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'offset_correlation.png'), dpi=300, bbox_inches='tight')
        plt.close()


def save_summary_statistics(all_h_offsets: List[float], all_n_offsets: List[float], 
                          h_offsets_by_target: Dict[str, float], 
                          n_offsets_by_target: Dict[str, float], 
                          output_dir: str) -> None:
    """
    Save summary statistics to a CSV file.
    
    Args:
        all_h_offsets: All H_offset values across all targets
        all_n_offsets: All N_offset values across all targets
        h_offsets_by_target: H_offset values grouped by target
        n_offsets_by_target: N_offset values grouped by target
        output_dir: Directory to save the statistics file
    """
    
    stats_data = []
    
    # Overall statistics
    if all_h_offsets:
        stats_data.append({
            'Target': 'ALL_TARGETS',
            'Offset_Type': 'H_offset',
            'Count': len(all_h_offsets),
            'Mean': np.mean(all_h_offsets),
            'Std': np.std(all_h_offsets),
            'Min': np.min(all_h_offsets),
            'Max': np.max(all_h_offsets),
            'Median': np.median(all_h_offsets),
            'Q25': np.percentile(all_h_offsets, 25),
            'Q75': np.percentile(all_h_offsets, 75)
        })
    
    if all_n_offsets:
        stats_data.append({
            'Target': 'ALL_TARGETS',
            'Offset_Type': 'N_offset',
            'Count': len(all_n_offsets),
            'Mean': np.mean(all_n_offsets),
            'Std': np.std(all_n_offsets),
            'Min': np.min(all_n_offsets),
            'Max': np.max(all_n_offsets),
            'Median': np.median(all_n_offsets),
            'Q25': np.percentile(all_n_offsets, 25),
            'Q75': np.percentile(all_n_offsets, 75)
        })
    
    # Per-target statistics
    for target in set(list(h_offsets_by_target.keys()) + list(n_offsets_by_target.keys())):
        h_value = h_offsets_by_target.get(target)
        n_value = n_offsets_by_target.get(target)
        
        if h_value is not None:
            stats_data.append({
                'Target': target,
                'Offset_Type': 'H_offset',
                'Count': 1,
                'Mean': h_value,
                'Std': 0.0,
                'Min': h_value,
                'Max': h_value,
                'Median': h_value,
                'Q25': h_value,
                'Q75': h_value
            })
        
        if n_value is not None:
            stats_data.append({
                'Target': target,
                'Offset_Type': 'N_offset',
                'Count': 1,
                'Mean': n_value,
                'Std': 0.0,
                'Min': n_value,
                'Max': n_value,
                'Median': n_value,
                'Q25': n_value,
                'Q75': n_value
            })
    
    # Save to CSV
    stats_file = os.path.join(output_dir, 'offset_statistics.csv')
    with open(stats_file, 'w', newline='') as f:
        if stats_data:
            writer = csv.DictWriter(f, fieldnames=stats_data[0].keys())
            writer.writeheader()
            writer.writerows(stats_data)
    print(f"Statistics saved to: {stats_file}")


def main():
    """Main function to analyze offset values and create histograms."""
    
    parser = argparse.ArgumentParser(description='Analyze H_offset and N_offset values from CSP table CSV files')
    parser.add_argument('--outputs-dir', default='outputs', 
                       help='Path to outputs directory containing CSP table files (default: outputs)')
    parser.add_argument('--output-dir', default='offset_analysis', 
                       help='Directory to save analysis results (default: offset_analysis)')
    
    args = parser.parse_args()
    
    print("Analyzing H_offset and N_offset values from CSP table files...")
    print(f"Outputs directory: {args.outputs_dir}")
    print(f"Results will be saved to: {args.output_dir}")
    
    # Check if outputs directory exists
    if not os.path.exists(args.outputs_dir):
        print(f"Error: Outputs directory '{args.outputs_dir}' does not exist")
        return
    
    # Collect best offsets from grid search CSV files
    heatmap_path: Optional[str] = None

    print("\nCollecting best offsets from grid search files...")
    (
        best_grid_h_offsets,
        best_grid_n_offsets,
        best_h_offsets_by_target,
        best_n_offsets_by_target,
        h_grid_values,
        n_grid_values,
    ) = collect_best_grid_offsets(args.outputs_dir)

    if best_grid_h_offsets or best_grid_n_offsets:
        print(f"Collected best offsets for {len(set(best_h_offsets_by_target.keys()).union(best_n_offsets_by_target.keys()))} targets")
        heatmap_path = create_grid_heatmap(
            best_grid_h_offsets,
            best_grid_n_offsets,
            args.output_dir,
            h_grid_values,
            n_grid_values,
        )
        if heatmap_path:
            print(f"Saved grid heatmap to: {heatmap_path}")
        else:
            print("No heatmap created (insufficient grid offset data).")
    else:
        print("No grid CSV files with best offsets were found.")

    # Collect all offset values from CSP tables
    all_h_offsets, all_n_offsets, h_offsets_by_target, n_offsets_by_target = collect_all_offsets(args.outputs_dir)
    
    if not all_h_offsets and not all_n_offsets:
        print("No offset values found in any CSP table files")
        return
    
    print(f"\nSummary:")
    print(f"Total grid-best H_offset values: {len(best_grid_h_offsets)}")
    print(f"Total grid-best N_offset values: {len(best_grid_n_offsets)}")
    print(f"Total H_offset values: {len(all_h_offsets)}")
    print(f"Total N_offset values: {len(all_n_offsets)}")
    print(f"Number of targets with data: {len(set(list(h_offsets_by_target.keys()) + list(n_offsets_by_target.keys())))}")
    
    # Create histograms
    print("\nCreating histograms...")
    create_histograms(all_h_offsets, all_n_offsets, h_offsets_by_target, n_offsets_by_target, args.output_dir)
    
    # Save summary statistics
    print("Saving summary statistics...")
    save_summary_statistics(all_h_offsets, all_n_offsets, h_offsets_by_target, n_offsets_by_target, args.output_dir)
    
    print(f"\nAnalysis complete! Results saved to: {args.output_dir}")
    print("Generated files:")
    if heatmap_path:
        print("- offset_grid_heatmap.png")
    print("- offset_distributions_all_targets.png")
    print("- offset_values_by_target.png") 
    print("- offset_correlation.png")
    print("- offset_statistics.csv")


if __name__ == "__main__":
    main()
