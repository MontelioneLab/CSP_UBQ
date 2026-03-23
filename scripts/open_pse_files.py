#!/usr/bin/env python3
"""
Script to open PSE files for a given target in separate PyMOL instances with automatic window positioning.

Usage:
    python scripts/open_pse_files.py <target_name> [options]
    
Examples:
    python scripts/open_pse_files.py 1vj6
    python scripts/open_pse_files.py 1vj6 --screen-size 2560x1440 --window-size 900x700

This script will open the following PSE files (if they exist) in separate PyMOL instances:
- csp_heatmap.pse
- delta_sasa_coloring.pse  
- csp_classification.pse

Windows are automatically positioned to avoid overlap and fit on screen.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Open PSE files for a target in separate PyMOL instances",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "target",
        help="Target name (e.g., 1vj6)"
    )
    parser.add_argument(
        "--screen-size",
        help="Screen resolution in format WIDTHxHEIGHT (default: 1920x1080)",
        default="1920x1080"
    )
    parser.add_argument(
        "--window-size", 
        help="Window size in format WIDTHxHEIGHT (default: 800x600)",
        default="800x600"
    )
    
    args = parser.parse_args()
    target = args.target
    
    # Parse screen and window dimensions
    try:
        screen_width, screen_height = map(int, args.screen_size.split('x'))
        window_width, window_height = map(int, args.window_size.split('x'))
    except ValueError:
        print("Error: Invalid format for screen-size or window-size. Use format WIDTHxHEIGHT (e.g., 1920x1080)")
        sys.exit(1)
    
    # Define the base output directory
    base_output_dir = Path("/Users/tiburon/Desktop/CSP_UBQ/outputs")
    target_dir = base_output_dir / target
    
    # Check if target directory exists
    if not target_dir.exists():
        print(f"Error: Target directory '{target_dir}' does not exist.")
        print(f"Available targets in {base_output_dir}:")
        if base_output_dir.exists():
            for item in sorted(base_output_dir.iterdir()):
                if item.is_dir():
                    print(f"  - {item.name}")
        sys.exit(1)
    
    # Define the PSE files to look for
    pse_files = [
        "csp_heatmap.pse",
        "delta_sasa_coloring.pse", 
        "csp_classification.pse"
    ]
    
    # Check which PSE files exist
    existing_files = []
    missing_files = []
    
    for pse_file in pse_files:
        pse_path = target_dir / pse_file
        if pse_path.exists():
            existing_files.append(pse_path)
        else:
            missing_files.append(pse_file)
    
    # Report status
    if not existing_files:
        print(f"Error: No PSE files found for target '{target}' in {target_dir}")
        print("Expected files:")
        for pse_file in pse_files:
            print(f"  - {pse_file}")
        sys.exit(1)
    
    # Print what we found
    print(f"Opening PSE files for target '{target}':")
    for pse_path in existing_files:
        print(f"  ✓ {pse_path.name}")
    
    if missing_files:
        print("Missing files:")
        for pse_file in missing_files:
            print(f"  ✗ {pse_file}")
    
    # Launch PyMOL instances with automatic window positioning
    print("\nLaunching PyMOL instances...")
    launched_count = 0
    
    # Calculate window positions to avoid overlap
    
    # Calculate positions for windows in a grid layout
    num_windows = len(existing_files)
    if num_windows == 1:
        positions = [(screen_width // 2 - window_width // 2, screen_height // 2 - window_height // 2)]
    elif num_windows == 2:
        # Side by side
        x_offset = (screen_width - 2 * window_width) // 3
        positions = [
            (x_offset, screen_height // 2 - window_height // 2),
            (x_offset * 2 + window_width, screen_height // 2 - window_height // 2)
        ]
    elif num_windows == 3:
        # 2 on top, 1 centered below
        x_offset = (screen_width - 2 * window_width) // 3
        y_offset = (screen_height - 2 * window_height) // 3
        positions = [
            (x_offset, y_offset * 2 + window_height),  # Top left
            (x_offset * 2 + window_width, y_offset * 2 + window_height),  # Top right
            (screen_width // 2 - window_width // 2, y_offset)  # Bottom center
        ]
    else:
        # Fallback: arrange in a grid
        cols = min(3, num_windows)
        rows = (num_windows + cols - 1) // cols
        x_spacing = screen_width // (cols + 1)
        y_spacing = screen_height // (rows + 1)
        positions = []
        for i in range(num_windows):
            col = i % cols
            row = i // cols
            x = x_spacing * (col + 1) - window_width // 2
            y = y_spacing * (row + 1) - window_height // 2
            positions.append((max(0, x), max(0, y)))
    
    for i, pse_path in enumerate(existing_files):
        try:
            x_pos, y_pos = positions[i]
            
            # Launch PyMOL with the PSE file and window geometry
            cmd = [
                "pymol", 
                str(pse_path),
                "-d", f"window,{x_pos},{y_pos},{window_width},{window_height}"
            ]
            subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"  Launched PyMOL with {pse_path.name} at position ({x_pos}, {y_pos})")
            launched_count += 1
        except FileNotFoundError:
            print(f"  Error: PyMOL not found. Please ensure PyMOL is installed and in your PATH.")
            sys.exit(1)
        except Exception as e:
            print(f"  Error launching PyMOL with {pse_path.name}: {e}")
    
    print(f"\nSuccessfully launched {launched_count} PyMOL instance(s).")


if __name__ == "__main__":
    main()
