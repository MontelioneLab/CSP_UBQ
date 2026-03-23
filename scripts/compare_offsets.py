#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Optional, Tuple
import math

try:
    from scripts.config import compute
except ImportError:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts.config import compute


def parse_grid_csv_best_offsets(csv_path: str) -> Optional[Tuple[float, float]]:
    best_n = None
    best_h = None
    try:
        with open(csv_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split(',')
                if parts[0] == 'best_n_offset' and len(parts) > 1:
                    best_n = float(parts[1])
                elif parts[0] == 'best_h_offset' and len(parts) > 1:
                    best_h = float(parts[1])
        if best_n is None or best_h is None:
            return None
        return (best_n, best_h)
    except Exception:
        return None


def pick_best_grid_csv(file_names: List[str]) -> Optional[str]:
    """Pick the grid CSV with the finest resolution based on step sizes in the filename."""
    best = None
    best_h_step = None
    best_n_step = None
    for name in file_names:
        # Expect pattern: offset_grid_H_<hmin>_<hmax>_<hstep>__N_<nmin>_<nmax>_<nstep>__C_<cut>.csv
        base = os.path.basename(name)
        try:
            parts = base.split('__')
            h_part = parts[0]  # H_<hmin>_<hmax>_<hstep>
            n_part = parts[1]  # N_<nmin>_<nmax>_<nstep>
            h_tokens = h_part.split('_')
            n_tokens = n_part.split('_')
            h_step = float(h_tokens[-1])
            n_step = float(n_tokens[-1])
            if best is None or (h_step < (best_h_step or 1e9)) or (h_step == best_h_step and n_step < (best_n_step or 1e9)):
                best = name
                best_h_step = h_step
                best_n_step = n_step
        except Exception:
            # Fallback: choose lexicographically smallest if parsing fails
            if best is None or base < os.path.basename(best):
                best = name
    return best


def compute_mean_offsets_from_csp_table(csp_table_csv: str, cutoff: float) -> Optional[Tuple[float, float]]:
    """Compute mean-based offsets using original holo shifts, referencing subset CSP < cutoff.
    Returns (n_offset, h_offset) or None if insufficient data.
    """
    n_diffs: List[float] = []
    h_diffs: List[float] = []
    try:
        with open(csp_table_csv, 'r', newline='') as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                try:
                    h_a = float(row['H_apo']) if row.get('H_apo') else None
                    n_a = float(row['N_apo']) if row.get('N_apo') else None
                    h_h_orig = float(row['H_holo_original']) if row.get('H_holo_original') else None
                    n_h_orig = float(row['N_holo_original']) if row.get('N_holo_original') else None
                except Exception:
                    continue
                if h_a is None or n_a is None or h_h_orig is None or n_h_orig is None:
                    continue
                dH = h_h_orig - h_a
                dN = n_h_orig - n_a
                wn = compute.csp_delta_n_scale
                csp_val = math.sqrt(0.5 * (dH * dH + (dN * wn) ** 2))
                if csp_val < cutoff:
                    # Mean method uses apo - holo differences
                    n_diffs.append(n_a - n_h_orig)
                    h_diffs.append(h_a - h_h_orig)
        if not n_diffs and not h_diffs:
            return None
        n_offset = sum(n_diffs) / len(n_diffs) if n_diffs else 0.0
        h_offset = sum(h_diffs) / len(h_diffs) if h_diffs else 0.0
        return (n_offset, h_offset)
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare N/H referencing offsets between grid and mean methods")
    parser.add_argument('--outputs-dir', default='outputs', help='Path to outputs directory')
    parser.add_argument('--cutoff', type=float, default=0.05, help='CSP cutoff for mean-based referencing subset (default: 0.05)')
    parser.add_argument('--summary-csv', default='outputs/offset_method_comparison.csv', help='Path to write summary CSV')
    parser.add_argument('--histogram-png', default='outputs/offset_method_histograms.png', help='Path to write histogram image')
    args = parser.parse_args()

    outputs_dir = args.outputs_dir
    targets = [d for d in os.listdir(outputs_dir) if os.path.isdir(os.path.join(outputs_dir, d))]

    rows: List[Dict[str, object]] = []
    diffs_h: List[float] = []
    diffs_n: List[float] = []
    grid_h_values: List[float] = []
    grid_n_values: List[float] = []
    mean_h_values: List[float] = []
    mean_n_values: List[float] = []

    for tgt in sorted(targets):
        tgt_dir = os.path.join(outputs_dir, tgt)
        # grid: pick best resolution CSV
        grid_csv_candidates = [os.path.join(tgt_dir, f) for f in os.listdir(tgt_dir) if f.startswith('offset_grid_') and f.endswith('.csv')]
        grid_n = None
        grid_h = None
        if grid_csv_candidates:
            chosen = pick_best_grid_csv(grid_csv_candidates)
            if chosen:
                best = parse_grid_csv_best_offsets(chosen)
                if best is not None:
                    grid_n, grid_h = best

        # mean: recompute from csp_table.csv using original holo shifts
        csp_table_csv = os.path.join(tgt_dir, 'csp_table.csv')
        mean_n = None
        mean_h = None
        if os.path.exists(csp_table_csv):
            mean = compute_mean_offsets_from_csp_table(csp_table_csv, cutoff=args.cutoff)
            if mean is not None:
                mean_n, mean_h = mean

        # Only include if we have both
        if grid_n is not None and grid_h is not None and mean_n is not None and mean_h is not None:
            delta_h = mean_h - grid_h
            delta_n = mean_n - grid_n
            rows.append({
                'target': tgt,
                'grid_n_offset': grid_n,
                'grid_h_offset': grid_h,
                'mean_n_offset': mean_n,
                'mean_h_offset': mean_h,
                'delta_n': delta_n,
                'delta_h': delta_h,
            })
            diffs_h.append(abs(delta_h))
            diffs_n.append(abs(delta_n))
            grid_h_values.append(grid_h)
            grid_n_values.append(grid_n)
            mean_h_values.append(mean_h)
            mean_n_values.append(mean_n)

    # Write summary CSV
    os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)
    with open(args.summary_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['target','grid_n_offset','grid_h_offset','mean_n_offset','mean_h_offset','delta_n','delta_h'])
        for r in rows:
            w.writerow([
                r['target'],
                f"{r['grid_n_offset']:.6f}", f"{r['grid_h_offset']:.6f}",
                f"{r['mean_n_offset']:.6f}", f"{r['mean_h_offset']:.6f}",
                f"{r['delta_n']:.6f}", f"{r['delta_h']:.6f}",
            ])
        # Append averages as a final row
        avg_abs_n = sum(diffs_n)/len(diffs_n) if diffs_n else 0.0
        avg_abs_h = sum(diffs_h)/len(diffs_h) if diffs_h else 0.0
        w.writerow(['AVERAGE_ABS_DIFF','','','','', f"{avg_abs_n:.6f}", f"{avg_abs_h:.6f}"])

    print(f"Compared {len(rows)} targets. Summary: {args.summary_csv}")
    if diffs_n and diffs_h:
        print(f"Average abs diff N: {sum(diffs_n)/len(diffs_n):.6f}")
        print(f"Average abs diff H: {sum(diffs_h)/len(diffs_h):.6f}")

    # Plot histograms
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        os.makedirs(os.path.dirname(args.histogram_png), exist_ok=True)
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))

        # Helper to plot one histogram
        def plot_hist(ax, data: List[float], title: str, bins=None):
            if not data:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax.set_title(title)
                return
            ax.hist(data, bins=(bins or 20), color='#4C72B0', edgecolor='black', alpha=0.8)
            ax.axvline(0.0, color='red', linestyle='--', linewidth=1)
            ax.set_title(title)
            ax.set_xlabel('Offset (ppm)')
            ax.set_ylabel('Count')

        # Compute common x-limits and bins for alignment
        # N offsets
        n_all = grid_n_values + mean_n_values
        if n_all:
            n_min, n_max = min(n_all), max(n_all)
            if n_min == n_max:
                n_min -= 0.01
                n_max += 0.01
            binsN = np.linspace(n_min, n_max, 21)
        else:
            n_min = n_max = 0.0
            binsN = 20
        # H offsets
        h_all = grid_h_values + mean_h_values
        if h_all:
            h_min, h_max = min(h_all), max(h_all)
            if h_min == h_max:
                h_min -= 0.01
                h_max += 0.01
            binsH = np.linspace(h_min, h_max, 21)
        else:
            h_min = h_max = 0.0
            binsH = 20

        plot_hist(axes[0,0], grid_n_values, 'Grid N offsets', bins=binsN)
        plot_hist(axes[1,0], mean_n_values, 'Mean N offsets', bins=binsN)
        axes[0,0].set_xlim(n_min, n_max)
        axes[1,0].set_xlim(n_min, n_max)

        plot_hist(axes[0,1], grid_h_values, 'Grid H offsets', bins=binsH)
        plot_hist(axes[1,1], mean_h_values, 'Mean H offsets', bins=binsH)
        axes[0,1].set_xlim(h_min, h_max)
        axes[1,1].set_xlim(h_min, h_max)

        fig.suptitle('Offset distributions by method')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(args.histogram_png, dpi=200)
        plt.close(fig)
        print(f"Histograms saved: {args.histogram_png}")
    except Exception as e:
        print(f"Warning: Could not generate histograms: {e}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())


