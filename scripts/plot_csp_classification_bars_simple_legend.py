#!/usr/bin/env python3
"""
Create CSP classification bar plot with a simpler legend.

Reads master_alignment.csv from a target output directory and generates the same
bar chart as csp_classification_bars_original.png, but with a simplified legend:
  - Color + quadrant (TP/FP/TN/FN) + count (e.g. "TP (13)")

Binding-site residues (TP and FN) are highlighted with a light gray background.

Usage:
  python scripts/plot_csp_classification_bars_simple_legend.py outputs/7jq8
  python scripts/plot_csp_classification_bars_simple_legend.py outputs/7jq8 -o outputs/7jq8/csp_classification_bars_simple.png
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    _HAS_PLT = True
except ImportError:
    _HAS_PLT = False

try:
    from .config import classification_colors
except Exception:
    import os as _os
    import sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors


def _get_colors():
    return {
        "TP": classification_colors.TP,
        "FP": classification_colors.FP,
        "TN": classification_colors.TN,
        "FN": classification_colors.FN,
    }


def plot_csp_classification_bars_simple_legend(
    target_dir: Path,
    output_path: Path,
) -> bool:
    """
    Create CSP classification bar plot with simple legend (color, quadrant, count).

    Returns True on success, False otherwise.
    """
    if not _HAS_PLT:
        print("matplotlib not available", file=sys.stderr)
        return False

    master_csv = target_dir / "master_alignment.csv"
    if not master_csv.exists():
        print(f"master_alignment.csv not found in {target_dir}", file=sys.stderr)
        return False

    residue_numbers = []
    csp_values = []
    classifications = []
    aa_labels = []

    with open(master_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cls = (row.get("classification") or "").strip()
            if cls not in ("TP", "FP", "TN", "FN"):
                continue
            csp_str = (row.get("csp_A") or "").strip()
            if not csp_str:
                continue
            try:
                csp_val = float(csp_str)
            except ValueError:
                continue
            holo_resi = row.get("holo_resi", "")
            holo_aa = (row.get("holo_aa") or "").strip()
            try:
                resi = int(float(holo_resi))
            except (ValueError, TypeError):
                continue

            residue_numbers.append(resi)
            csp_values.append(csp_val)
            classifications.append(cls)
            aa_labels.append(holo_aa[0] if holo_aa else "")

    if not residue_numbers:
        print("No valid rows with csp_A and classification", file=sys.stderr)
        return False

    colors = _get_colors()
    bar_colors = [colors[cls] for cls in classifications]

    # Threshold: min CSP among significant (for display)
    sig_csps = [csp_values[i] for i, c in enumerate(classifications) if c in ("TP", "FP")]
    threshold = min(sig_csps) if sig_csps else 0.0

    plt.figure(figsize=(14, 8))
    ax = plt.gca()

    # Light gray backing for binding-site residues (TP and FN)
    bar_width = 0.8
    for x, cls in zip(residue_numbers, classifications):
        if cls in ("TP", "FN"):
            ax.axvspan(x - bar_width / 2, x + bar_width / 2, color="lightgray", alpha=0.5, zorder=0)

    plt.bar(residue_numbers, csp_values, color=bar_colors, alpha=0.8, edgecolor="black", linewidth=0.5, zorder=1)
    plt.axhline(y=threshold, color="black", linestyle="--", linewidth=2, alpha=0.8, label=f"Threshold = {threshold:.3f}")

    tp_count = classifications.count("TP")
    fp_count = classifications.count("FP")
    tn_count = classifications.count("TN")
    fn_count = classifications.count("FN")

    # Simple legend: color + quadrant + count
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors["TP"], alpha=0.8, label=f"TP ({tp_count})"),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors["FP"], alpha=0.8, label=f"FP ({fp_count})"),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors["TN"], alpha=0.8, label=f"TN ({tn_count})"),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors["FN"], alpha=0.8, label=f"FN ({fn_count})"),
    ]
    plt.legend(handles=legend_elements, loc="upper left", fontsize=28)

    ax = plt.gca()
    ax.set_xticks(residue_numbers)
    ax.set_xticklabels(aa_labels, fontsize=14)
    ax.tick_params(axis="x", which="major", pad=2, labelsize=14)
    ax.tick_params(axis="y", which="major", labelsize=22)
    plt.grid(True, alpha=0.3, axis="y")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_path}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create CSP classification bar plot with simple legend (TP/FP/TN/FN + count)."
    )
    parser.add_argument(
        "target_dir",
        type=Path,
        help="Target output directory (e.g. outputs/7jq8) containing master_alignment.csv",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output path (default: <target_dir>/csp_classification_bars_simple_legend.png)",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    target_dir = args.target_dir if args.target_dir.is_absolute() else project_root / args.target_dir

    if args.output:
        output_path = args.output if args.output.is_absolute() else project_root / args.output
    else:
        output_path = target_dir / "csp_classification_bars_simple_legend.png"

    success = plot_csp_classification_bars_simple_legend(target_dir, output_path)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
