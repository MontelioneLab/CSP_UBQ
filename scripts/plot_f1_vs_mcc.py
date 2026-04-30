#!/usr/bin/env python3
"""
Plot F1 score vs MCC for per-system confusion data.

Reads outputs/confusion_matrix_per_system.csv and creates a scatter plot
of F1 score vs MCC for each system_id.

By default, saves the figure to outputs/f1_vs_mcc.png.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot F1 score vs MCC from confusion_matrix_per_system.csv."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv (default: outputs/confusion_matrix_per_system.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs") / "f1_vs_mcc.png",
        help="Path to output image file (default: outputs/f1_vs_mcc.png).",
    )
    return parser.parse_args(argv)


def load_holo_pdb_ids_from_targets_csv(path: Path) -> Set[str]:
    """Lowercase holo_pdb values from a CSP_UBQ-style CSV."""
    holo: Set[str] = set()
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or "holo_pdb" not in reader.fieldnames:
            raise ValueError(f"Missing 'holo_pdb' column in {path}")
        for row in reader:
            h = (row.get("holo_pdb") or "").strip().lower()
            if h:
                holo.add(h)
    return holo


def system_id_matches_holo_csv(
    system_id: str, holo_pdb_lower: Set[str]
) -> bool:
    """True if system_id (outputs folder style, e.g. 2mur_1) matches holo_pdb set from CSP CSV."""
    s = system_id.strip().lower()
    if not s:
        return False
    base = s.split("_", 1)[0]
    return s in holo_pdb_lower or base in holo_pdb_lower


def load_f1_mcc(
    path: Path,
    holo_pdb_filter: Optional[Set[str]] = None,
    *,
    allowed_system_ids: Optional[Set[str]] = None,
) -> Tuple[List[float], List[float]]:
    """Load F1 and MCC columns as floats from the confusion-matrix CSV.

    If ``allowed_system_ids`` is set, only rows whose ``system_id`` (case-insensitive)
    matches one of those strings are included — use resolved ``outputs/<dir>``
    basenames from :mod:`scripts.target_resolution`.

    Otherwise if ``holo_pdb_filter`` is set, only rows whose system_id matches an entry
    in the set (or shares the same base PDB before '_') are included.
    """
    f1_values: List[float] = []
    mcc_values: List[float] = []

    allowed_lower: Optional[Set[str]] = None
    if allowed_system_ids is not None:
        allowed_lower = {s.strip().lower() for s in allowed_system_ids if s.strip()}

    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if allowed_lower is not None:
                sid = (row.get("system_id") or "").strip()
                if not sid or sid.lower() not in allowed_lower:
                    continue
            elif holo_pdb_filter is not None:
                sid = (row.get("system_id") or "").strip()
                if not system_id_matches_holo_csv(sid, holo_pdb_filter):
                    continue
            f1_str = row.get("f1_score")
            mcc_str = row.get("mcc")
            if not f1_str or not mcc_str:
                continue
            try:
                f1 = float(f1_str)
                mcc = float(mcc_str)
            except ValueError:
                continue
            f1_values.append(f1)
            mcc_values.append(mcc)

    return f1_values, mcc_values


def plot_f1_vs_mcc(f1_values: List[float], mcc_values: List[float], output: Path) -> None:
    """Create and save a scatter plot of F1 vs MCC with linear fit and R²."""
    if not f1_values or not mcc_values:
        raise ValueError("No valid F1/MCC values found to plot.")

    output.parent.mkdir(parents=True, exist_ok=True)

    x = np.array(mcc_values)
    y = np.array(f1_values)

    # Linear fit: y = slope * x + intercept
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = slope * x_line + intercept

    plt.figure(figsize=(6, 5))
    plt.scatter(mcc_values, f1_values, alpha=0.7)
    plt.plot(x_line, y_line, "r--", linewidth=2, label=f"Linear fit (R² = {r_squared:.3f})")
    plt.xlabel("MCC")
    plt.ylabel("F1 score")
    plt.title("Per-system F1 vs MCC")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main(argv=None) -> int:
    args = parse_args(argv)
    f1_vals, mcc_vals = load_f1_mcc(args.input)
    plot_f1_vs_mcc(f1_vals, mcc_vals, args.output)
    print(f"Wrote F1 vs MCC plot to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

