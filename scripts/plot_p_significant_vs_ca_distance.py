#!/usr/bin/env python3
"""P(significant CSP | min interchain CA–CA distance) stacked histogram.

Residue-pooled estimate across selected targets. Bar height is

    100 * N_significant_in_bin / N_recorded_CSP_in_bin

with the significant numerator split into CSP z-score stacks:

    0 ≤ z < 1,  1 ≤ z < 2,  2 ≤ z < 3,  z ≥ 3

Default targets: data/CSP_UBQ_ph0.5_temp5C.csv (140 systems).
Default significance: primary ``significant`` (max(cleaned mean, 0.05)).

Usage:
    python scripts/plot_p_significant_vs_ca_distance.py \\
      --outputs outputs \\
      --targets-csv data/CSP_UBQ_ph0.5_temp5C.csv \\
      --output outputs/p_significant_vs_ca_distance_ph05.png
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.merge_csv import filter_recorded_csp_dataframe, parse_optional_bool
    from scripts.target_resolution import load_target_rows, resolve_target_rows

CSP_COLUMN = "csp_A"
CSP_Z_COLUMN = "csp_z"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
DEFAULT_SIGNIFICANT_COLUMN = "significant"
DISTANCE_XLABEL = r"Minimum interchain $C_\alpha$–$C_\alpha$ distance (Å)"
PROB_YLABEL = r"$P$(significant CSP $\mid$ distance) (%)"
DEFAULT_MAX_DISTANCE_A = 40.0

# Stack order bottom → top
Z_BANDS: Sequence[Tuple[str, str, Optional[float], Optional[float]]] = (
    # key, legend label, lo (inclusive), hi (exclusive); hi=None means ≥ lo
    ("z_0_1", r"$0 \leq z < 1$", 0.0, 1.0),
    ("z_1_2", r"$1 \leq z < 2$", 1.0, 2.0),
    ("z_2_3", r"$2 \leq z < 3$", 2.0, 3.0),
    ("z_ge_3", r"$z \geq 3$", 3.0, None),
)

Z_BAND_COLORS = {
    "z_0_1": "#c6dbef",
    "z_1_2": "#6baed6",
    "z_2_3": "#2171b5",
    "z_ge_3": "#08306b",
}


def _as_bool(value: object) -> bool:
    parsed = parse_optional_bool(value)
    return bool(parsed)


def _resolve_selected_dirs(
    outputs_dir: Path,
    targets_csv: Optional[Path],
    targets_str: Optional[str],
) -> Optional[Set[str]]:
    if targets_csv is None and not targets_str:
        return None
    extras: List[str] = []
    if targets_str:
        extras = [t for t in targets_str.split(",") if t.strip()]
    rows = load_target_rows(targets_csv, extra_holo_pdbs=extras)
    if not rows:
        return set()
    paths = resolve_target_rows(rows, outputs_dir)
    return {p.name for p in paths}


def _z_band_key(z: float) -> Optional[str]:
    if not np.isfinite(z):
        return None
    for key, _label, lo, hi in Z_BANDS:
        assert lo is not None
        if hi is None:
            if z >= lo:
                return key
        elif lo <= z < hi:
            return key
    return None


def collect_distance_significance(
    outputs_dir: Path,
    selected_dir_names: Optional[Set[str]] = None,
    *,
    significant_column: str = DEFAULT_SIGNIFICANT_COLUMN,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (distances, is_significant, csp_z, n_targets) for recorded-CSP residues.

    ``csp_z`` is NaN when missing.
    """
    distances: List[float] = []
    significance: List[bool] = []
    csp_z_vals: List[float] = []
    n_targets = 0

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if selected_dir_names is not None and target_name not in selected_dir_names:
            continue

        df = pd.read_csv(alignment_path)
        required = {significant_column, CA_DISTANCE_COLUMN}
        missing = [c for c in required if c not in df.columns]
        if missing:
            print(f"[WARN] Skipping {alignment_path}: missing {missing}")
            continue

        df = filter_recorded_csp_dataframe(
            df, csp_column=CSP_COLUMN, significant_column=significant_column
        )
        if df.empty:
            continue

        dist = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")
        sig = df[significant_column].map(_as_bool)
        if CSP_Z_COLUMN in df.columns:
            z = pd.to_numeric(df[CSP_Z_COLUMN], errors="coerce")
        else:
            z = pd.Series(np.nan, index=df.index)

        valid = dist.notna()
        if not valid.any():
            continue

        distances.extend(dist[valid].astype(float).tolist())
        significance.extend(sig[valid].astype(bool).tolist())
        csp_z_vals.extend(z[valid].astype(float).tolist())
        n_targets += 1

    return (
        np.asarray(distances, dtype=float),
        np.asarray(significance, dtype=bool),
        np.asarray(csp_z_vals, dtype=float),
        n_targets,
    )


def compute_stacked_bin_probabilities(
    distances: np.ndarray,
    significance: np.ndarray,
    csp_z: np.ndarray,
    *,
    bin_width: float = 1.0,
    max_distance: Optional[float] = None,
) -> pd.DataFrame:
    """Per distance bin: N_total and significant counts/percents by z-band."""
    if distances.size == 0:
        raise ValueError("No CA-distance / significance data found for selected targets.")
    if bin_width <= 0:
        raise ValueError("bin_width must be positive.")

    data_max = float(np.max(distances)) if max_distance is None else float(max_distance)
    n_bins = max(1, int(math.ceil((data_max + 1e-12) / bin_width)))
    edges = np.arange(0.0, (n_bins + 1) * bin_width, bin_width, dtype=float)

    if max_distance is not None:
        keep = distances <= max_distance
        distances = distances[keep]
        significance = significance[keep]
        csp_z = csp_z[keep]

    bin_idx = np.digitize(distances, edges[1:-1], right=False)
    rows = []
    for i in range(n_bins):
        left = float(edges[i])
        right = float(edges[i + 1])
        mask = bin_idx == i
        n_total = int(mask.sum())
        if n_total == 0:
            continue

        counts = {key: 0 for key, *_ in Z_BANDS}
        sig_mask = mask & significance
        for z in csp_z[sig_mask]:
            key = _z_band_key(float(z))
            if key is not None:
                counts[key] += 1

        n_sig = sum(counts.values())
        row = {
            "bin_left": left,
            "bin_right": right,
            "bin_center": 0.5 * (left + right),
            "n_total": n_total,
            "n_significant": n_sig,
            "p_significant": n_sig / n_total,
        }
        for key, *_ in Z_BANDS:
            row[f"n_{key}"] = counts[key]
            row[f"pct_{key}"] = 100.0 * counts[key] / n_total
        rows.append(row)
    return pd.DataFrame(rows)


def plot_p_significant_vs_distance(
    bin_df: pd.DataFrame,
    out_path: Path,
    *,
    n_residues: int,
    n_targets: int,
    significant_column: str,
    bin_width: float,
    dpi: int = 300,
    fig_width: float = 10.0,
    fig_height: float = 5.5,
) -> None:
    if bin_df.empty:
        raise ValueError("No non-empty distance bins to plot.")

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    bottoms = np.zeros(len(bin_df), dtype=float)
    x = bin_df["bin_left"].to_numpy(dtype=float)
    width = bin_width * 0.95

    for key, label, _lo, _hi in Z_BANDS:
        heights = bin_df[f"pct_{key}"].to_numpy(dtype=float)
        ax.bar(
            x,
            heights,
            width=width,
            bottom=bottoms,
            align="edge",
            color=Z_BAND_COLORS[key],
            edgecolor="black",
            linewidth=0.4,
            alpha=0.95,
            label=label,
        )
        bottoms = bottoms + heights

    total_pct = bottoms
    ax.set_xlabel(DISTANCE_XLABEL, fontsize=14)
    ax.set_ylabel(PROB_YLABEL, fontsize=14)
    y_max = float(total_pct.max()) if total_pct.size else 0.0
    ax.set_ylim(0.0, min(100.0, y_max * 1.15 + 2.0))
    ax.set_xlim(0.0, float(bin_df["bin_right"].max()))
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, axis="y", alpha=0.35, linestyle="--", linewidth=0.6)
    ax.legend(title="Significant CSP z-score", fontsize=10, title_fontsize=11, loc="best")
    ax.set_title(
        f"Pooled over {n_targets} targets, {n_residues} residues "
        f"(significance={significant_column})",
        fontsize=12,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_bin_csv(bin_df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bin_df.to_csv(out_path, index=False)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stacked histogram of P(significant CSP | min interchain CA–CA distance), "
            "with stacks by CSP z-score band."
        )
    )
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_REPO_ROOT / "outputs",
        help="Outputs directory containing */master_alignment.csv "
        "(default: outputs).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=_REPO_ROOT / "data/CSP_UBQ_ph0.5_temp5C.csv",
        help="Targets CSV (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        help="Optional comma-separated holo_pdb list (merged with --targets-csv).",
    )
    parser.add_argument(
        "--significant-column",
        default=DEFAULT_SIGNIFICANT_COLUMN,
        help="Significance column in master_alignment.csv (default: significant).",
    )
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument(
        "--max-distance",
        type=float,
        default=DEFAULT_MAX_DISTANCE_A,
        help="Cap X-axis / included distances at this Å value (default: 40).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: <outputs>/p_significant_vs_ca_distance_ph05.png).",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Optional bin-stats CSV (default: sibling of --output with .csv).",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--fig-width", type=float, default=10.0)
    parser.add_argument("--fig-height", type=float, default=5.5)
    return parser.parse_args(list(argv))


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs if args.outputs.is_absolute() else _REPO_ROOT / args.outputs
    targets_csv = args.targets_csv
    if targets_csv is not None and not targets_csv.is_absolute():
        targets_csv = _REPO_ROOT / targets_csv

    selected = _resolve_selected_dirs(outputs_dir, targets_csv, args.targets)
    if selected is not None:
        print(f"Resolved {len(selected)} target directories from targets list.")

    distances, significance, csp_z, n_targets = collect_distance_significance(
        outputs_dir,
        selected,
        significant_column=args.significant_column,
    )
    print(
        f"Collected {len(distances)} residues with recorded CSP from {n_targets} targets."
    )

    bin_df = compute_stacked_bin_probabilities(
        distances,
        significance,
        csp_z,
        bin_width=args.bin_width,
        max_distance=args.max_distance,
    )

    out_png = args.output or (
        outputs_dir / "p_significant_vs_ca_distance_ph05.png"
    )
    if not out_png.is_absolute():
        out_png = _REPO_ROOT / out_png
    out_csv = args.output_csv
    if out_csv is None:
        out_csv = out_png.with_suffix(".csv")
    elif not out_csv.is_absolute():
        out_csv = _REPO_ROOT / out_csv

    # Count residues used in the plotted distance range for the title.
    if args.max_distance is not None:
        n_residues_plot = int(np.sum(distances <= args.max_distance))
    else:
        n_residues_plot = int(len(distances))

    plot_p_significant_vs_distance(
        bin_df,
        out_png,
        n_residues=n_residues_plot,
        n_targets=n_targets,
        significant_column=args.significant_column,
        bin_width=args.bin_width,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )
    write_bin_csv(bin_df, out_csv)

    print(f"Wrote {out_png}")
    print(f"Wrote {out_csv}")
    print(
        "bin_left  n_total  n_sig  pct_total  "
        "pct_0_1  pct_1_2  pct_2_3  pct_ge_3"
    )
    for row in bin_df.itertuples(index=False):
        stack_sum = row.pct_z_0_1 + row.pct_z_1_2 + row.pct_z_2_3 + row.pct_z_ge_3
        print(
            f"{row.bin_left:7.1f}  {row.n_total:7d}  {row.n_significant:5d}  "
            f"{stack_sum:8.2f}  "
            f"{row.pct_z_0_1:7.2f}  {row.pct_z_1_2:7.2f}  "
            f"{row.pct_z_2_3:7.2f}  {row.pct_z_ge_3:7.2f}"
        )
        # Stack percents must equal 100 * p_significant
        expected = 100.0 * row.p_significant
        if abs(stack_sum - expected) > 1e-6:
            raise AssertionError(
                f"Stack percent mismatch at bin {row.bin_left}: "
                f"{stack_sum} vs {expected}"
            )

    near = bin_df[bin_df["bin_right"] <= 6.0]
    far = bin_df[bin_df["bin_left"] >= 20.0]
    if not near.empty and not far.empty:
        p_near = (
            near["n_significant"].sum() / near["n_total"].sum()
            if near["n_total"].sum()
            else float("nan")
        )
        p_far = (
            far["n_significant"].sum() / far["n_total"].sum()
            if far["n_total"].sum()
            else float("nan")
        )
        print(f"Pooled P(sig | d<6 Å) = {p_near:.3f}")
        print(f"Pooled P(sig | d≥20 Å) = {p_far:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
