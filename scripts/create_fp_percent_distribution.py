#!/usr/bin/env python3
"""
Per-target FP% distributions using the primary CSP significance cutoff.

Two complementary rates (same FP residue definition):

1. **All-residue FP%** = 100 × FP / (TP+FP+TN+FN)
2. **Among-significant FP%** = 100 × FP / (FP+TP)
   (fraction of significant CSPs that are allosteric; targets with TP+FP = 0
    are omitted from that metric)

Significant CSPs use the ``significant`` column only (iterative outlier
removal, then ``max(cleaned mean, 0.05 ppm)`` with default
``significance_z=0``). Does **not** use 1SD/2SD, percentile, or fixed-ppm
cutoffs.

The top three and bottom three targets per metric are labeled as
``holo_pdb/apo_bmrb``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from .analyze_targets import (
        AlignmentParsingError,
        PREDICTOR_COLUMNS,
        load_alignment,
    )
    from .config import classification_colors
    from .target_resolution import load_target_rows, resolve_target_rows
except ImportError:
    from analyze_targets import (  # type: ignore
        AlignmentParsingError,
        PREDICTOR_COLUMNS,
        load_alignment,
    )
    from config import classification_colors  # type: ignore
    from target_resolution import load_target_rows, resolve_target_rows  # type: ignore

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SIGNIFICANT_COL = "significant"  # max(cleaned mean, 0.05 ppm) primary cutoff

# Metric column → axis / print labels
_METRICS = {
    "fp_pct": {
        "short": "FP / all residues",
        "xlabel": r"100 $\times$ FP / (TP+FP+TN+FN)",
        "title": "All-residue FP%",
    },
    "fp_of_sig_pct": {
        "short": "FP / (FP+TP)",
        "xlabel": r"100 $\times$ FP / (FP+TP)",
        "title": "Among-significant FP%",
    },
}


def _pair_label(apo_pdb: str, apo_bmrb: str, holo_pdb: str) -> str:
    """Format target label as ``holo_pdb/apo_bmrb``."""
    holo = (holo_pdb or "").strip().upper()
    apo = (apo_bmrb or "").strip()
    return f"{holo}/{apo}"


def _load_targets_meta(targets_csv: Path) -> Dict[str, dict]:
    """Map system_id → apo/holo metadata from the targets CSV."""
    meta: Dict[str, dict] = {}
    with targets_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = {fn.lower(): fn for fn in (reader.fieldnames or [])}
        for raw in reader:
            holo = (raw.get(fields.get("holo_pdb", "holo_pdb"), "") or "").strip()
            apo_bmrb = (raw.get(fields.get("apo_bmrb", "apo_bmrb"), "") or "").strip()
            if not holo or not apo_bmrb:
                continue
            system_id = f"{holo.upper()}_{apo_bmrb}"
            apo_pdb = ""
            if "apo_pdb" in fields:
                apo_pdb = (raw.get(fields["apo_pdb"], "") or "").strip()
            meta[system_id] = {
                "apo_pdb": apo_pdb,
                "apo_bmrb": apo_bmrb,
                "holo_pdb": holo.upper(),
                "holo_bmrb": (raw.get(fields.get("holo_bmrb", "holo_bmrb"), "") or "").strip()
                if "holo_bmrb" in fields
                else "",
            }
    return meta


def compute_fp_pct_cleaned_mean(target_dir: Path) -> Optional[dict]:
    """Both FP rates for one target using primary ``significant`` only."""
    alignment_path = target_dir / "master_alignment.csv"
    if not alignment_path.is_file():
        return None
    try:
        df = load_alignment(alignment_path)
    except AlignmentParsingError as exc:
        print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
        return None

    if _SIGNIFICANT_COL not in df.columns:
        return None
    available = [c for c in PREDICTOR_COLUMNS if c in df.columns]
    if not available:
        return None

    predicted = df[available].any(axis=1)
    actual = df[_SIGNIFICANT_COL]
    total = int(len(df))
    if total == 0:
        return None

    n_fp = int((actual & ~predicted).sum())
    n_tp = int((actual & predicted).sum())
    n_fn = int((~actual & predicted).sum())
    n_tn = total - n_tp - n_fp - n_fn
    fp_pct = 100.0 * n_fp / total
    sig = n_tp + n_fp
    fp_of_sig_pct = (100.0 * n_fp / sig) if sig > 0 else float("nan")

    return {
        "system_id": target_dir.name,
        "n_tp": n_tp,
        "n_fp": n_fp,
        "n_tn": n_tn,
        "n_fn": n_fn,
        "n_total": total,
        "n_significant": sig,
        "fp_pct": fp_pct,
        "fp_of_sig_pct": fp_of_sig_pct,
    }


def collect_fp_pct_table(
    outputs_dir: Path,
    targets_csv: Path,
) -> pd.DataFrame:
    """Resolve targets and compute primary-cutoff FP rates per system."""
    rows = load_target_rows(targets_csv)
    target_paths = resolve_target_rows(rows, outputs_dir)
    meta = _load_targets_meta(targets_csv)

    records: List[dict] = []
    for target_dir in target_paths:
        metrics = compute_fp_pct_cleaned_mean(target_dir)
        if metrics is None:
            continue
        info = meta.get(target_dir.name, {})
        apo_pdb = info.get("apo_pdb", "")
        apo_bmrb = info.get("apo_bmrb", "")
        holo_pdb = info.get("holo_pdb", target_dir.name.split("_")[0])
        if not apo_bmrb and "_" in target_dir.name:
            apo_bmrb = target_dir.name.split("_", 1)[1]
        records.append(
            {
                **metrics,
                "apo_pdb": apo_pdb,
                "apo_bmrb": apo_bmrb,
                "holo_pdb": holo_pdb,
                "pair_label": _pair_label(apo_pdb, apo_bmrb, holo_pdb),
            }
        )

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values("fp_pct", ascending=False).reset_index(drop=True)


def summarize(series: pd.Series) -> Dict[str, float]:
    """Summary stats over finite values only."""
    s = series.dropna().astype(float)
    if s.empty:
        return {
            "n": 0.0,
            "mean": float("nan"),
            "std": float("nan"),
            "median": float("nan"),
            "q25": float("nan"),
            "q75": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
        }
    return {
        "n": float(len(s)),
        "mean": float(s.mean()),
        "std": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
        "median": float(s.median()),
        "q25": float(s.quantile(0.25)),
        "q75": float(s.quantile(0.75)),
        "min": float(s.min()),
        "max": float(s.max()),
    }


def _top_bottom_mask(values: pd.Series, *, n_extreme: int = 3) -> pd.Series:
    """Boolean mask for the top and bottom ``n_extreme`` finite values."""
    valid = values.dropna()
    mask = pd.Series(False, index=values.index)
    if valid.empty:
        return mask
    ranked_desc = valid.rank(method="first", ascending=False)
    n = len(valid)
    n_keep = min(n_extreme, n)
    mask.loc[valid.index] = (ranked_desc <= n_keep) | (ranked_desc > n - n_keep)
    return mask


def _annotate_top_bottom(
    ax: plt.Axes,
    df: pd.DataFrame,
    col: str,
    label_mask: pd.Series,
    jitter: np.ndarray,
) -> None:
    """Label top/bottom extremes as holo_pdb/apo_bmrb without overlapping text."""
    from matplotlib.transforms import blended_transform_factory

    labeled = df.loc[label_mask].dropna(subset=[col]).copy()
    if labeled.empty:
        return

    n_keep = min(3, len(labeled))
    bottom = labeled.nsmallest(n_keep, col).sort_values(col, ascending=True)
    top = labeled.nlargest(n_keep, col).sort_values(col, ascending=False)
    top = top.loc[~top.index.isin(bottom.index)]

    # x in data coords, y in axes fraction → evenly spaced labels that cannot collide.
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    def _draw(group: pd.DataFrame, *, x_text: float, y_fracs: Sequence[float], ha: str) -> None:
        for (idx, row), y_frac in zip(group.iterrows(), y_fracs):
            loc = df.index.get_loc(idx)
            x = float(jitter[loc]) if isinstance(loc, (int, np.integer)) else 0.0
            y = float(row[col])
            ax.annotate(
                row["pair_label"],
                xy=(x, y),
                xytext=(x_text, y_frac),
                textcoords=trans,
                fontsize=7.5,
                ha=ha,
                va="center",
                arrowprops=dict(
                    arrowstyle="-",
                    color="0.45",
                    lw=0.6,
                    shrinkA=0,
                    shrinkB=2,
                    connectionstyle="arc3,rad=0",
                ),
                clip_on=False,
            )

    # Bottom trio: left side, low → high in axes fraction.
    n_b = len(bottom)
    if n_b:
        y_fracs_b = [0.10 + i * 0.14 for i in range(n_b)]
        _draw(bottom, x_text=-0.72, y_fracs=y_fracs_b, ha="right")

    # Top trio: right side, high → low in axes fraction.
    n_t = len(top)
    if n_t:
        y_fracs_t = [0.90 - i * 0.14 for i in range(n_t)]
        _draw(top, x_text=0.72, y_fracs=y_fracs_t, ha="left")


def _plot_metric_row(
    axes: Sequence[plt.Axes],
    df: pd.DataFrame,
    col: str,
    *,
    n_extreme: int,
    color: str,
    legend_loc: str = "upper right",
) -> pd.DataFrame:
    """Draw histogram + labeled strip for one metric column; return labeled rows."""
    meta = _METRICS[col]
    valid = df.dropna(subset=[col]).copy().reset_index(drop=True)
    values = valid[col].to_numpy(dtype=float)
    stats = summarize(valid[col])
    label_mask = _top_bottom_mask(valid[col], n_extreme=n_extreme)

    ax_hist, ax_strip = axes

    xmax = max(55.0, float(np.ceil(np.nanmax(values) / 5.0) * 5.0)) if len(values) else 100.0
    bins = np.linspace(0, xmax, 13)
    ax_hist.hist(values, bins=bins, color=color, edgecolor="white", linewidth=0.8, alpha=0.85)
    if np.isfinite(stats["mean"]):
        ax_hist.axvline(
            stats["mean"], color="black", linestyle="-", linewidth=1.4,
            label=f"mean = {stats['mean']:.1f}%",
        )
    if np.isfinite(stats["median"]):
        ax_hist.axvline(
            stats["median"], color="black", linestyle="--", linewidth=1.2,
            label=f"median = {stats['median']:.1f}%",
        )
    ax_hist.set_xlabel(meta["xlabel"], fontsize=10)
    ax_hist.set_ylabel("Number of Targets", fontsize=10)
    ax_hist.legend(frameon=False, fontsize=8, loc=legend_loc)
    ax_hist.spines["top"].set_visible(False)
    ax_hist.spines["right"].set_visible(False)

    rng = np.random.default_rng(0 if col == "fp_pct" else 1)
    jitter = rng.uniform(-0.12, 0.12, size=len(values))
    ax_strip.boxplot(
        [values],
        positions=[0],
        widths=0.35,
        vert=True,
        whis=1.5,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(facecolor=color, alpha=0.25, edgecolor=color),
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="0.35"),
        capprops=dict(color="0.35"),
    )
    not_labeled = ~label_mask.to_numpy()
    ax_strip.scatter(
        jitter[not_labeled],
        values[not_labeled],
        s=26,
        c=color,
        alpha=0.55,
        edgecolors="none",
        zorder=2,
    )
    ax_strip.scatter(
        jitter[label_mask.to_numpy()],
        values[label_mask.to_numpy()],
        s=46,
        c=color,
        alpha=1.0,
        edgecolors="black",
        linewidths=0.7,
        zorder=3,
    )
    _annotate_top_bottom(ax_strip, valid, col, label_mask, jitter)
    # Room for left/right labels outside the strip cloud.
    ax_strip.set_xlim(-1.15, 1.15)
    y0, y1 = float(np.nanmin(values)), float(np.nanmax(values))
    pad = max(0.08 * (y1 - y0), 3.0)
    ax_strip.set_ylim(y0 - pad, y1 + pad)
    ax_strip.set_xticks([])
    ax_strip.set_ylabel(meta["xlabel"], fontsize=9)
    ax_strip.spines["top"].set_visible(False)
    ax_strip.spines["right"].set_visible(False)
    ax_strip.spines["bottom"].set_visible(False)

    return valid.loc[label_mask].sort_values(col, ascending=False)


def plot_fp_percent_distribution(
    df: pd.DataFrame,
    out_path: Path,
    *,
    n_extreme: int = 3,
) -> Dict[str, pd.DataFrame]:
    """Two-row figure: FP/(FP+TP) on top, all-residue FP% below. Return labeled rows per metric."""
    fp_color = classification_colors.FP
    among_color = "#7d3c98"

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 9.2),
        gridspec_kw={"width_ratios": [1.25, 1.15], "wspace": 0.32, "hspace": 0.32},
    )

    notable_sig = _plot_metric_row(
        axes[0],
        df,
        "fp_of_sig_pct",
        n_extreme=n_extreme,
        color=among_color,
        legend_loc="upper left",
    )
    notable_all = _plot_metric_row(
        axes[1], df, "fp_pct", n_extreme=n_extreme, color=fp_color, legend_loc="upper right",
    )

    for ax, label in zip(axes.flat, ("a", "b", "c", "d")):
        ax.text(
            -0.10,
            1.02,
            label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            color="black",
            clip_on=False,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return {"fp_pct": notable_all, "fp_of_sig_pct": notable_sig}


def _print_metric_block(df: pd.DataFrame, col: str, notable: pd.DataFrame) -> None:
    meta = _METRICS[col]
    stats = summarize(df[col])
    print(f"\n{meta['title']}  [{meta['short']}]")
    print("-" * 52)
    print(f"Targets: {int(stats['n'])}")
    print(f"Mean:     {stats['mean']:.2f}%")
    print(f"Median:   {stats['median']:.2f}%")
    print(f"Std:      {stats['std']:.2f}%")
    print(f"IQR:      {stats['q25']:.2f}% – {stats['q75']:.2f}%")
    if stats["n"] > 0:
        imin = df[col].idxmin()
        imax = df[col].idxmax()
        print(f"Min:      {stats['min']:.2f}%  ({df.loc[imin, 'pair_label']})")
        print(f"Max:      {stats['max']:.2f}%  ({df.loc[imax, 'pair_label']})")
    print("Notable apo/holo pairs:")
    for _, row in notable.iterrows():
        if col == "fp_pct":
            detail = f"{int(row['n_fp'])} FP / {int(row['n_total'])} residues"
        else:
            detail = f"{int(row['n_fp'])} FP / {int(row['n_significant'])} significant"
        print(f"  {row['pair_label']:28s}  {row[col]:5.2f}%  ({detail})")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Plot and tabulate per-target FP rates using the primary CSP cutoff "
            "(max of cleaned mean and 0.05 ppm): FP/(all residues) and FP/(FP+TP)."
        )
    )
    p.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs_cached_pre_nature_revisions"),
        help="Root with per-target master_alignment.csv (default: %(default)s).",
    )
    p.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="Targets CSV with apo_bmrb / holo_pdb (default: %(default)s).",
    )
    p.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "fp_percent_distribution_cleaned_mean.png",
        help="Output PNG path (PDF written alongside).",
    )
    p.add_argument(
        "--output-csv",
        type=Path,
        default=Path("figures") / "fp_percent_distribution_cleaned_mean.csv",
        help="Per-target FP rates table CSV.",
    )
    p.add_argument(
        "--n-extreme",
        type=int,
        default=3,
        help="Also label top/bottom N extremes in addition to IQR outliers (default: %(default)s).",
    )
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else _REPO_ROOT / args.targets_csv
    out_img = args.output_image if args.output_image.is_absolute() else _REPO_ROOT / args.output_image
    out_csv = args.output_csv if args.output_csv.is_absolute() else _REPO_ROOT / args.output_csv

    if not targets_csv.is_file():
        print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
        return 1

    df = collect_fp_pct_table(outputs_dir, targets_csv)
    if df.empty:
        print("No receptors with primary-cutoff classification data found.", file=sys.stderr)
        return 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    notable = plot_fp_percent_distribution(df, out_img, n_extreme=args.n_extreme)

    print("FP rate distributions — max(cleaned mean, 0.05 ppm) CSP cutoff")
    print("=" * 52)
    _print_metric_block(df, "fp_pct", notable["fp_pct"])
    _print_metric_block(df, "fp_of_sig_pct", notable["fp_of_sig_pct"])
    print()
    print(f"Wrote {out_img}")
    print(f"Wrote {out_img.with_suffix('.pdf')}")
    print(f"Wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
