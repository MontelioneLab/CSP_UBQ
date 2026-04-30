#!/usr/bin/env python3
"""
Boxplot of per-target F1 scores for 1D H/N/CA CSPs with paired significance
annotations.

For each target directory under ``outputs/``, the F1 scores for 1D H, N, and
CA CSPs are computed (reusing :func:`collect_1d_f1_results` from
:mod:`analyze_targets_single_atom_shifts`). The figure is then restricted to
the intersection of targets that yield a valid F1 for all three atom types -
which coincides with the subset of "targets with Ca shifts", since CA F1
requires ``CA_apo``/``CA_holo`` to be populated.

Inter-group comparisons are performed with paired Wilcoxon signed-rank tests
(H vs N, H vs CA, N vs CA) and the three raw p-values are adjusted with
Holm-Bonferroni. Significance is overlaid on the boxplot (ns / * / ** / ***)
and the pairwise results are written to a companion CSV.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Support running as a script or module.
try:
    from .analyze_targets_single_atom_shifts import (
        TargetResult,
        collect_1d_f1_results,
        load_allowed_targets,
    )
except Exception:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from scripts.analyze_targets_single_atom_shifts import (  # type: ignore
        TargetResult,
        collect_1d_f1_results,
        load_allowed_targets,
    )


ATOM_ORDER: Tuple[str, str, str] = ("H", "N", "CA")
ATOM_COLORS = {
    "H": "#66c2a5",
    "N": "#fc8d62",
    "CA": "#8da0cb",
}


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot F1 score distributions for 1D H/N/CA CSPs on the subset of "
            "targets with CA shifts, annotated with paired Wilcoxon signed-rank "
            "significance (Holm-Bonferroni corrected)."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=None,
        help=(
            "Optional CSV containing a 'holo_pdb' column to filter targets "
            "before computing F1 scores."
        ),
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF_f1_1d_boxplot.png",
        help="Destination for the boxplot figure (default: %(default)s).",
    )
    parser.add_argument(
        "--stats-csv",
        type=Path,
        default=None,
        help=(
            "Destination for the pairwise Wilcoxon stats CSV. "
            "If omitted, defaults to sibling of --output-image with the same "
            "stem plus '_stats.csv'."
        ),
    )
    parser.add_argument(
        "--min-ca-coverage",
        type=float,
        default=0.5,
        help=(
            "Minimum fraction of residues in 1d_analysis.csv with both CA_apo "
            "and CA_holo populated for a target to be considered as having "
            "CA shifts (default: %(default)s, i.e. >50%% coverage)."
        ),
    )
    return parser.parse_args(list(argv))


def _ca_coverage_for_target(target_dir: Path) -> Optional[float]:
    """Return the fraction of residues with both CA_apo and CA_holo populated.

    Returns ``None`` if ``1d_analysis.csv`` is missing, unreadable, empty, or
    lacks the required CA columns.
    """
    one_d_path = target_dir / "1d_analysis.csv"
    if not one_d_path.exists():
        return None
    try:
        df = pd.read_csv(one_d_path)
    except Exception:
        return None
    if df.empty or "CA_apo" not in df.columns or "CA_holo" not in df.columns:
        return None
    total = int(len(df))
    if total == 0:
        return None
    both = int((df["CA_apo"].notna() & df["CA_holo"].notna()).sum())
    return both / total


def _targets_with_ca_coverage(
    outputs_dir: Path,
    min_coverage: float,
    allowed_targets: Optional[dict],
) -> Tuple[set, dict]:
    """Enumerate target subdirectories and return the set meeting ``min_coverage``.

    Also returns a mapping of ``target -> coverage`` for reporting.
    """
    passing: set = set()
    coverages: dict = {}
    if not outputs_dir.exists():
        return passing, coverages
    for path in sorted(outputs_dir.iterdir()):
        if not path.is_dir() or path.name.startswith("."):
            continue
        if allowed_targets is not None and path.name not in allowed_targets:
            continue
        cov = _ca_coverage_for_target(path)
        if cov is None:
            continue
        coverages[path.name] = cov
        if cov > min_coverage:
            passing.add(path.name)
    return passing, coverages


def _build_paired_arrays(
    results_H: Sequence[TargetResult],
    results_N: Sequence[TargetResult],
    results_CA: Sequence[TargetResult],
) -> Tuple[List[str], np.ndarray, np.ndarray, np.ndarray]:
    """Intersect on target names and return paired F1 arrays for H, N, CA."""
    by_target_h = {r.target: r.f1 for r in results_H}
    by_target_n = {r.target: r.f1 for r in results_N}
    by_target_ca = {r.target: r.f1 for r in results_CA}

    common = sorted(set(by_target_h) & set(by_target_n) & set(by_target_ca))
    h = np.array([by_target_h[t] for t in common], dtype=float)
    n = np.array([by_target_n[t] for t in common], dtype=float)
    ca = np.array([by_target_ca[t] for t in common], dtype=float)
    return common, h, n, ca


def _holm_bonferroni(pvals: Sequence[float]) -> List[float]:
    """Holm-Bonferroni step-down correction for a small set of p-values."""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running_max = 0.0
    for rank, idx in enumerate(order):
        scaled = (m - rank) * pvals[idx]
        scaled = min(scaled, 1.0)
        running_max = max(running_max, scaled)
        adj[idx] = running_max
    return adj


def _significance_label(p: float) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _paired_wilcoxon_table(
    h: np.ndarray,
    n: np.ndarray,
    ca: np.ndarray,
) -> pd.DataFrame:
    """Compute paired Wilcoxon signed-rank for the three atom-type pairs."""
    pairs: List[Tuple[str, str, np.ndarray, np.ndarray]] = [
        ("H", "N", h, n),
        ("H", "CA", h, ca),
        ("N", "CA", n, ca),
    ]

    rows: List[dict] = []
    raw_pvalues: List[float] = []
    for a, b, x, y in pairs:
        diff = x - y
        n_pairs = int(np.sum(np.isfinite(diff)))
        # If all paired differences are zero, wilcoxon raises; guard it.
        if n_pairs == 0 or np.all(diff[np.isfinite(diff)] == 0):
            stat = float("nan")
            pval = 1.0
        else:
            try:
                res = wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
                stat = float(res.statistic)
                pval = float(res.pvalue)
            except ValueError:
                stat = float("nan")
                pval = 1.0
        raw_pvalues.append(pval)
        rows.append(
            {
                "group_a": a,
                "group_b": b,
                "n": n_pairs,
                "statistic": stat,
                "p_raw": pval,
            }
        )

    adj = _holm_bonferroni(raw_pvalues)
    for row, p_adj in zip(rows, adj):
        row["p_adj_holm"] = float(p_adj)
        row["signif"] = _significance_label(p_adj)

    return pd.DataFrame(rows, columns=["group_a", "group_b", "n", "statistic", "p_raw", "p_adj_holm", "signif"])


def _atom_positions() -> dict:
    return {atom: i + 1 for i, atom in enumerate(ATOM_ORDER)}


def _render_boxplot(
    h: np.ndarray,
    n: np.ndarray,
    ca: np.ndarray,
    stats_df: pd.DataFrame,
    output_image: Path,
    n_targets: int,
) -> None:
    """Draw the annotated boxplot and save to ``output_image``."""
    data_for_plot = [h, n, ca]
    positions = list(_atom_positions().values())

    fig, ax = plt.subplots(figsize=(7, 6))

    bp = ax.boxplot(
        data_for_plot,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        showmeans=True,
        meanline=True,
        medianprops=dict(color="black", linewidth=1.5),
        meanprops=dict(color="black", linestyle="--", linewidth=1.2),
    )
    for patch, atom in zip(bp["boxes"], ATOM_ORDER):
        patch.set_facecolor(ATOM_COLORS[atom])
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(seed=0)
    for pos, values, atom in zip(positions, data_for_plot, ATOM_ORDER):
        if len(values) == 0:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(
            np.full_like(values, pos, dtype=float) + jitter,
            values,
            color=ATOM_COLORS[atom],
            edgecolor="black",
            linewidth=0.3,
            alpha=0.55,
            s=16,
            zorder=2,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(list(ATOM_ORDER))
    ax.set_ylabel("F1 Score", fontsize=12)
    ax.set_xlabel("Atom Type (1D CSP)", fontsize=12)
    ax.set_title(
        f"F1 score distributions for 1D single-atom CSPs (n = {n_targets} targets)",
        fontsize=13,
        fontweight="bold",
        pad=28,
    )
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.set_ylim(0.0, 1.0)

    _annotate_significance(ax, stats_df)

    fig.tight_layout()
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _annotate_significance(ax: plt.Axes, stats_df: pd.DataFrame) -> None:
    """Draw stacked significance bars over the boxes (within y in [0, 1])."""
    pos = _atom_positions()
    # Bracket stack fits below y=1 when axis is [0, 1] (was above 1 for ylim 1.3).
    base_y = 0.82
    step = 0.045
    bar_height = 0.012

    # Order bars so shorter spans sit beneath the longer one.
    ordering = [("H", "N"), ("N", "CA"), ("H", "CA")]

    for level, (a, b) in enumerate(ordering):
        row = stats_df[(stats_df["group_a"] == a) & (stats_df["group_b"] == b)]
        if row.empty:
            continue
        p_adj = float(row["p_adj_holm"].iloc[0])
        label = _significance_label(p_adj)
        if label == "ns":
            annotation = f"ns (p = {p_adj:.2g})"
        else:
            annotation = f"{label} (p = {p_adj:.2g})"

        x1 = pos[a]
        x2 = pos[b]
        y = base_y + level * step
        ax.plot(
            [x1, x1, x2, x2],
            [y, y + bar_height, y + bar_height, y],
            lw=1.1,
            color="black",
            clip_on=False,
        )
        ax.text(
            (x1 + x2) / 2.0,
            y + bar_height + 0.005,
            annotation,
            ha="center",
            va="bottom",
            fontsize=9,
            clip_on=False,
        )


def _resolve_stats_csv(output_image: Path, stats_csv: Optional[Path]) -> Path:
    if stats_csv is not None:
        return stats_csv
    return output_image.with_name(f"{output_image.stem}_stats.csv")


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    outputs_dir: Path = args.outputs_dir.resolve()
    if not outputs_dir.exists():
        print(f"Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1

    allowed_targets = load_allowed_targets(args.targets_csv.resolve() if args.targets_csv else None)

    ca_passing, ca_coverages = _targets_with_ca_coverage(
        outputs_dir,
        min_coverage=float(args.min_ca_coverage),
        allowed_targets=allowed_targets,
    )
    print(
        f"{len(ca_passing)} targets pass CA coverage > {args.min_ca_coverage:.0%} "
        f"(of {len(ca_coverages)} candidates inspected)"
    )
    if not ca_passing:
        print(
            "No targets satisfy the CA coverage threshold; nothing to plot.",
            file=sys.stderr,
        )
        return 1

    # Restrict F1 collection to the CA-coverage-passing subset so all three
    # atom-type distributions are computed on the same set of targets.
    ca_allowed = {t: True for t in ca_passing}
    results_H, results_N, results_CA = collect_1d_f1_results(outputs_dir, ca_allowed)
    if not (results_H and results_N and results_CA):
        print(
            "No 1D F1 results found for one or more atom types; cannot render boxplot.",
            file=sys.stderr,
        )
        return 1

    common, h, n, ca = _build_paired_arrays(results_H, results_N, results_CA)
    if len(common) == 0:
        print(
            "No targets produce F1 scores for all of H, N, and CA; nothing to plot.",
            file=sys.stderr,
        )
        return 1

    stats_df = _paired_wilcoxon_table(h, n, ca)

    output_image: Path = args.output_image
    stats_csv = _resolve_stats_csv(output_image, args.stats_csv)

    _render_boxplot(h, n, ca, stats_df, output_image, n_targets=len(common))

    stats_csv.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(stats_csv, index=False)

    print(f"n = {len(common)} targets with paired H/N/CA F1 scores")
    print(f"Wrote {output_image}")
    print(f"Wrote {stats_csv}")
    for _, row in stats_df.iterrows():
        print(
            f"  {row['group_a']} vs {row['group_b']}: "
            f"W={row['statistic']:.3g}, p_raw={row['p_raw']:.3g}, "
            f"p_adj(Holm)={row['p_adj_holm']:.3g} ({row['signif']})"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
