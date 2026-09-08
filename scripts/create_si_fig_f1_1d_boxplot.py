#!/usr/bin/env python3
"""
Boxplot of per-target F1 scores for configurable 1D atom CSPs (default H/N/CA/HA)
with paired significance annotations.

For each target directory under ``outputs/``, the F1 scores for 1D H, N, CA,
and optionally HA CSPs are computed (reusing :func:`collect_1d_f1_results` from
:mod:`analyze_targets_single_atom_shifts`). Targets are gated by CA row
coverage (:func:`target_basenames_passing_ca_shift_coverage`, same rule as SI
Fig. S11 / S12). The plot uses the intersection of targets that yield a valid
F1 for every atom in the chosen ``atom_order`` (SI Fig. S15 uses H/N/CA only so
*n* matches the CA-gated cohort without requiring HA F1).

Inter-group comparisons use paired Wilcoxon signed-rank tests for all
unordered atom pairs, with Holm-Bonferroni–adjusted *p*-values overlaid on the
figure and written to a companion CSV.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Support running as a script or module.
try:
    from .analyze_targets_single_atom_shifts import (
        DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        TargetResult,
        collect_1d_f1_results,
        target_basenames_passing_ca_shift_coverage,
    )
    from .target_resolution import load_target_rows, resolve_target_rows
except Exception:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from scripts.analyze_targets_single_atom_shifts import (  # type: ignore
        DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        TargetResult,
        collect_1d_f1_results,
        target_basenames_passing_ca_shift_coverage,
    )
    from scripts.target_resolution import load_target_rows, resolve_target_rows  # type: ignore


ATOM_ORDER: Tuple[str, ...] = ("H", "N", "CA", "HA")
ATOM_COLORS = {
    "H": "#66c2a5",
    "N": "#fc8d62",
    "CA": "#8da0cb",
    "HA": "#e78ac3",
}

_VALID_ATOMS = frozenset(ATOM_ORDER)

# SI Fig. S15: N/H/Cα only (same CA gate as S11/S12; omit Hα so n matches CA cohort).
SF15_ATOM_ORDER: Tuple[str, ...] = ("H", "N", "CA")


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot F1 score distributions for 1D H/N/CA/HA CSPs on the subset of "
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
        default=DEFAULT_MIN_CA_SHIFT_ROW_COVERAGE,
        help=(
            "Minimum fraction of residues in 1d_analysis.csv with both CA_apo "
            "and CA_holo populated for a target to be considered as having "
            "analyzable CA shifts (default: %(default)s, i.e. strictly >50%% when 0.5)."
        ),
    )
    return parser.parse_args(list(argv))


def _allowed_basenames_from_targets_csv(
    targets_csv: Optional[Path],
    outputs_dir: Path,
) -> Optional[Dict[str, bool]]:
    """Map CSP_UBQ-style rows to canonical ``outputs/{HOLO}_{apo_bmrb}/`` basenames."""
    if targets_csv is None:
        return None
    path = targets_csv.resolve()
    if not path.exists():
        raise FileNotFoundError(f"Targets CSV not found: {path}")
    rows = load_target_rows(path)
    resolved = resolve_target_rows(rows, outputs_dir.resolve(), log_warnings=False)
    return {p.name: True for p in resolved}


def _results_by_atom(
    results_H: Sequence[TargetResult],
    results_N: Sequence[TargetResult],
    results_CA: Sequence[TargetResult],
    results_HA: Sequence[TargetResult],
) -> Dict[str, Sequence[TargetResult]]:
    return {"H": results_H, "N": results_N, "CA": results_CA, "HA": results_HA}


def _build_paired_arrays(
    results_by_atom: Dict[str, Sequence[TargetResult]],
    atom_order: Tuple[str, ...],
) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Intersect on target names; one F1 vector per atom in ``atom_order``."""
    by_target = {atom: {r.target: r.f1 for r in results_by_atom[atom]} for atom in atom_order}

    common: set[str] = set(by_target[atom_order[0]])
    for atom in atom_order[1:]:
        common &= set(by_target[atom])

    sorted_targets = sorted(common)
    arrays = {
        atom: np.array([by_target[atom][t] for t in sorted_targets], dtype=float)
        for atom in atom_order
    }
    return sorted_targets, arrays


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


def _paired_wilcoxon_table(atom_order: Tuple[str, ...], arrays: Dict[str, np.ndarray]) -> pd.DataFrame:
    """Compute paired Wilcoxon signed-rank for all unordered atom-type pairs."""
    pairs: List[Tuple[str, str, np.ndarray, np.ndarray]] = [
        (a, b, arrays[a], arrays[b]) for a, b in combinations(atom_order, 2)
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


def _atom_positions(atom_order: Tuple[str, ...]) -> Dict[str, int]:
    return {atom: i + 1 for i, atom in enumerate(atom_order)}


def _render_boxplot(
    atom_order: Tuple[str, ...],
    arrays: Dict[str, np.ndarray],
    stats_df: pd.DataFrame,
    output_image: Path,
    n_targets: int,
) -> None:
    """Draw the annotated boxplot and save to ``output_image``."""
    data_for_plot = [arrays[atom] for atom in atom_order]
    positions = list(_atom_positions(atom_order).values())

    fig_w = 6.0 + 0.7 * (len(atom_order) - 3)
    fig, ax = plt.subplots(figsize=(fig_w, 6.5))

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
    for patch, atom in zip(bp["boxes"], atom_order):
        patch.set_facecolor(ATOM_COLORS[atom])
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")

    rng = np.random.default_rng(seed=0)
    for pos, values, atom in zip(positions, data_for_plot, atom_order):
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
    ax.set_xticklabels(list(atom_order))
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

    _annotate_significance(ax, stats_df, atom_order)

    fig.tight_layout()
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _annotate_significance(
    ax: plt.Axes,
    stats_df: pd.DataFrame,
    atom_order: Tuple[str, ...],
) -> None:
    """Draw stacked significance bars over the boxes (within y in [0, 1])."""
    pos = _atom_positions(atom_order)
    base_y = 0.68
    step = 0.04
    bar_height = 0.012

    pairs_all = list(combinations(atom_order, 2))
    ordering = sorted(
        pairs_all,
        key=lambda ab: (abs(pos[ab[0]] - pos[ab[1]]), ab[0], ab[1]),
    )

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


def run_f1_1d_boxplot(
    outputs_dir: Path,
    *,
    targets_csv: Optional[Path],
    output_image: Path,
    stats_csv: Optional[Path],
    min_ca_coverage: float,
    atom_order: Tuple[str, ...] = ATOM_ORDER,
) -> int:
    """
    Pairwise Wilcoxon-annotated boxplot of per-target F1 scores for ``atom_order`` 1D CSPs.

    When ``targets_csv`` is given, rows are mapped to canonical ``outputs/`` subdirectory
    names via :mod:`scripts.target_resolution`.
    """
    if not atom_order:
        print("atom_order must be non-empty.", file=sys.stderr)
        return 1
    if not frozenset(atom_order) <= _VALID_ATOMS or len(atom_order) != len(set(atom_order)):
        print(
            f"Invalid atom_order={atom_order!r}; "
            "use a tuple of distinct keys from {'H','N','CA','HA'}.",
            file=sys.stderr,
        )
        return 1

    outputs_dir = outputs_dir.resolve()
    if not outputs_dir.exists():
        print(f"Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1

    try:
        allowed_targets = _allowed_basenames_from_targets_csv(targets_csv, outputs_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    ca_passing, ca_coverages = target_basenames_passing_ca_shift_coverage(
        outputs_dir,
        min_coverage=min_ca_coverage,
        allowed_basenames=allowed_targets,
    )
    thr = float(min_ca_coverage)
    print(
        f"{len(ca_passing)} targets pass CA coverage > {thr:.0%} "
        f"(of {len(ca_coverages)} candidates inspected)"
    )
    if not ca_passing:
        print(
            "No targets satisfy the CA coverage threshold; nothing to plot.",
            file=sys.stderr,
        )
        return 1

    ca_allowed = {t: True for t in ca_passing}
    results_H, results_N, results_CA, results_HA = collect_1d_f1_results(outputs_dir, ca_allowed)
    results_by = _results_by_atom(results_H, results_N, results_CA, results_HA)
    for atom in atom_order:
        if not results_by[atom]:
            print(
                f"No 1D F1 results for atom {atom!r}; cannot render boxplot.",
                file=sys.stderr,
            )
            return 1

    common, arrays = _build_paired_arrays(results_by, atom_order)
    if len(common) == 0:
        labels = "/".join(atom_order)
        print(f"No targets produce F1 scores for all of {labels}; nothing to plot.", file=sys.stderr)
        return 1

    stats_df = _paired_wilcoxon_table(atom_order, arrays)

    resolved_stats = stats_csv.resolve() if stats_csv is not None else _resolve_stats_csv(
        output_image, None
    )

    output_image.parent.mkdir(parents=True, exist_ok=True)
    _render_boxplot(atom_order, arrays, stats_df, output_image, n_targets=len(common))

    resolved_stats.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(resolved_stats, index=False)

    joined = "/".join(atom_order)
    print(f"n = {len(common)} targets with paired {joined} F1 scores")
    print(f"Wrote {output_image}")
    print(f"Wrote {resolved_stats}")
    for _, row in stats_df.iterrows():
        print(
            f"  {row['group_a']} vs {row['group_b']}: "
            f"W={row['statistic']:.3g}, p_raw={row['p_raw']:.3g}, "
            f"p_adj(Holm)={row['p_adj_holm']:.3g} ({row['signif']})"
        )

    return 0


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    output_image = args.output_image
    tc = args.targets_csv.resolve() if args.targets_csv else None

    return run_f1_1d_boxplot(
        args.outputs_dir.resolve(),
        targets_csv=tc,
        output_image=output_image,
        stats_csv=args.stats_csv,
        min_ca_coverage=float(args.min_ca_coverage),
        atom_order=ATOM_ORDER,
    )


if __name__ == "__main__":
    raise SystemExit(main())
