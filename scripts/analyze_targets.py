#!/usr/bin/env python3
"""
Analyze master alignment outputs and visualize per-target F1 scores.

When invoked, this script now also orchestrates additional analyses:
- N/H CSP F1 summary from `master_alignment.csv` (this module).
- CA-inclusive CSP F1 summary from `csp_table_CA.csv` (via analyze_targets_ca).
- CA-inclusive CSP F1 summary restricted to same-author pairs (via analyze_targets_same_author).
- Optional single-atom (H, N, CA) 1D shift analysis helper (via analyze_targets_single_atom_shifts).

All top-level summary artifacts are written into `outputs/summary_statistics`
by default (configurable via --summary-dir).
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
import math
from math import ceil
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

try:
    from .config import classification_colors
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors


SIGNIFICANT_COLUMN = "significant"
CA_DISTANCE_COLUMN = "min_ca_distance_distance"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


class AlignmentParsingError(RuntimeError):
    """Raised when a `master_alignment.csv` file cannot be parsed."""


@dataclass
class TargetResult:
    target: str
    f1: float
    mcc: float
    true_positives: int
    false_positives: int
    false_negatives: int
    total_rows: int


@dataclass
class DistanceRecord:
    distance: float
    is_predicted_positive: bool


@dataclass
class ConfusionRecord:
    """Record for confusion matrix: ground truth = binding site, prediction = CSP significance."""
    distance: float
    is_binding: bool  # ground truth: residue in binding site (predictor columns)
    is_significant: bool  # prediction: CSP significant


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-target F1 scores from master_alignment.csv files "
            "and render them as a heatmap."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--output-image",
        type=Path,
        default=Path("outputs") / "f1_heatmap.png",
        help="Destination for the generated heatmap image (default: %(default)s).",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write the per-target F1 summary table as CSV.",
    )
    parser.add_argument(
        "--histogram-image",
        type=Path,
        default=Path("outputs") / "significant_ca_distance_hist.png",
        help=(
            "Destination for the histogram of significant residues by CA distance "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--stacked-histogram-image",
        type=Path,
        default=Path("outputs") / "significant_ca_distance_stacked_hist.png",
        help=(
            "Destination for the stacked histogram splitting predicted positives "
            "and negatives among significant residues (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--confusion-matrix-stacked-histogram-image",
        type=Path,
        default=Path("outputs") / "confusion_matrix_stacked_histogram.png",
        help="Destination for confusion matrix stacked histogram.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        help=(
            "Optional CSV file containing target IDs to filter. "
            "Should have a 'holo_pdb' column. Only targets listed in this file will be included."
        ),
    )
    parser.add_argument(
        "--targets",
        type=str,
        metavar="IDS",
        help=(
            "Comma-separated list of holo_pdb IDs to include (e.g. 1cf4,1d5g,1l8c). "
            "If this set matches a previously run selection, that name is reused. Otherwise prompts for a name."
        ),
    )
    parser.add_argument(
        "--targets-from",
        type=str,
        metavar="NAME",
        help="Use targets from an existing selection by name (e.g. CBP, TFIIH). Looks for outputs/targets_{NAME}.csv.",
    )
    parser.add_argument(
        "--selection-name",
        type=str,
        help="Name for a new target selection (used with --targets). If omitted and no matching selection exists, prompts interactively.",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=Path("outputs") / "summary_statistics",
        help="Directory to store all summary statistics outputs (default: %(default)s).",
    )
    parser.add_argument(
        "--no-plot-titles",
        action="store_true",
        help="Hide titles on stacked histogram plots.",
    )
    parser.add_argument(
        "--axis-fontsize",
        type=float,
        default=None,
        help="Font size for axis labels on stacked histogram plots.",
    )
    parser.add_argument(
        "--tick-fontsize",
        type=float,
        default=None,
        help="Font size for axis tick labels on stacked histogram plots.",
    )
    parser.add_argument(
        "--legend-fontsize",
        type=float,
        default=None,
        help="Font size for legend text on stacked histogram plots.",
    )
    parser.add_argument(
        "--legend-title-fontsize",
        type=float,
        default=None,
        help="Font size for legend title on confusion stacked histogram plots.",
    )
    parser.add_argument(
        "--compact-legend-labels",
        action="store_true",
        help="Use compact legend labels (e.g., TP (N), FP (N), ...).",
    )
    return parser.parse_args(list(argv))


def _run_ca_analysis(summary_root: Path, outputs_dir: Path, targets_csv: Optional[Path]) -> None:
    """Invoke analyze_targets_ca.main with outputs routed into summary_root."""
    if targets_csv is None:
        print("[ORCH] Skipping CA-inclusive F1 analysis (analyze_targets_ca): --targets-csv not provided")
        return
    try:
        try:
            from . import analyze_targets_ca as at_ca  # type: ignore
        except Exception:
            import analyze_targets_ca as at_ca  # type: ignore
        ca_args = [
            "--outputs-dir",
            str(outputs_dir),
            "--targets-csv",
            str(targets_csv),
            "--output-image",
            str(summary_root / "f1_heatmap_ca.png"),
            "--histogram-image",
            str(summary_root / "significant_ca_distance_hist_ca.png"),
            "--stacked-histogram-image",
            str(summary_root / "significant_ca_distance_stacked_hist_ca.png"),
            "--scatterplot-image",
            str(summary_root / "f1_comparison_scatterplot_ca_vs_nh.png"),
            "--summary-csv",
            str(summary_root / "f1_summary_ca.csv"),
        ]
        print("[ORCH] Running CA-inclusive F1 analysis (analyze_targets_ca)")
        at_ca.main(ca_args)
    except Exception as exc:  # pragma: no cover - orchestration best-effort
        print(f"[ORCH] WARNING: analyze_targets_ca failed: {exc}", file=sys.stderr)


def _run_same_author_analysis(summary_root: Path, outputs_dir: Path) -> None:
    """Invoke analyze_targets_same_author.main with outputs routed into summary_root."""
    try:
        try:
            from . import analyze_targets_same_author as at_same  # type: ignore
        except Exception:
            import analyze_targets_same_author as at_same  # type: ignore
        same_args = [
            "--outputs-dir",
            str(outputs_dir),
            "--csprank-csv",
            "CSP_UBQ.csv",
            "--output-image",
            str(summary_root / "f1_heatmap_same_author.png"),
            "--histogram-image",
            str(summary_root / "significant_ca_distance_hist_same_author.png"),
            "--stacked-histogram-image",
            str(summary_root / "significant_ca_distance_stacked_hist_same_author.png"),
            "--scatterplot-image",
            str(summary_root / "f1_comparison_scatterplot_same_author_ca_vs_nh.png"),
            "--summary-csv",
            str(summary_root / "f1_summary_same_author.csv"),
        ]
        print("[ORCH] Running same-author CA-inclusive F1 analysis (analyze_targets_same_author)")
        at_same.main(same_args)
    except Exception as exc:  # pragma: no cover
        print(f"[ORCH] WARNING: analyze_targets_same_author failed: {exc}", file=sys.stderr)


def _run_single_atom_analysis(summary_root: Path, outputs_dir: Path, targets_csv: Optional[Path]) -> None:
    """Invoke analyze_targets_single_atom_shifts.main if available."""
    try:
        try:
            from . import analyze_targets_single_atom_shifts as at_1d  # type: ignore
        except Exception:
            import analyze_targets_single_atom_shifts as at_1d  # type: ignore
        one_d_args = [
            "--outputs-dir",
            str(outputs_dir),
            "--output-image",
            str(summary_root / "f1_heatmap_1d_single_atom.png"),
            "--summary-csv",
            str(summary_root / "f1_summary_1d_single_atom.csv"),
        ]
        if targets_csv is not None:
            one_d_args.extend(["--targets-csv", str(targets_csv)])
        print("[ORCH] Running single-atom 1D analysis helper (analyze_targets_single_atom_shifts)")
        at_1d.main(one_d_args)
    except Exception as exc:  # pragma: no cover
        print(f"[ORCH] WARNING: analyze_targets_single_atom_shifts failed: {exc}", file=sys.stderr)


def load_targets_from_csv(csv_path: Path) -> Set[str]:
    """
    Load allowed target IDs from a CSV file.
    
    Args:
        csv_path: Path to CSV file with 'holo_pdb' column
        
    Returns:
        Set of target IDs (holo_pdb values)
    """
    try:
        df = pd.read_csv(csv_path)
        if 'holo_pdb' not in df.columns:
            raise ValueError(f"CSV file {csv_path} must have a 'holo_pdb' column")
        # Convert to set and strip whitespace
        targets = set(df['holo_pdb'].astype(str).str.strip())
        return targets
    except Exception as exc:
        raise ValueError(f"Failed to load targets from {csv_path}: {exc}") from exc


def to_bool(value) -> bool:
    """Coerce assorted representations into booleans."""
    if pd.isna(value):
        return False
    if isinstance(value, (bool, int)):
        return bool(value)
    if isinstance(value, float):
        return bool(round(value))
    if isinstance(value, str):
        text = value.strip().lower()
        if not text:
            return False
        if text in {"true", "t", "yes", "y", "1"}:
            return True
        if text in {"false", "f", "no", "n", "0"}:
            return False
    raise AlignmentParsingError(f"Unable to interpret value {value!r} as boolean.")


def load_alignment(alignment_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(alignment_path)
    except Exception as exc:
        raise AlignmentParsingError(
            f"Failed to read alignment CSV: {alignment_path}"
        ) from exc

    required_columns = (SIGNIFICANT_COLUMN, CA_DISTANCE_COLUMN, *PREDICTOR_COLUMNS)
    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]
    if missing_columns:
        raise AlignmentParsingError(
            f"Alignment file {alignment_path} is missing required columns: "
            f"{', '.join(missing_columns)}"
        )

    for column in (SIGNIFICANT_COLUMN, *PREDICTOR_COLUMNS):
        df[column] = df[column].apply(to_bool)

    df[CA_DISTANCE_COLUMN] = pd.to_numeric(df[CA_DISTANCE_COLUMN], errors="coerce")

    return df


def compute_f1_score(df: pd.DataFrame, predicted: pd.Series | None = None) -> TargetResult:
    actual = df[SIGNIFICANT_COLUMN]
    if predicted is None:
        predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)

    true_positives = int((actual & predicted).sum())
    false_positives = int((~actual & predicted).sum())
    false_negatives = int((actual & ~predicted).sum())
    total_rows = int(len(df))
    true_negatives = total_rows - true_positives - false_positives - false_negatives

    denominator = 2 * true_positives + false_positives + false_negatives
    f1 = (2 * true_positives / denominator) if denominator else 0.0

    # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    mcc_denom = (
        (true_positives + false_positives)
        * (true_positives + false_negatives)
        * (true_negatives + false_positives)
        * (true_negatives + false_negatives)
    )
    if mcc_denom <= 0:
        mcc = 0.0
    else:
        mcc = (true_positives * true_negatives - false_positives * false_negatives) / math.sqrt(
            mcc_denom
        )
        mcc = max(-1.0, min(1.0, mcc))

    return TargetResult(
        target="",
        f1=f1,
        mcc=mcc,
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        total_rows=total_rows,
    )


def discover_targets(outputs_dir: Path) -> List[Path]:
    if not outputs_dir.exists():
        raise FileNotFoundError(f"Outputs directory not found: {outputs_dir}")

    return sorted(
        path for path in outputs_dir.iterdir() if path.is_dir() and not path.name.startswith(".")
    )


def collect_results(
    outputs_dir: Path, 
    allowed_targets: Optional[Set[str]] = None
) -> Tuple[List[TargetResult], List[DistanceRecord], List[ConfusionRecord]]:
    results: List[TargetResult] = []
    distances: List[DistanceRecord] = []
    confusion_records: List[ConfusionRecord] = []

    for target_dir in discover_targets(outputs_dir):
        # Filter by allowed targets if provided
        if allowed_targets is not None and target_dir.name not in allowed_targets:
            continue
            
        alignment_path = target_dir / "master_alignment.csv"
        if not alignment_path.exists():
            continue

        try:
            df = load_alignment(alignment_path)
        except AlignmentParsingError as exc:
            print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
            continue

        predicted = df[list(PREDICTOR_COLUMNS)].any(axis=1)
        actual = df[SIGNIFICANT_COLUMN]
        metrics = compute_f1_score(df, predicted)
        metrics.target = target_dir.name
        results.append(metrics)

        significant_mask = df[SIGNIFICANT_COLUMN] & df[CA_DISTANCE_COLUMN].notna()
        if significant_mask.any():
            selected_distances = df.loc[significant_mask, CA_DISTANCE_COLUMN]
            selected_predictions = predicted.loc[significant_mask]
            for distance_value, predicted_value in zip(
                selected_distances.tolist(), selected_predictions.tolist()
            ):
                distances.append(
                    DistanceRecord(
                        distance=float(distance_value),
                        is_predicted_positive=bool(predicted_value),
                    )
                )

        # Collect all residues with CA distances for confusion matrix histogram.
        # Uses same definitions as plot_csp_classification_bars: ground truth = binding site,
        # prediction = CSP significance. TP=Sig. CSP in binding site, FP=Sig. CSP allosteric,
        # FN=low CSP in binding site, TN=low CSP allosteric.
        ca_distance_mask = df[CA_DISTANCE_COLUMN].notna()
        if ca_distance_mask.any():
            subset = df.loc[ca_distance_mask]
            for idx, row in subset.iterrows():
                confusion_records.append(
                    ConfusionRecord(
                        distance=float(row[CA_DISTANCE_COLUMN]),
                        is_binding=bool(predicted.loc[idx]),  # ground truth: in binding site
                        is_significant=bool(row[SIGNIFICANT_COLUMN]),  # prediction: CSP significant
                    )
                )

    return results, distances, confusion_records


def render_heatmap(results: Sequence[TargetResult], output_image: Path) -> None:
    dataframe = pd.DataFrame([result.__dict__ for result in results]).set_index("target")
    dataframe = dataframe.sort_values("f1", ascending=False)
    heatmap_data = dataframe[["f1", "mcc"]].rename(columns={"f1": "F1", "mcc": "MCC"})

    plt.figure(figsize=(6, max(3, len(heatmap_data) * 0.3)))
    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        vmin=-1.0,
        vmax=1.0,
        cbar_kws={"label": "F1 / MCC"},
    )
    # set larger font sizes for labels and title
    ax.set_xlabel("Metric", fontsize=16)
    ax.set_ylabel("Target", fontsize=16)
    ax.set_title("Per-target F1 and MCC", fontsize=18)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label("F1 / MCC", fontsize=16)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_histogram(distance_records: Sequence[DistanceRecord], output_image: Path) -> None:
    distances = [record.distance for record in distance_records]
    if not distances:
        print("No significant residues with CA distances available; skipping histogram.", file=sys.stderr)
        return

    max_distance = max(distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=bins, edgecolor="black", color="#4c72b0")
    plt.xlabel("Minimum CA Distance (Å)")
    plt.ylabel("Number of Significant Residues")
    plt.title("Distribution of Significant Residues by Minimum CA Distance")
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_stacked_histogram(
    distance_records: Sequence[DistanceRecord],
    output_image: Path,
    positive_count: int,
    negative_count: int,
    *,
    show_title: bool = True,
    axis_fontsize: Optional[float] = None,
    tick_fontsize: Optional[float] = None,
    legend_fontsize: Optional[float] = None,
    compact_legend_labels: bool = False,
) -> None:
    if not distance_records:
        print(
            "No significant residues with CA distances available; skipping stacked histogram.",
            file=sys.stderr,
        )
        return

    positive_distances = [record.distance for record in distance_records if record.is_predicted_positive]
    negative_distances = [record.distance for record in distance_records if not record.is_predicted_positive]

    if not positive_distances and not negative_distances:
        print(
            "No significant residues with CA distances available; skipping stacked histogram.",
            file=sys.stderr,
        )
        return

    max_distance = max(positive_distances + negative_distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    if compact_legend_labels:
        labels = [f"TP ({positive_count})", f"FP ({negative_count})"]
    else:
        labels = [
            f"(TP) Sig. CSP in Binding Site ({positive_count})",
            f"(FP) Sig. CSP -- Allosteric ({negative_count})",
        ]

    ax.hist(
        [positive_distances, negative_distances],
        bins=bins,
        stacked=True,
        color=[classification_colors.TP, classification_colors.FP],
        edgecolor="black",
        label=labels,
    )
    plt.xlabel("Minimum CA Distance (Å)", fontsize=axis_fontsize)
    plt.ylabel("Number of Residues", fontsize=axis_fontsize)
    # if show_title:
    #     plt.title("Predicted Outcomes for Significant Residues by Minimum CA Distance")
    if tick_fontsize is not None:
        ax.tick_params(axis="both", labelsize=tick_fontsize)
    legend_kwargs = {}
    if legend_fontsize is not None:
        legend_kwargs["fontsize"] = legend_fontsize
    plt.legend(**legend_kwargs)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def render_confusion_matrix_stacked_histogram(
    confusion_records: Sequence[ConfusionRecord],
    output_image: Path,
    *,
    show_title: bool = True,
    axis_fontsize: Optional[float] = None,
    tick_fontsize: Optional[float] = None,
    legend_fontsize: Optional[float] = None,
    legend_title_fontsize: Optional[float] = None,
    compact_legend_labels: bool = False,
) -> None:
    """Render stacked histogram of TP, FP, FN, TN by Minimum CA Distance.

    Uses same definitions as plot_csp_classification_bars: ground truth = binding site,
    prediction = CSP significance.
    """
    if not confusion_records:
        print(
            "No residues with CA distances available; skipping confusion matrix stacked histogram.",
            file=sys.stderr,
        )
        return

    # Same classification as plot_csp_classification_bars: TP=sig+binding, FP=sig+allosteric,
    # FN=low+binding, TN=low+allosteric
    tp_distances = [
        r.distance for r in confusion_records
        if r.is_significant and r.is_binding
    ]
    fp_distances = [
        r.distance for r in confusion_records
        if r.is_significant and not r.is_binding
    ]
    fn_distances = [
        r.distance for r in confusion_records
        if not r.is_significant and r.is_binding
    ]
    tn_distances = [
        r.distance for r in confusion_records
        if not r.is_significant and not r.is_binding
    ]

    all_distances = tp_distances + fp_distances + fn_distances + tn_distances
    max_distance = max(all_distances)
    upper_edge = ceil(max_distance) + 1
    bins = list(range(0, upper_edge + 1))

    # Colors from config: TN, FP, FN, TP. Stack order bottom to top: TN, FP, FN, TP
    colors = [
        classification_colors.TN,
        classification_colors.FP,
        classification_colors.FN,
        classification_colors.TP,
    ]
    tp_count = len(tp_distances)
    fp_count = len(fp_distances)
    fn_count = len(fn_distances)
    tn_count = len(tn_distances)
    if compact_legend_labels:
        labels = [
            f"TN ({tn_count})",
            f"FP ({fp_count})",
            f"FN ({fn_count})",
            f"TP ({tp_count})",
        ]
    else:
        labels = [
            f"(TN) Small CSP -- allosteric ({tn_count})",
            f"(FP) Sig. CSP -- allosteric ({fp_count})",
            f"(FN) Small CSP in Binding Site ({fn_count})",
            f"(TP) Sig. CSP in Binding Site ({tp_count})",
        ]

    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    ax.hist(
        [tn_distances, fp_distances, fn_distances, tp_distances],
        bins=bins,
        stacked=True,
        color=colors,
        edgecolor="black",
        label=labels,
    )
    plt.xlabel("Minimum CA Distance (Å)", fontsize=axis_fontsize)
    plt.ylabel("Number of Residues", fontsize=axis_fontsize)
    # if show_title:
    #     plt.title("Confusion Matrix by Minimum CA Distance")
    if tick_fontsize is not None:
        ax.tick_params(axis="both", labelsize=tick_fontsize)
    legend_kwargs = {"title": "Confusion Matrix"}
    if legend_fontsize is not None:
        legend_kwargs["fontsize"] = legend_fontsize
    if legend_title_fontsize is not None:
        legend_kwargs["title_fontsize"] = legend_title_fontsize
    plt.legend(**legend_kwargs)
    plt.tight_layout()

    output_image.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_image, dpi=300)
    plt.close()


def _parse_targets_list(targets_str: str) -> Set[str]:
    """Parse comma-separated holo_pdb IDs into a set."""
    ids = [x.strip() for x in targets_str.split(",") if x.strip()]
    return set(ids)


def _discover_existing_selections(outputs_dir: Path) -> Dict[str, Set[str]]:
    """Load existing selections from targets_*.csv files in outputs_dir."""
    selections: Dict[str, Set[str]] = {}
    for path in outputs_dir.glob("targets_*.csv"):
        name = path.stem.replace("targets_", "", 1)
        try:
            df = pd.read_csv(path)
            if "holo_pdb" in df.columns:
                targets = set(df["holo_pdb"].astype(str).str.strip())
                selections[name] = targets
        except Exception:
            continue
    return selections


def _find_matching_selection(
    target_set: Set[str], existing: Dict[str, Set[str]]
) -> Optional[str]:
    """Return the name of an existing selection with the same target set, or None."""
    for name, targets in existing.items():
        if targets == target_set:
            return name
    return None


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir.resolve()
    output_image = args.output_image.resolve()
    summary_csv = args.summary_csv.resolve() if args.summary_csv else None
    histogram_image = args.histogram_image.resolve()
    stacked_histogram_image = args.stacked_histogram_image.resolve()
    confusion_matrix_stacked_histogram_image = (
        args.confusion_matrix_stacked_histogram_image.resolve()
    )
    summary_root = args.summary_dir.resolve()
    summary_root.mkdir(parents=True, exist_ok=True)

    # Load allowed targets: --targets-from > --targets > --targets-csv
    allowed_targets = None
    targets_csv_path = args.targets_csv.resolve() if args.targets_csv else None
    selection_name: Optional[str] = None

    if args.targets_from:
        # Load targets from existing selection by name
        targets_file = outputs_dir / f"targets_{args.targets_from}.csv"
        if not targets_file.exists():
            print(
                f"Selection '{args.targets_from}' not found: {targets_file} does not exist.",
                file=sys.stderr,
            )
            return 1
        allowed_targets = load_targets_from_csv(targets_file)
        selection_name = args.targets_from
        print(f"Using selection '{selection_name}' ({len(allowed_targets)} targets)")
        outputs_base = outputs_dir
        output_image = outputs_base / f"f1_heatmap_{selection_name}.png"
        histogram_image = outputs_base / f"significant_ca_distance_hist_{selection_name}.png"
        stacked_histogram_image = outputs_base / f"significant_ca_distance_stacked_hist_{selection_name}.png"
        confusion_matrix_stacked_histogram_image = (
            outputs_base / f"confusion_matrix_stacked_histogram_{selection_name}.png"
        )
        summary_csv = outputs_base / f"f1_summary_{selection_name}.csv"
        summary_root = outputs_base / f"summary_statistics_{selection_name}"
        summary_root.mkdir(parents=True, exist_ok=True)
        targets_csv_path = targets_file
    elif args.targets:
        allowed_targets = _parse_targets_list(args.targets)
        existing = _discover_existing_selections(outputs_dir)
        match = _find_matching_selection(allowed_targets, existing)
        if match is not None:
            selection_name = match
            print(f"Target set matches existing selection '{selection_name}' ({len(allowed_targets)} targets)")
        else:
            selection_name = args.selection_name
            if not selection_name:
                selection_name = input("Enter a name for this target selection: ").strip()
            if not selection_name:
                print("No name provided; using 'selection'", file=sys.stderr)
                selection_name = "selection"
            print(f"Filtering to {len(allowed_targets)} targets from --targets: {sorted(allowed_targets)}")
        outputs_base = outputs_dir
        output_image = outputs_base / f"f1_heatmap_{selection_name}.png"
        histogram_image = outputs_base / f"significant_ca_distance_hist_{selection_name}.png"
        stacked_histogram_image = outputs_base / f"significant_ca_distance_stacked_hist_{selection_name}.png"
        confusion_matrix_stacked_histogram_image = (
            outputs_base / f"confusion_matrix_stacked_histogram_{selection_name}.png"
        )
        summary_csv = outputs_base / f"f1_summary_{selection_name}.csv"
        summary_root = outputs_base / f"summary_statistics_{selection_name}"
        summary_root.mkdir(parents=True, exist_ok=True)
        targets_csv_path = outputs_base / f"targets_{selection_name}.csv"
        pd.DataFrame({"holo_pdb": sorted(allowed_targets)}).to_csv(
            targets_csv_path, index=False
        )
    elif args.targets_csv:
        allowed_targets = load_targets_from_csv(args.targets_csv)
        print(f"Filtering to {len(allowed_targets)} targets from {args.targets_csv}")

    results, distances, confusion_records = collect_results(outputs_dir, allowed_targets)
    if not results:
        print(
            f"No master_alignment.csv files found beneath {outputs_dir}",
            file=sys.stderr,
        )
        return 1

    if summary_csv:
        summary_df = pd.DataFrame([result.__dict__ for result in results]).sort_values("target")
        summary_csv.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_csv, index=False)

    positive_count = sum(1 for record in distances if record.is_predicted_positive)
    negative_count = sum(1 for record in distances if not record.is_predicted_positive)

    print(f"Significant CSPs predicted significant: {positive_count}")
    print(f"Significant CSPs missed: {negative_count}")

    # When using a target list, report aggregate F1 and MCC from accumulated confusion matrix
    if allowed_targets is not None:
        agg_tp = sum(r.true_positives for r in results)
        agg_fp = sum(r.false_positives for r in results)
        agg_fn = sum(r.false_negatives for r in results)
        agg_tn = sum(r.total_rows - r.true_positives - r.false_positives - r.false_negatives for r in results)
        denom = 2 * agg_tp + agg_fp + agg_fn
        agg_f1 = (2 * agg_tp / denom) if denom else 0.0
        agg_mcc_denom = (
            (agg_tp + agg_fp) * (agg_tp + agg_fn) * (agg_tn + agg_fp) * (agg_tn + agg_fn)
        )
        if agg_mcc_denom <= 0:
            agg_mcc = 0.0
        else:
            agg_mcc = (agg_tp * agg_tn - agg_fp * agg_fn) / math.sqrt(agg_mcc_denom)
            agg_mcc = max(-1.0, min(1.0, agg_mcc))
        print(f"\nAggregate F1 (accumulated across {len(results)} targets): {agg_f1:.4f}")
        print(f"Aggregate MCC: {agg_mcc:.4f}")
        print(f"  TP={agg_tp}, FP={agg_fp}, FN={agg_fn}, TN={agg_tn}")

    render_heatmap(results, output_image)
    render_histogram(distances, histogram_image)
    render_stacked_histogram(
        distances,
        stacked_histogram_image,
        positive_count,
        negative_count,
        show_title=not args.no_plot_titles,
        axis_fontsize=args.axis_fontsize,
        tick_fontsize=args.tick_fontsize,
        legend_fontsize=args.legend_fontsize,
        compact_legend_labels=args.compact_legend_labels,
    )
    render_confusion_matrix_stacked_histogram(
        confusion_records,
        confusion_matrix_stacked_histogram_image,
        show_title=not args.no_plot_titles,
        axis_fontsize=args.axis_fontsize,
        tick_fontsize=args.tick_fontsize,
        legend_fontsize=args.legend_fontsize,
        legend_title_fontsize=args.legend_title_fontsize,
        compact_legend_labels=args.compact_legend_labels,
    )
    print(f"Heatmap saved to {output_image}")
    print(f"Histogram saved to {histogram_image}")
    print(f"Stacked histogram saved to {stacked_histogram_image}")
    print(f"Confusion matrix stacked histogram saved to {confusion_matrix_stacked_histogram_image}")
    if summary_csv:
        print(f"Summary table written to {summary_csv}")

    # Canonical NH summary artifacts in summary_root
    nh_heatmap = summary_root / "f1_heatmap_nh.png"
    nh_hist = summary_root / "significant_ca_distance_hist_nh.png"
    nh_stacked = summary_root / "significant_ca_distance_stacked_hist_nh.png"
    nh_confusion_stacked = summary_root / "confusion_matrix_stacked_histogram_nh.png"
    nh_summary = summary_root / "f1_summary_nh.csv"

    print(f"[ORCH] Writing NH F1 summaries into {summary_root}")
    render_heatmap(results, nh_heatmap)
    render_histogram(distances, nh_hist)
    render_stacked_histogram(
        distances,
        nh_stacked,
        positive_count,
        negative_count,
        show_title=not args.no_plot_titles,
        axis_fontsize=args.axis_fontsize,
        tick_fontsize=args.tick_fontsize,
        legend_fontsize=args.legend_fontsize,
        compact_legend_labels=args.compact_legend_labels,
    )
    render_confusion_matrix_stacked_histogram(
        confusion_records,
        nh_confusion_stacked,
        show_title=not args.no_plot_titles,
        axis_fontsize=args.axis_fontsize,
        tick_fontsize=args.tick_fontsize,
        legend_fontsize=args.legend_fontsize,
        legend_title_fontsize=args.legend_title_fontsize,
        compact_legend_labels=args.compact_legend_labels,
    )
    summary_df = pd.DataFrame([result.__dict__ for result in results]).sort_values("target")
    nh_summary.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(nh_summary, index=False)

    # Run additional analyses and route their outputs into summary_root
    # targets_csv_path was set above when using --targets-from or --targets
    if targets_csv_path is None and args.targets_csv:
        targets_csv_path = args.targets_csv.resolve()
    _run_ca_analysis(summary_root, outputs_dir, targets_csv_path)
    _run_same_author_analysis(summary_root, outputs_dir)
    _run_single_atom_analysis(summary_root, outputs_dir, targets_csv_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

