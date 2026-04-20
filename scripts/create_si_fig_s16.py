#!/usr/bin/env python3
"""
Generator script filename (s16) ≠ SI index: this produces SI Fig. S13.

SI Fig. S13 — CSP DB confusion matrix histograms by closest interchain atom–atom distance.

Two-panel figure:
- Panel A: stacked histogram of significant residues (TP vs FP only)
- Panel B: stacked confusion-matrix histogram (TN, FP, FN, TP)

Output: figures/SF13_any_atom_distance.png

Default targets list: CSP_UBQ_ph0.5_temp5C.csv (buffer-filtered subset). Override with --targets-csv.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

try:
    from .config import classification_colors, paths
    from .interaction_analysis import compute_min_atom_distance_filter
    from .rcsb_io import fetch_pdb
except Exception:
    import os as _os
    import sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import classification_colors, paths
    from scripts.interaction_analysis import compute_min_atom_distance_filter
    from scripts.rcsb_io import fetch_pdb


SIGNIFICANT_COLUMN = "significant"
DISTANCE_COLUMN = "min_any_atom_distance_any_atom"
PREDICTOR_COLUMNS: Sequence[str] = (
    "passes_filter_distance",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "has_hbond_interaction",
    "is_occluded_occlusion",
)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _parse_targets(targets_csv: Optional[Path], targets_str: Optional[str]) -> Optional[Set[str]]:
    selected: Set[str] = set()
    if targets_csv:
        df = pd.read_csv(targets_csv)
        if "holo_pdb" not in df.columns:
            raise ValueError(f"Missing 'holo_pdb' column in {targets_csv}")
        selected.update(df["holo_pdb"].astype(str).str.strip().str.lower().tolist())
    if targets_str:
        selected.update(t.strip().lower() for t in targets_str.split(",") if t.strip())
    return selected or None


def _target_is_selected(target_name: str, selected: Optional[Set[str]]) -> bool:
    if selected is None:
        return True
    target_lower = target_name.lower()
    base_target = target_lower.split("_")[0]
    return target_lower in selected or base_target in selected


def _resolve_pdb_path(holo_pdb: str) -> Optional[str]:
    import os
    holo_pdb = str(holo_pdb).strip().lower()
    if not holo_pdb:
        return None
    local_path = os.path.join(paths.pdb_cache_dir, f"{holo_pdb}.pdb")
    if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
        return local_path
    try:
        return fetch_pdb(holo_pdb, cache_dir=paths.pdb_cache_dir)
    except Exception:
        return None


def _compute_distances_for_target(
    alignment_path: Path,
    distance_key: str,
) -> Optional[dict]:
    df = pd.read_csv(alignment_path)
    if "holo_pdb" not in df.columns or "pdb_residue_number" not in df.columns:
        return None
    holo_pdb = df["holo_pdb"].iloc[0]
    if pd.isna(holo_pdb) or not str(holo_pdb).strip():
        return None
    holo_pdb = str(holo_pdb).strip().lower()
    pdb_path = _resolve_pdb_path(holo_pdb)
    if not pdb_path:
        return None
    result = compute_min_atom_distance_filter(pdb_path, distance_threshold=6.0)
    if "error" in result or not result.get("residue_info"):
        return None
    return {info["residue_number"]: info[distance_key] for info in result["residue_info"]}


def _get_bins(data: Sequence[float], bin_width: float, max_distance: Optional[float]) -> tuple[List[float], float]:
    if not data:
        raise ValueError("No inter-atomic distance data found for selected targets.")
    data_max = max(data) if max_distance is None else max_distance
    n_bins = max(1, int((data_max + bin_width) / bin_width))
    bins = [i * bin_width for i in range(n_bins + 1)]
    return bins, data_max


def collect_distance_categories(
    outputs_dir: Path,
    selected_targets: Optional[Set[str]] = None,
) -> tuple[List[float], List[float], List[float], List[float]]:
    tp_distances: List[float] = []
    fp_distances: List[float] = []
    fn_distances: List[float] = []
    tn_distances: List[float] = []

    base_required = {SIGNIFICANT_COLUMN, *PREDICTOR_COLUMNS, "pdb_residue_number", "holo_pdb"}

    for alignment_path in sorted(outputs_dir.glob("*/master_alignment.csv")):
        target_name = alignment_path.parent.name
        if not _target_is_selected(target_name, selected_targets):
            continue

        df = pd.read_csv(alignment_path)
        missing_base = [c for c in base_required if c not in df.columns]
        if missing_base:
            continue

        is_significant = df[SIGNIFICANT_COLUMN].map(_as_bool)
        is_binding = pd.DataFrame({c: df[c].map(_as_bool) for c in PREDICTOR_COLUMNS}).any(axis=1)

        if DISTANCE_COLUMN in df.columns:
            distances = pd.to_numeric(df[DISTANCE_COLUMN], errors="coerce")
        else:
            distance_map = _compute_distances_for_target(alignment_path, "min_any_atom_distance")
            if distance_map is None:
                continue
            pdb_resi = pd.to_numeric(df["pdb_residue_number"], errors="coerce")
            distances = pdb_resi.map(lambda x: distance_map.get(int(x), float("nan")) if pd.notna(x) else float("nan"))

        valid = distances.notna()
        if not valid.any():
            continue

        sig = is_significant[valid]
        bind = is_binding[valid]
        dist = distances[valid].astype(float)

        tp_distances.extend(dist[sig & bind].tolist())
        fp_distances.extend(dist[sig & ~bind].tolist())
        fn_distances.extend(dist[~sig & bind].tolist())
        tn_distances.extend(dist[~sig & ~bind].tolist())

    return tp_distances, fp_distances, fn_distances, tn_distances


def _set_plot_style(dpi: int) -> None:
    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 13,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
            "figure.dpi": dpi,
        }
    )


def plot_panel_a(
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    output_path: Path,
    *,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    all_distances = list(tp_distances) + list(fp_distances)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)
    _set_plot_style(dpi)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(
        [tp_distances, fp_distances],
        bins=bins,
        stacked=True,
        color=[classification_colors.TP, classification_colors.FP],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TP) CSP -- in binding site ({len(tp_distances)})",
            f"(FP) CSP -- not in binding site ({len(fp_distances)})",
        ],
    )
    ax = plt.gca()
    plt.xlabel("Minimum Inter-Atomic Distance (Å)")
    plt.ylabel("Number of Residues")
    ax.tick_params(axis="both", labelsize=14)
    plt.legend(frameon=True, fontsize=14, title_fontsize=14)
    plt.xlim(0, x_max)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_panel_b(
    tp_distances: Sequence[float],
    fp_distances: Sequence[float],
    fn_distances: Sequence[float],
    tn_distances: Sequence[float],
    output_path: Path,
    *,
    bin_width: float,
    max_distance: Optional[float],
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    all_distances = list(tp_distances) + list(fp_distances) + list(fn_distances) + list(tn_distances)
    bins, x_max = _get_bins(all_distances, bin_width, max_distance)
    _set_plot_style(dpi)

    plt.figure(figsize=(fig_width, fig_height))
    plt.hist(
        [tn_distances, fp_distances, fn_distances, tp_distances],
        bins=bins,
        stacked=True,
        color=[
            classification_colors.TN,
            classification_colors.FP,
            classification_colors.FN,
            classification_colors.TP,
        ],
        edgecolor="black",
        linewidth=0.8,
        label=[
            f"(TN) No CSP -- not in binding site ({len(tn_distances)})",
            f"(FP) CSP -- not in binding site ({len(fp_distances)})",
            f"(FN) No CSP -- in binding site ({len(fn_distances)})",
            f"(TP) CSP -- in binding site ({len(tp_distances)})",
        ],
    )
    ax = plt.gca()
    plt.xlabel("Minimum Inter-Atomic Distance (Å)")
    plt.ylabel("Number of Residues")
    ax.tick_params(axis="both", labelsize=14)
    plt.legend(frameon=True, fontsize=14, title_fontsize=14)
    plt.xlim(0, x_max)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def compose_two_panel_figure(panel_a_path: Path, panel_b_path: Path, output_image: Path) -> None:
    image_a = mpimg.imread(panel_a_path)
    image_b = mpimg.imread(panel_b_path)

    fig, axes = plt.subplots(2, 1, figsize=(10, 12))
    for ax in axes:
        ax.set_axis_off()

    axes[0].imshow(image_a)
    axes[1].imshow(image_b)

    axes[0].text(0.01, 0.99, "A.", transform=axes[0].transAxes, va="top", ha="left", fontsize=20, fontweight="bold")
    axes[1].text(0.01, 0.99, "B.", transform=axes[1].transAxes, va="top", ha="left", fontsize=20, fontweight="bold")

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.02)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300)
    plt.close(fig)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create SI Fig. S13 (atom–atom distance confusion matrix histograms).")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--output", type=Path, default=None, help="Output path for SF13_any_atom_distance.png")
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="CSV with 'holo_pdb' column (default: data/CSP_UBQ_ph0.5_temp5C.csv).",
    )
    parser.add_argument("--targets", type=str, help="Optional comma-separated holo_pdb list.")
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--max-distance", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--fig-width", type=float, default=8.0)
    parser.add_argument("--fig-height", type=float, default=5.0)
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    project_root = Path(__file__).resolve().parent.parent
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else project_root / args.outputs_dir
    figures_dir = args.figures_dir if args.figures_dir.is_absolute() else project_root / args.figures_dir
    output_image = args.output or (figures_dir / "SF13_any_atom_distance.png")
    if not output_image.is_absolute():
        output_image = project_root / output_image

    targets_csv = args.targets_csv
    if not targets_csv.is_absolute():
        targets_csv = project_root / targets_csv
    selected_targets = _parse_targets(targets_csv, args.targets)
    tp, fp, fn, tn = collect_distance_categories(outputs_dir, selected_targets)
    if not (tp or fp or fn or tn):
        print(
            "No inter-atomic distance data found for selected targets. "
            "Run the pipeline so outputs/<id>/master_alignment.csv exists.",
            file=sys.stderr,
        )
        return 1
    if not (tp or fp):
        print(
            "No inter-atomic distance data for significant residues (TP/FP) for selected targets.",
            file=sys.stderr,
        )
        return 1

    with tempfile.TemporaryDirectory(prefix="si_fig_s16_") as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        panel_a_path = tmp_dir_path / "panel_a.png"
        panel_b_path = tmp_dir_path / "panel_b.png"

        plot_panel_a(
            tp, fp, panel_a_path,
            bin_width=args.bin_width, max_distance=args.max_distance,
            dpi=args.dpi, fig_width=args.fig_width, fig_height=args.fig_height,
        )
        plot_panel_b(
            tp, fp, fn, tn, panel_b_path,
            bin_width=args.bin_width, max_distance=args.max_distance,
            dpi=args.dpi, fig_width=args.fig_width, fig_height=args.fig_height,
        )
        compose_two_panel_figure(panel_a_path, panel_b_path, output_image)

    print(f"SI Fig. S13 written to {output_image.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
