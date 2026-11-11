#!/usr/bin/env python3
"""Render the 2×2 SI Fig. S17 from the frozen CSP threshold analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    sys.path.insert(0, str(REPO))

from scripts.create_si_fig_s27 import (
    _load,
    draw_pooled_intensity,
    draw_roc,
    set_style,
    validate_analysis,
)
from scripts.plot_csp_z_vs_ca_distance import (
    ANY_ATOM_DISTANCE_XLABEL,
    _draw_csp_z_vs_distance_scatter,
    collect_points_from_master_csv,
)
from scripts.plot_p_significant_vs_ca_distance import (
    DEFAULT_MAX_DISTANCE_A,
    DISTANCE_XLABEL,
    collect_distance_significance,
    compute_stacked_bin_probabilities,
    draw_stacked_bin_probabilities,
)

import matplotlib.pyplot as plt

DEFAULT_ANALYSIS = REPO / "outputs" / "2026-11-04_csp-threshold-analysis_002"
DEFAULT_OUTPUT = REPO / "figures" / "SF17_csp_distance_and_binding_site.pdf"
S17D_KL_TEST_DATA_SHA256 = "9ef05b0c00bbc2373f56a1ff5d9fa4f7eea76eb60775cc58a651db0390277026"


def render(analysis_dir: Path, output: Path, *, intensity_scale: str = "linear") -> None:
    manifest = validate_analysis(analysis_dir)
    sources = manifest["sources"]
    if any(not source.get("snapshot_path") for source in sources):
        raise ValueError("S17 requires frozen master alignment tables in the analysis snapshot.")
    selected = {source["target_id"] for source in sources}
    distances, significance, z_scores, n_targets = collect_distance_significance(
        analysis_dir / "inputs", selected
    )
    if n_targets != manifest["n_targets"]:
        raise ValueError(f"S17A distance data cover {n_targets} of {manifest['n_targets']} targets.")
    bin_width = 1.0
    bins = compute_stacked_bin_probabilities(
        distances, significance, z_scores, bin_width=bin_width, max_distance=DEFAULT_MAX_DISTANCE_A
    )
    points = [
        point
        for source in sorted(sources, key=lambda entry: entry["target_id"])
        for point in collect_points_from_master_csv(str(analysis_dir / source["snapshot_path"]))
    ]

    set_style()
    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(7.3, 5.85))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.105, hspace=0.43, wspace=0.33)

    draw_stacked_bin_probabilities(ax_a, bins, bin_width=bin_width)
    ax_a.set(xlabel=DISTANCE_XLABEL, ylabel="Significant CSPs (%)")
    ax_a.set_xticks([0, 10, 20, 30, 40])
    ax_a.legend(loc="upper right", title="Significant CSP z-score", title_fontsize=7)

    _draw_csp_z_vs_distance_scatter(
        ax_b, points, xlabel=ANY_ATOM_DISTANCE_XLABEL, legend_fontsize=6.5,
        label_fontsize=8.5, marker_size=6, fp_text_fontsize=7, compact_legend=True,
    )

    draw_roc(
        ax_c, _load(analysis_dir, "roc_curves.csv"), _load(analysis_dir, "macro_roc.csv"),
        _load(analysis_dir, "aggregate_metrics.csv"), _load(analysis_dir, "target_summary.csv"),
    )
    draw_pooled_intensity(ax_d, _load(analysis_dir, "violin_residues.csv"))
    ax_d.set_xticks([0, 1], ["TP\nBinding site", "FP\nDistal"])
    if intensity_scale == "linear":
        ax_d.set(yscale="linear", ylim=(0, 2.5))
        ax_d.set_yticks([0, 0.5, 1, 1.5, 2, 2.5], ["0", "0.5", "1", "1.5", "2", "2.5"])
    else:
        ax_d.set_yticks([0.1, 0.2, 0.5, 1, 2], ["0.1", "0.2", "0.5", "1", "2"])
    # The permutation result applies to this exact frozen TP/FP dataset.
    if manifest["data_sha256"].get("violin_residues.csv") == S17D_KL_TEST_DATA_SHA256:
        ax_d.text(
            0.5, 0.84, r"$p = 1.0 \times 10^{-5}$",
            transform=ax_d.transAxes, ha="center", va="center",
        )

    for ax, title in zip(
        (ax_a, ax_b, ax_c, ax_d),
        ("A.  Significant CSPs by distance", "B.  CSP z-score by distance",
         "C.  Binding-site localization", "D.  Significant CSP magnitudes"),
    ):
        ax.set_title(title, loc="left", pad=7)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Title": "SI Fig. S17 - CSP distance and binding-site localization"})
    plt.close(fig)
    print(f"Wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--intensity-scale", choices=("linear", "log"), default="linear",
                        help="Display scale for panel D (default: linear).")
    args = parser.parse_args()
    render(args.analysis_dir, args.output, intensity_scale=args.intensity_scale)


if __name__ == "__main__":
    main()
