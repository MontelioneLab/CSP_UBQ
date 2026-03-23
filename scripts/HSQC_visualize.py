"""
HSQC-style scatter plot visualizations for apo and holo chemical shifts.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _HAS_PLT = True
except Exception:
    _HAS_PLT = False

try:
    from .csp import CSPResult, run_offset_grid_search
    from .config import Referencing as _Referencing, compute as _compute
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.csp import CSPResult, run_offset_grid_search  # type: ignore
    from scripts.config import Referencing as _Referencing, compute as _compute  # type: ignore


def _extract_shift_pairs(
    results: Sequence[CSPResult],
    h_attr: str,
    n_attr: str,
) -> List[Tuple[float, float]]:
    """Collect (N, H) shift pairs from CSP results for plotting."""
    pairs: List[Tuple[float, float]] = []
    for res in results:
        h_val = getattr(res, h_attr, None)
        n_val = getattr(res, n_attr, None)
        if h_val is not None and n_val is not None:
            pairs.append((float(n_val), float(h_val)))
    return pairs


def _compute_grid_offsets(
    results: Sequence[CSPResult],
    grid_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float]:
    """Run the grid search to obtain best (N, H) offsets."""
    params = grid_params or {}
    ref_cfg = _Referencing()
    h_min = params.get("h_min", ref_cfg.grid_h_min)
    h_max = params.get("h_max", ref_cfg.grid_h_max)
    h_step = params.get("h_step", ref_cfg.grid_h_step)
    n_min = params.get("n_min", ref_cfg.grid_n_min)
    n_max = params.get("n_max", ref_cfg.grid_n_max)
    n_step = params.get("n_step", ref_cfg.grid_n_step)
    cutoff = float(params.get("cutoff", ref_cfg.grid_cutoff))

    grid_result = run_offset_grid_search(
        list(results),
        h_min=h_min,
        h_max=h_max,
        h_step=h_step,
        n_min=n_min,
        n_max=n_max,
        n_step=n_step,
        cutoff=cutoff,
    )
    return (
        float(grid_result.get("best_n_offset", 0.0)),
        float(grid_result.get("best_h_offset", 0.0)),
    )


def plot_hsqc_variants(
    results: Sequence[CSPResult],
    out_path: str,
    *,
    title: Optional[str] = None,
    grid_params: Optional[Dict[str, float]] = None,
    apo_label: str = "Apo",
    holo_label: str = "Holo",
) -> None:
    """
    Generate a four-panel HSQC-like visualization comparing apo and holo shifts.

    Panels:
      1. Apo HSQC (apo shifts only)
      2. Holo HSQC (holo shifts only, original values)
      3. Apo vs Holo (original shifts overlay)
      4. Apo vs Holo after applying best grid offsets to holo shifts
    """
    if not _HAS_PLT:
        raise RuntimeError("matplotlib is required to generate HSQC plots.")
    if not results:
        raise ValueError("No CSP results were provided for HSQC plotting.")

    apo_pairs = _extract_shift_pairs(results, "H_apo", "N_apo")
    holo_pairs_raw = _extract_shift_pairs(results, "H_holo_original", "N_holo_original")

    if not holo_pairs_raw:
        # Fall back to referenced shifts if originals are absent
        holo_pairs_raw = _extract_shift_pairs(results, "H_holo", "N_holo")

    # Determine offsets via grid search when possible
    n_offset = 0.0
    h_offset = 0.0
    apply_offsets = bool(apo_pairs and holo_pairs_raw)
    if apply_offsets:
        n_offset, h_offset = _compute_grid_offsets(results, grid_params)
    holo_pairs_offset: List[Tuple[float, float]] = []
    if apply_offsets:
        for n_val, h_val in holo_pairs_raw:
            holo_pairs_offset.append((n_val + n_offset, h_val + h_offset))

    # Shared axis limits for consistent comparison across panels
    axis_points: List[Tuple[float, float]] = []
    axis_points.extend(apo_pairs)
    axis_points.extend(holo_pairs_raw)
    axis_points.extend(holo_pairs_offset)

    if axis_points:
        n_values = [n for n, _ in axis_points]
        h_values = [h for _, h in axis_points]
        n_margin = 0.05 * (max(n_values) - min(n_values) or 1.0)
        h_margin = 0.05 * (max(h_values) - min(h_values) or 1.0)
        n_limits = (min(n_values) - n_margin, max(n_values) + n_margin)
        h_limits = (min(h_values) - h_margin, max(h_values) + h_margin)
    else:
        n_limits = (0.0, 1.0)
        h_limits = (0.0, 1.0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    panels = axes.flatten()

    def _setup_axis(ax):
        ax.set_xlabel(r"$^{15}$N $\delta$ (ppm)", fontsize=14)
        ax.set_ylabel(r"$^{1}$H $\delta$ (ppm)", fontsize=14)
        ax.tick_params(axis="both", labelsize=12)
        ax.set_xlim(*n_limits)
        ax.set_ylim(*h_limits)
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)

    # Panel 1: Apo only
    _setup_axis(panels[0])
    if apo_pairs:
        panels[0].scatter(
            [n for n, _ in apo_pairs],
            [h for _, h in apo_pairs],
            s=16,
            c="tab:blue",
            alpha=0.8,
            edgecolors="none",
        )
    else:
        panels[0].text(
            0.5,
            0.5,
            "No apo shifts",
            transform=panels[0].transAxes,
            ha="center",
            va="center",
        )
    panels[0].set_title(f"{apo_label} HSQC")

    # Panel 2: Holo only
    _setup_axis(panels[1])
    if holo_pairs_raw:
        panels[1].scatter(
            [n for n, _ in holo_pairs_raw],
            [h for _, h in holo_pairs_raw],
            s=16,
            c="tab:orange",
            alpha=0.8,
            edgecolors="none",
        )
    else:
        panels[1].text(
            0.5,
            0.5,
            "No holo shifts",
            transform=panels[1].transAxes,
            ha="center",
            va="center",
        )
    panels[1].set_title(f"{holo_label} HSQC")

    # Panel 3: Combined raw
    _setup_axis(panels[2])
    if apo_pairs or holo_pairs_raw:
        if apo_pairs:
            panels[2].scatter(
                [n for n, _ in apo_pairs],
                [h for _, h in apo_pairs],
                s=18,
                c="tab:blue",
                alpha=0.6,
                edgecolors="none",
                label=apo_label,
            )
        if holo_pairs_raw:
            panels[2].scatter(
                [n for n, _ in holo_pairs_raw],
                [h for _, h in holo_pairs_raw],
                s=18,
                c="tab:orange",
                alpha=0.6,
                edgecolors="none",
                label=holo_label,
            )
        panels[2].legend()
    else:
        panels[2].text(
            0.5,
            0.5,
            "No shifts to display",
            transform=panels[2].transAxes,
            ha="center",
            va="center",
        )
    panels[2].set_title("Apo vs Holo (raw)")

    # Panel 4: Combined with offsets
    _setup_axis(panels[3])
    if apo_pairs or holo_pairs_offset:
        if apo_pairs:
            panels[3].scatter(
                [n for n, _ in apo_pairs],
                [h for _, h in apo_pairs],
                s=18,
                c="tab:blue",
                alpha=0.6,
                edgecolors="none",
                label=apo_label,
            )
        if holo_pairs_offset:
            panels[3].scatter(
                [n for n, _ in holo_pairs_offset],
                [h for _, h in holo_pairs_offset],
                s=18,
                c="tab:orange",
                alpha=0.6,
                edgecolors="none",
                label=f"{holo_label} (offset)"
                if (h_offset or n_offset)
                else holo_label,
            )
        if panels[3].has_data():
            legend = panels[3].legend()
            if legend and (h_offset or n_offset):
                legend.set_title(
                    f"Offsets: ΔH={h_offset:.3f}, ΔN={n_offset:.3f}"
                )
    else:
        panels[3].text(
            0.5,
            0.5,
            "No shifts to display",
            transform=panels[3].transAxes,
            ha="center",
            va="center",
        )
    # This panel is reused in case-study figures; keep it title-free.
    panels[3].set_title("")

    if title:
        fig.suptitle(title, fontsize=16)

    fig.tight_layout(rect=[0, 0, 1, 0.97] if title else None)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


