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
    from .csp import CSPResult, run_offset_grid_search, run_offset_grid_search_ha_ca
    from .config import Referencing as _Referencing, compute as _compute
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.csp import CSPResult, run_offset_grid_search, run_offset_grid_search_ha_ca  # type: ignore
    from scripts.config import Referencing as _Referencing, compute as _compute  # type: ignore


def _extract_shift_pairs(
    results: Sequence[CSPResult],
    x_attr: str,
    y_attr: str,
) -> List[Tuple[float, float]]:
    """Collect (x, y) shift pairs from CSP results for plotting."""
    pairs: List[Tuple[float, float]] = []
    for res in results:
        x_val = getattr(res, x_attr, None)
        y_val = getattr(res, y_attr, None)
        if x_val is not None and y_val is not None:
            pairs.append((float(x_val), float(y_val)))
    return pairs


def _compute_grid_offsets(
    results: Sequence[CSPResult],
    grid_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float]:
    """Run the N/H grid search to obtain best (N, H) offsets."""
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


def _compute_grid_offsets_ha_ca(
    results: Sequence[CSPResult],
    grid_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float]:
    """Run the HA/CA grid search to obtain best (CA, HA) offsets."""
    params = grid_params or {}
    ref_cfg = _Referencing()
    ha_min = params.get("ha_min", ref_cfg.grid_ha_min)
    ha_max = params.get("ha_max", ref_cfg.grid_ha_max)
    ha_step = params.get("ha_step", ref_cfg.grid_ha_step)
    ca_min = params.get("ca_min", ref_cfg.grid_ca_min)
    ca_max = params.get("ca_max", ref_cfg.grid_ca_max)
    ca_step = params.get("ca_step", ref_cfg.grid_ca_step)
    cutoff = float(params.get("cutoff", ref_cfg.grid_cutoff))

    points: List[Tuple[float, float, float, float]] = []
    for res in results:
        if (
            getattr(res, "HA_apo", None) is not None
            and getattr(res, "CA_apo", None) is not None
            and getattr(res, "HA_holo_original", None) is not None
            and getattr(res, "CA_holo_original", None) is not None
        ):
            points.append(
                (
                    float(res.HA_apo),  # type: ignore[arg-type]
                    float(res.CA_apo),  # type: ignore[arg-type]
                    float(res.HA_holo_original),  # type: ignore[arg-type]
                    float(res.CA_holo_original),  # type: ignore[arg-type]
                )
            )

    grid_result = run_offset_grid_search_ha_ca(
        points,
        ha_min=ha_min,
        ha_max=ha_max,
        ha_step=ha_step,
        ca_min=ca_min,
        ca_max=ca_max,
        ca_step=ca_step,
        cutoff=cutoff,
    )
    return (
        float(grid_result.get("best_ca_offset", 0.0)),
        float(grid_result.get("best_ha_offset", 0.0)),
    )


def plot_hsqc_variants(
    results: Sequence[CSPResult],
    out_path: str,
    *,
    title: Optional[str] = None,
    grid_params: Optional[Dict[str, float]] = None,
    apo_label: str = "Apo",
    holo_label: str = "Holo",
    x_atom: str = "N",
    y_atom: str = "H",
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

    x_atom = x_atom.upper()
    y_atom = y_atom.upper()
    apo_pairs = _extract_shift_pairs(results, f"{x_atom}_apo", f"{y_atom}_apo")
    holo_pairs_raw = _extract_shift_pairs(results, f"{x_atom}_holo_original", f"{y_atom}_holo_original")

    if not holo_pairs_raw:
        # Fall back to referenced shifts if originals are absent
        holo_pairs_raw = _extract_shift_pairs(results, f"{x_atom}_holo", f"{y_atom}_holo")

    # Determine offsets via grid search when possible
    x_offset = 0.0
    y_offset = 0.0
    apply_offsets = bool(apo_pairs and holo_pairs_raw)
    if apply_offsets:
        if (x_atom, y_atom) == ("N", "H"):
            x_offset, y_offset = _compute_grid_offsets(results, grid_params)
        elif (x_atom, y_atom) == ("CA", "HA"):
            x_offset, y_offset = _compute_grid_offsets_ha_ca(results, grid_params)
    holo_pairs_offset: List[Tuple[float, float]] = []
    if apply_offsets:
        for x_val, y_val in holo_pairs_raw:
            holo_pairs_offset.append((x_val + x_offset, y_val + y_offset))

    # Shared axis limits for consistent comparison across panels
    axis_points: List[Tuple[float, float]] = []
    axis_points.extend(apo_pairs)
    axis_points.extend(holo_pairs_raw)
    axis_points.extend(holo_pairs_offset)

    if axis_points:
        x_values = [x for x, _ in axis_points]
        y_values = [y for _, y in axis_points]
        x_margin = 0.05 * (max(x_values) - min(x_values) or 1.0)
        y_margin = 0.05 * (max(y_values) - min(y_values) or 1.0)
        n_limits = (min(x_values) - x_margin, max(x_values) + x_margin)
        h_limits = (min(y_values) - y_margin, max(y_values) + y_margin)
    else:
        n_limits = (0.0, 1.0)
        h_limits = (0.0, 1.0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    panels = axes.flatten()

    def _setup_axis(ax):
        labels = {
            "N": r"$^{15}$N $\delta$ (ppm)",
            "H": r"$^{1}$H $\delta$ (ppm)",
            "CA": r"$C_\alpha$ $\delta$ (ppm)",
            "HA": r"$H_\alpha$ $\delta$ (ppm)",
        }
        ax.set_xlabel(labels.get(x_atom, f"{x_atom} shift (ppm)"), fontsize=14)
        ax.set_ylabel(labels.get(y_atom, f"{y_atom} shift (ppm)"), fontsize=14)
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
                if (x_offset or y_offset)
                else holo_label,
            )
        if panels[3].has_data():
            legend = panels[3].legend()
            if legend and (x_offset or y_offset):
                legend.set_title(
                    f"Offsets: Δ{y_atom}={y_offset:.3f}, Δ{x_atom}={x_offset:.3f}"
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

