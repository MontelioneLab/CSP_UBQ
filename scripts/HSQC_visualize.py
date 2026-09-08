"""
HSQC-style scatter plot visualizations for apo and holo chemical shifts.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    _HAS_PLT = True
except Exception:
    _HAS_PLT = False

try:
    from .csp import CSPResult, run_offset_grid_search, run_offset_grid_search_ha_ca
    from .config import Referencing as _Referencing, classification_colors
    from .merge_csv import exceeds_max_classification_csp_z
except Exception:
    import os as _os
    import sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.csp import CSPResult, run_offset_grid_search, run_offset_grid_search_ha_ca  # type: ignore
    from scripts.config import Referencing as _Referencing, classification_colors  # type: ignore
    from scripts.merge_csv import exceeds_max_classification_csp_z  # type: ignore

# Apo / holo point colors for overlay panels
_APO_COLOR = "#7f7f7f"
_HOLO_COLOR = "#e74c3c"


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


def _extract_paired_shifts(
    results: Sequence[CSPResult],
    x_atom: str,
    y_atom: str,
) -> List[Tuple[CSPResult, Tuple[float, float], Tuple[float, float]]]:
    """
    Collect per-residue apo/holo (x, y) pairs when both ends are present.

    Prefers holo *_holo_original attributes, falling back to referenced *_holo.
    """
    apo_x_attr = f"{x_atom}_apo"
    apo_y_attr = f"{y_atom}_apo"
    holo_x_orig = f"{x_atom}_holo_original"
    holo_y_orig = f"{y_atom}_holo_original"
    holo_x_ref = f"{x_atom}_holo"
    holo_y_ref = f"{y_atom}_holo"

    paired: List[Tuple[CSPResult, Tuple[float, float], Tuple[float, float]]] = []
    for res in results:
        apo_x = getattr(res, apo_x_attr, None)
        apo_y = getattr(res, apo_y_attr, None)
        if apo_x is None or apo_y is None:
            continue

        holo_x = getattr(res, holo_x_orig, None)
        holo_y = getattr(res, holo_y_orig, None)
        if holo_x is None or holo_y is None:
            holo_x = getattr(res, holo_x_ref, None)
            holo_y = getattr(res, holo_y_ref, None)
        if holo_x is None or holo_y is None:
            continue

        paired.append(
            (
                res,
                (float(apo_x), float(apo_y)),
                (float(holo_x), float(holo_y)),
            )
        )
    return paired


# Neutral color for apo→holo lines when CSP significance was not recorded
# (e.g. |ΔN| > max_abs_delta_n_ppm exclusions). Not a TP/FP/TN/FN label.
_EXCLUDED_LINE_COLOR = "#808080"


def _classify_residue(
    res: CSPResult,
    *,
    position_map: Dict[int, int],
    binding_lookup: Dict[int, bool],
    significance_field: str,
) -> str:
    """Return TP/FP/TN/FN, or EXCLUDED when significance was not recorded.

    Matches merge_csv.compute_classification: missing CSP significance must not
    receive a confusion-matrix label (including FN for binding residues).
    Extreme CSP z-scores (> max_classification_csp_z) also return EXCLUDED so
    apo→holo connectors are drawn gray.
    """
    if exceeds_max_classification_csp_z(getattr(res, "z_score", None)):
        return "EXCLUDED"

    significance = getattr(res, significance_field, None)
    if significance is None:
        return "EXCLUDED"

    pdb_residue_number = position_map.get(res.holo_index)
    is_binding = (
        False
        if pdb_residue_number is None
        else binding_lookup.get(pdb_residue_number, False)
    )
    is_significant = bool(significance)

    if is_significant and is_binding:
        return "TP"
    if is_significant and not is_binding:
        return "FP"
    if not is_significant and not is_binding:
        return "TN"
    return "FN"


def _classification_color_map() -> Dict[str, str]:
    return {
        "TP": classification_colors.TP,
        "FP": classification_colors.FP,
        "TN": classification_colors.TN,
        "FN": classification_colors.FN,
        "EXCLUDED": _EXCLUDED_LINE_COLOR,
    }


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


def _first_applied_offset_pair(
    results: Sequence[CSPResult],
    x_attr: str,
    y_attr: str,
) -> Optional[Tuple[float, float]]:
    """Return the first non-None (x_offset, y_offset) pair from CSPResult attrs."""
    for res in results:
        x_off = getattr(res, x_attr, None)
        y_off = getattr(res, y_attr, None)
        if x_off is not None and y_off is not None:
            return float(x_off), float(y_off)
    return None


def format_offset_annotation(x_atom: str, x_offset: float, y_atom: str, y_offset: float) -> str:
    """Signed ppm label for applied referencing offsets (e.g. ΔH=+0.010 ppm, ΔN=-0.050 ppm)."""
    return f"Δ{x_atom}={x_offset:+.3f} ppm, Δ{y_atom}={y_offset:+.3f} ppm"


def resolve_hsqc_offsets(
    results: Sequence[CSPResult],
    x_atom: str,
    y_atom: str,
    grid_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float]:
    """
    Resolve (x_offset, y_offset) for the plotted atom pair.

    Prefers applied offsets already stored on CSPResult; falls back to grid search.
    """
    x_atom = x_atom.upper()
    y_atom = y_atom.upper()

    if (x_atom, y_atom) == ("H", "N"):
        applied = _first_applied_offset_pair(results, "H_offset", "N_offset")
        if applied is not None:
            return applied
        n_offset, h_offset = _compute_grid_offsets(results, grid_params)
        return h_offset, n_offset

    if (x_atom, y_atom) == ("N", "H"):
        applied = _first_applied_offset_pair(results, "N_offset", "H_offset")
        if applied is not None:
            return applied
        return _compute_grid_offsets(results, grid_params)

    if (x_atom, y_atom) == ("HA", "CA"):
        applied = _first_applied_offset_pair(results, "HA_offset", "CA_offset")
        if applied is not None:
            return applied
        ca_offset, ha_offset = _compute_grid_offsets_ha_ca(results, grid_params)
        return ha_offset, ca_offset

    if (x_atom, y_atom) == ("CA", "HA"):
        applied = _first_applied_offset_pair(results, "CA_offset", "HA_offset")
        if applied is not None:
            return applied
        return _compute_grid_offsets_ha_ca(results, grid_params)

    # Unknown atom pair: try matching *_offset attrs, else zero.
    applied = _first_applied_offset_pair(results, f"{x_atom}_offset", f"{y_atom}_offset")
    if applied is not None:
        return applied
    return 0.0, 0.0


def _draw_overlay_panel(
    ax,
    *,
    apo_points: List[Tuple[float, float]],
    holo_points: List[Tuple[float, float]],
    segments: List[Tuple[Tuple[float, float], Tuple[float, float]]],
    segment_colors: List[str],
    apo_label: str,
    holo_label: str,
    legend_title: Optional[str] = None,
) -> None:
    """Draw classification-colored apo→holo lines under gray apo / red holo points."""
    if segments and segment_colors:
        lc = LineCollection(
            segments,
            colors=segment_colors,
            linewidths=1.0,
            alpha=0.75,
            zorder=1,
        )
        ax.add_collection(lc)

    if apo_points:
        ax.scatter(
            [x for x, _ in apo_points],
            [y for _, y in apo_points],
            s=18,
            c=_APO_COLOR,
            alpha=0.8,
            edgecolors="none",
            label=apo_label,
            zorder=2,
        )
    if holo_points:
        ax.scatter(
            [x for x, _ in holo_points],
            [y for _, y in holo_points],
            s=18,
            c=_HOLO_COLOR,
            alpha=0.8,
            edgecolors="none",
            label=holo_label,
            zorder=2,
        )

    if apo_points or holo_points or segments:
        legend = ax.legend()
        if legend and legend_title:
            legend.set_title(legend_title)


def plot_hsqc_variants(
    results: Sequence[CSPResult],
    out_path: str,
    *,
    title: Optional[str] = None,
    grid_params: Optional[Dict[str, float]] = None,
    apo_label: str = "Apo",
    holo_label: str = "Holo",
    x_atom: str = "H",
    y_atom: str = "N",
    binding_results: Optional[dict] = None,
    significance_field: str = "significant",
) -> Tuple[float, float]:
    """
    Generate a four-panel HSQC-like visualization comparing apo and holo shifts.

    Axes follow conventional 2D NMR orientation (descending ppm):
      X = direct-dimension nucleus (default ¹H), high ppm left → low right
      Y = indirect-dimension nucleus (default ¹⁵N), high ppm bottom → low top

    Panels:
      1. Apo HSQC (apo shifts only)
      2. Holo HSQC (holo shifts only, original values)
      3. Apo vs Holo (original shifts overlay) with classification-colored lines
      4. Apo vs Holo after applying referencing offsets to holo shifts
         (same styling; reused in case-study figures)

    Returns:
        (x_offset, y_offset) applied to holo shifts in panel 4.
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

    paired = _extract_paired_shifts(results, x_atom, y_atom)

    # Classification lookup (optional; default non-binding when absent)
    position_map: Dict[int, int] = {}
    binding_lookup: Dict[int, bool] = {}
    if binding_results:
        try:
            from .visualize import (
                create_sequence_alignment_map_from_results,
                _build_binding_lookup,
            )
        except Exception:
            from scripts.visualize import (  # type: ignore
                create_sequence_alignment_map_from_results,
                _build_binding_lookup,
            )
        position_map = create_sequence_alignment_map_from_results(list(results), binding_results)
        binding_lookup, _ = _build_binding_lookup(binding_results)

    class_colors = _classification_color_map()

    # Prefer applied CSPResult offsets; fall back to grid search when absent.
    x_offset = 0.0
    y_offset = 0.0
    if apo_pairs and holo_pairs_raw:
        x_offset, y_offset = resolve_hsqc_offsets(
            results, x_atom, y_atom, grid_params=grid_params
        )

    holo_pairs_offset: List[Tuple[float, float]] = [
        (x_val + x_offset, y_val + y_offset) for x_val, y_val in holo_pairs_raw
    ]

    # Shared axis limits for consistent comparison across panels.
    # Descending ppm: high left / high bottom (conventional NMR HSQC).
    axis_points: List[Tuple[float, float]] = []
    axis_points.extend(apo_pairs)
    axis_points.extend(holo_pairs_raw)
    axis_points.extend(holo_pairs_offset)

    if axis_points:
        x_values = [x for x, _ in axis_points]
        y_values = [y for _, y in axis_points]
        x_margin = 0.05 * (max(x_values) - min(x_values) or 1.0)
        y_margin = 0.05 * (max(y_values) - min(y_values) or 1.0)
        x_limits = (max(x_values) + x_margin, min(x_values) - x_margin)
        y_limits = (max(y_values) + y_margin, min(y_values) - y_margin)
    else:
        x_limits = (1.0, 0.0)
        y_limits = (1.0, 0.0)

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
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)

    # Panel 1: Apo only
    _setup_axis(panels[0])
    if apo_pairs:
        panels[0].scatter(
            [n for n, _ in apo_pairs],
            [h for _, h in apo_pairs],
            s=16,
            c=_APO_COLOR,
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
            c=_HOLO_COLOR,
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

    # Build overlay segment data once; panel 4 applies offsets to holo ends.
    raw_apo_pts: List[Tuple[float, float]] = []
    raw_holo_pts: List[Tuple[float, float]] = []
    raw_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    offset_apo_pts: List[Tuple[float, float]] = []
    offset_holo_pts: List[Tuple[float, float]] = []
    offset_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    segment_colors: List[str] = []

    for res, apo_xy, holo_xy in paired:
        cls = _classify_residue(
            res,
            position_map=position_map,
            binding_lookup=binding_lookup,
            significance_field=significance_field,
        )
        color = class_colors[cls]
        holo_offset_xy = (holo_xy[0] + x_offset, holo_xy[1] + y_offset)

        raw_apo_pts.append(apo_xy)
        raw_holo_pts.append(holo_xy)
        raw_segments.append((apo_xy, holo_xy))

        offset_apo_pts.append(apo_xy)
        offset_holo_pts.append(holo_offset_xy)
        offset_segments.append((apo_xy, holo_offset_xy))
        segment_colors.append(color)

    # Also show unpaired apo/holo points on overlays (no connecting line).
    paired_apo_set = set(raw_apo_pts)
    paired_holo_set = set(raw_holo_pts)
    for pt in apo_pairs:
        if pt not in paired_apo_set:
            raw_apo_pts.append(pt)
            offset_apo_pts.append(pt)
    for pt in holo_pairs_raw:
        if pt not in paired_holo_set:
            raw_holo_pts.append(pt)
            offset_holo_pts.append((pt[0] + x_offset, pt[1] + y_offset))

    # Panel 3: Combined raw
    _setup_axis(panels[2])
    if raw_apo_pts or raw_holo_pts or raw_segments:
        _draw_overlay_panel(
            panels[2],
            apo_points=raw_apo_pts,
            holo_points=raw_holo_pts,
            segments=raw_segments,
            segment_colors=segment_colors,
            apo_label=apo_label,
            holo_label=holo_label,
        )
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
    offset_label = format_offset_annotation(x_atom, x_offset, y_atom, y_offset)
    offset_holo_label = f"{holo_label} (offset)"
    if offset_apo_pts or offset_holo_pts or offset_segments:
        _draw_overlay_panel(
            panels[3],
            apo_points=offset_apo_pts,
            holo_points=offset_holo_pts,
            segments=offset_segments,
            segment_colors=segment_colors,
            apo_label=apo_label,
            holo_label=offset_holo_label,
            legend_title=offset_label,
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
    # Applied offsets are shown only in the legend title.
    panels[3].set_title("")

    if title:
        fig.suptitle(title, fontsize=16)

    fig.tight_layout(rect=[0, 0, 1, 0.97] if title else None)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return x_offset, y_offset
