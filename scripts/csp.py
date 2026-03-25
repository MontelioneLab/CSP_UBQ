"""
Compute chemical shift perturbations (CSP) using equation A and summarize.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import os
import math
import statistics

try:
    from .config import compute, thresholds, paths, Referencing as _Referencing  # type: ignore
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import compute, thresholds, paths, Referencing as _Referencing  # type: ignore


def _csp_n_term(d_n: float) -> float:
    w = compute.csp_delta_n_scale
    return (d_n * w) ** 2


def _csp_ca_term(d_ca: float) -> float:
    w = compute.csp_delta_ca_scale
    return (d_ca * w) ** 2


def compute_threshold_with_outlier_removal(
    values: List[float],
    outlier_z: float,
    significance_z: float,
    max_iterations: int,
    max_outlier_fraction: float
) -> ThresholdComputation:
    """
    Compute CSP significance threshold using iterative outlier removal.
    
    Parameters
    ----------
    values : List[float]
        List of CSP values
    outlier_z : float
        Z-score threshold for outlier detection
    significance_z : float
        Z-score for final threshold calculation
    max_iterations : int
        Maximum number of iterations
    max_outlier_fraction : float
        Maximum fraction of values that can be removed as outliers
        
    Returns
    -------
    Tuple[float, int, int]
        (threshold, num_iterations, num_outliers_removed)
    """
    if not values:
        return ThresholdComputation(
            threshold=0.0,
            iterations=0,
            outliers_removed=0,
            cleaned_values=[],
            mean=0.0,
            sd=0.0,
        )
    
    # Start with all values
    current_values = values.copy()
    original_count = len(values)
    max_outliers = int(original_count * max_outlier_fraction)
    
    iteration_count = 0
    for iteration in range(max_iterations):
        if len(current_values) <= 1:
            break
        iteration_count = iteration + 1
            
        # Calculate mean and standard deviation
        mean_val = statistics.mean(current_values)
        sd_val = statistics.pstdev(current_values) if len(current_values) > 1 else 0.0
        
        # Find outliers (values > mean + outlier_z * SD)
        outlier_threshold = mean_val + outlier_z * sd_val
        outliers = [v for v in current_values if v > outlier_threshold]
        
        # Check if we should stop
        if not outliers:
            # No outliers found, convergence achieved
            break
            
        # Check if removing outliers would exceed max fraction
        if len(outliers) > max_outliers:
            # Stop iteration to avoid removing too many values
            break
            
        # Remove outliers for next iteration
        current_values = [v for v in current_values if v <= outlier_threshold]
    
    # Calculate final threshold from cleaned dataset
    final_values = current_values.copy()
    final_mean, final_sd = _compute_basic_stats(final_values)
    threshold = final_mean + significance_z * final_sd
    
    outliers_removed = original_count - len(current_values)
    
    return ThresholdComputation(
        threshold=threshold,
        iterations=iteration_count,
        outliers_removed=outliers_removed,
        cleaned_values=final_values,
        mean=final_mean,
        sd=final_sd,
    )


@dataclass
class ThresholdComputation:
    threshold: float
    iterations: int
    outliers_removed: int
    cleaned_values: List[float]
    mean: float
    sd: float


def _compute_basic_stats(values: List[float]) -> Tuple[float, float]:
    if not values:
        return 0.0, 0.0
    mean_val = statistics.mean(values)
    sd_val = statistics.pstdev(values) if len(values) > 1 else 0.0
    return mean_val, sd_val


def _make_absolute_threshold(values: List[float], cutoff: float) -> ThresholdComputation:
    cleaned = values.copy()
    mean_val, sd_val = _compute_basic_stats(cleaned)
    return ThresholdComputation(
        threshold=float(cutoff),
        iterations=0,
        outliers_removed=0,
        cleaned_values=cleaned,
        mean=mean_val,
        sd=sd_val,
    )


@dataclass
class CSPResult:
    # per-residue mapped values; indices are 1-based within their sequences
    apo_index: int
    holo_index: int
    apo_aa: str
    holo_aa: str
    H_apo: Optional[float]
    N_apo: Optional[float]
    H_holo: Optional[float]
    N_holo: Optional[float]
    dH: Optional[float]
    dN: Optional[float]
    csp_A: Optional[float]
    significant: Optional[bool]
    # Additional significance thresholds
    significant_1sd: Optional[bool] = None
    significant_2sd: Optional[bool] = None
    # Referencing fields
    H_holo_original: Optional[float] = None
    N_holo_original: Optional[float] = None
    H_offset: Optional[float] = None
    N_offset: Optional[float] = None
    z_score: Optional[float] = None
    # Optional CA-related fields (used for CA-inclusive CSP analysis)
    CA_apo: Optional[float] = None
    CA_holo: Optional[float] = None
    CA_holo_original: Optional[float] = None
    CA_offset: Optional[float] = None
    dCA: Optional[float] = None


@dataclass
class AtomDeltaResult:
    """
    Container for per-atom perturbation analysis derived from CSPResult.

    This is used for 1D offset grid searches and per-atom thresholding that
    treat each nucleus (H, N, CA) independently from the combined CSP(A) metric.
    """

    atom: str
    offset: float
    cutoff: float
    values_abs: List[float]
    per_residue: List[Dict[str, object]]


def _frange(min_value: float, max_value: float, step: float, decimals: int) -> List[float]:
    """Generate an inclusive floating range with stable rounding."""
    if step <= 0:
        return []
    values: List[float] = []
    # Use integer steps to avoid accumulation error
    count = int(round((max_value - min_value) / step))
    for i in range(count + 1):
        v = min_value + i * step
        values.append(round(v, decimals))
    # Ensure max is included exactly
    if values and values[-1] != round(max_value, decimals):
        values.append(round(max_value, decimals))
    return values


def run_offset_grid_search_1d(
    points: List[Tuple[float, float]],
    *,
    min_offset: float,
    max_offset: float,
    step: float,
    cutoff: float,
) -> Dict[str, object]:
    """
    1D grid search for a single atom type.

    Parameters
    ----------
    points : List[Tuple[float, float]]
        List of (apo_shift, holo_shift_original) tuples.
    min_offset, max_offset, step : float
        Offset grid parameters. The search is inclusive of the end points.
    cutoff : float
        Absolute delta cutoff; counts residues with |(holo+offset) - apo| < cutoff.
    """
    if not points or step <= 0:
        return {
            "offset_values": [],
            "counts": [],
            "best_offset": 0.0,
            "best_count": 0,
            "n_points": len(points),
        }

    decimals = max(0, min(6, len(str(step).split(".")[-1]) if "." in str(step) else 0))
    offset_values = _frange(min_offset, max_offset, step, decimals)

    best_count = -1
    best_offsets: List[float] = []
    counts: List[int] = []

    for off in offset_values:
        count_lt = 0
        for apo_shift, holo_orig in points:
            delta = (holo_orig + off) - apo_shift
            if abs(delta) < cutoff:
                count_lt += 1
        counts.append(count_lt)
        if count_lt > best_count:
            best_count = count_lt
            best_offsets = [off]
        elif count_lt == best_count:
            best_offsets.append(off)

    if best_offsets:
        best_offset = min(best_offsets, key=lambda o: (abs(o), o))
    else:
        best_offset = 0.0

    return {
        "offset_values": offset_values,
        "counts": counts,
        "best_offset": best_offset,
        "best_count": best_count if best_count >= 0 else 0,
        "n_points": len(points),
    }


def run_offset_grid_search(
    results: List[CSPResult],
    *,
    h_min: float,
    h_max: float,
    h_step: float,
    n_min: float,
    n_max: float,
    n_step: float,
    cutoff: float,
) -> Dict[str, object]:
    """
    Grid search for N/H offsets that maximize the number of CSPs < cutoff.

    Returns dict with keys:
      - h_values: List[float]
      - n_values: List[float]
      - count_matrix: List[List[int]] with shape (len(h_values), len(n_values))
      - best_h_offset: float
      - best_n_offset: float
      - best_count: int
    """
    # Determine reasonable rounding for axes
    h_decimals = max(0, min(6, len(str(h_step).split(".")[-1]) if "." in str(h_step) else 0))
    n_decimals = max(0, min(6, len(str(n_step).split(".")[-1]) if "." in str(n_step) else 0))

    h_values = _frange(h_min, h_max, h_step, h_decimals)
    n_values = _frange(n_min, n_max, n_step, n_decimals)

    count_matrix: List[List[int]] = []

    best_count = -1
    best_pairs: List[Tuple[float, float]] = []

    for h_off in h_values:
        row: List[int] = []
        for n_off in n_values:
            count_lt = 0
            for r in results:
                if (
                    r.H_apo is not None and r.N_apo is not None and
                    r.H_holo_original is not None and r.N_holo_original is not None
                ):
                    dH = (r.H_holo_original + h_off) - r.H_apo
                    dN = (r.N_holo_original + n_off) - r.N_apo
                    csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
                    if csp_val < cutoff:
                        count_lt += 1
            row.append(count_lt)

            if count_lt > best_count:
                best_count = count_lt
                best_pairs = [(n_off, h_off)]
            elif count_lt == best_count:
                best_pairs.append((n_off, h_off))
        count_matrix.append(row)

    # Tie-breaker: pick pair with smallest |H| + |N|
    if best_pairs:
        best_n_offset, best_h_offset = min(best_pairs, key=lambda p: (abs(p[0]) + abs(p[1]), abs(p[0]), abs(p[1])))
    else:
        best_n_offset, best_h_offset = 0.0, 0.0

    return {
        "h_values": h_values,
        "n_values": n_values,
        "count_matrix": count_matrix,
        "best_h_offset": best_h_offset,
        "best_n_offset": best_n_offset,
        "best_count": best_count if best_count >= 0 else 0,
    }


def run_offset_grid_search_3d(
    points: List[Tuple[float, float, float, float, float, float]],
    *,
    h_min: float,
    h_max: float,
    h_step: float,
    n_min: float,
    n_max: float,
    n_step: float,
    ca_min: float,
    ca_max: float,
    ca_step: float,
    cutoff: float,
) -> Dict[str, object]:
    """
    Grid search for H/N/CA offsets that maximize the number of CA-inclusive CSPs < cutoff.

    Each point is a tuple:
      (H_apo, N_apo, CA_apo, H_holo_orig, N_holo_orig, CA_holo_orig)
    using the CA-inclusive CSP formula:
      CSP = sqrt(1/3*(dH^2+(wN*dN)^2+(wCA*dCA)^2)), wN/wCA from config Compute
    """
    if not points:
        return {
            "best_h_offset": 0.0,
            "best_n_offset": 0.0,
            "best_ca_offset": 0.0,
            "best_count": 0,
        }

    # Determine reasonable rounding for axes
    h_decimals = max(0, min(6, len(str(h_step).split(".")[-1]) if "." in str(h_step) else 0))
    n_decimals = max(0, min(6, len(str(n_step).split(".")[-1]) if "." in str(n_step) else 0))
    ca_decimals = max(0, min(6, len(str(ca_step).split(".")[-1]) if "." in str(ca_step) else 0))

    h_values = _frange(h_min, h_max, h_step, h_decimals)
    n_values = _frange(n_min, n_max, n_step, n_decimals)
    ca_values = _frange(ca_min, ca_max, ca_step, ca_decimals)

    best_count = -1
    best_triplets: List[Tuple[float, float, float]] = []

    cutoff_sq = cutoff * cutoff

    for h_off in h_values:
        for n_off in n_values:
            for ca_off in ca_values:
                count_lt = 0
                for h_a, n_a, ca_a, h_h_orig, n_h_orig, ca_h_orig in points:
                    dH = (h_h_orig + h_off) - h_a
                    dN = (n_h_orig + n_off) - n_a
                    dCA = (ca_h_orig + ca_off) - ca_a
                    # CSP = sqrt(1/3*(dH^2+(wN*dN)^2+(wCA*dCA)^2)), so CSP < cutoff iff inner < 3*cutoff^2
                    inner = dH * dH + _csp_n_term(dN) + _csp_ca_term(dCA)
                    if inner < 3.0 * cutoff_sq:
                        count_lt += 1
                if count_lt > best_count:
                    best_count = count_lt
                    best_triplets = [(n_off, h_off, ca_off)]
                elif count_lt == best_count:
                    best_triplets.append((n_off, h_off, ca_off))

    # Tie-breaker: pick triplet with smallest |H| + |N| + |CA|
    if best_triplets:
        best_n_offset, best_h_offset, best_ca_offset = min(
            best_triplets,
            key=lambda p: (abs(p[1]) + abs(p[0]) + abs(p[2]), abs(p[0]), abs(p[1]), abs(p[2])),
        )
    else:
        best_n_offset, best_h_offset, best_ca_offset = 0.0, 0.0, 0.0

    return {
        "best_h_offset": best_h_offset,
        "best_n_offset": best_n_offset,
        "best_ca_offset": best_ca_offset,
        "best_count": best_count if best_count >= 0 else 0,
    }


def compute_atom_deltas_with_offset(
    results: List[CSPResult],
    atom: str,
    threshold_config: Optional[Dict] = None,
    grid_params: Optional[Dict] = None,
) -> AtomDeltaResult:
    """
    Compute per-residue 1D atom-specific offsets and significance.

    This treats each nucleus independently using a 1D grid search over the
    configured referencing ranges (from Referencing) and then applies the
    standard outlier-based thresholding on |Δatom|.
    """
    atom_upper = atom.upper()
    if atom_upper not in ("H", "N", "CA"):
        raise ValueError(f"Unsupported atom type for per-atom analysis: {atom}")

    # Collect apo/original holo shift pairs for this atom
    points: List[Tuple[float, float]] = []
    for r in results:
        if atom_upper == "H":
            apo = r.H_apo
            holo_orig = r.H_holo_original
        elif atom_upper == "N":
            apo = r.N_apo
            holo_orig = r.N_holo_original
        else:  # "CA"
            apo = getattr(r, "CA_apo", None)
            holo_orig = getattr(r, "CA_holo_original", None)
        if apo is not None and holo_orig is not None:
            points.append((float(apo), float(holo_orig)))

    # If no valid points, return an empty result
    if not points:
        return AtomDeltaResult(
            atom=atom_upper,
            offset=0.0,
            cutoff=0.0,
            values_abs=[],
            per_residue=[],
        )

    cfg = _Referencing()
    if atom_upper == "H":
        min_offset = (grid_params or {}).get("h_min", cfg.grid_h_min)
        max_offset = (grid_params or {}).get("h_max", cfg.grid_h_max)
        step = (grid_params or {}).get("h_step", cfg.grid_h_step)
    elif atom_upper == "N":
        min_offset = (grid_params or {}).get("n_min", cfg.grid_n_min)
        max_offset = (grid_params or {}).get("n_max", cfg.grid_n_max)
        step = (grid_params or {}).get("n_step", cfg.grid_n_step)
    else:  # "CA"
        min_offset = (grid_params or {}).get("ca_min", cfg.grid_ca_min)
        max_offset = (grid_params or {}).get("ca_max", cfg.grid_ca_max)
        step = (grid_params or {}).get("ca_step", cfg.grid_ca_step)

    gs_cutoff = (grid_params or {}).get("cutoff", cfg.grid_cutoff)

    grid_result = run_offset_grid_search_1d(
        points,
        min_offset=min_offset,
        max_offset=max_offset,
        step=step,
        cutoff=float(gs_cutoff),
    )
    best_offset = float(grid_result["best_offset"])  # type: ignore

    # Compute per-residue deltas using the best offset
    values_abs: List[float] = []
    per_residue: List[Dict[str, object]] = []

    for r in results:
        if atom_upper == "H":
            apo = r.H_apo
            holo_orig = r.H_holo_original
        elif atom_upper == "N":
            apo = r.N_apo
            holo_orig = r.N_holo_original
        else:
            apo = getattr(r, "CA_apo", None)
            holo_orig = getattr(r, "CA_holo_original", None)

        if apo is None or holo_orig is None:
            continue

        delta = (float(holo_orig) + best_offset) - float(apo)
        abs_delta = abs(delta)
        values_abs.append(abs_delta)
        per_residue.append(
            {
                "holo_index": r.holo_index,
                "apo_index": r.apo_index,
                "apo_aa": r.apo_aa,
                "holo_aa": r.holo_aa,
                "delta": delta,
                "delta_abs": abs_delta,
            }
        )

    if not values_abs:
        return AtomDeltaResult(
            atom=atom_upper,
            offset=best_offset,
            cutoff=0.0,
            values_abs=[],
            per_residue=[],
        )

    # Determine thresholds using the same logic as for CSP(A)
    if threshold_config:
        outlier_z = threshold_config.get("outlier_z_score", thresholds.outlier_z_score)
        significance_z = threshold_config.get("significance_z_score", thresholds.significance_z_score)
        max_iterations = threshold_config.get("max_outlier_iterations", thresholds.max_outlier_iterations)
        max_fraction = threshold_config.get("max_outlier_fraction", thresholds.max_outlier_fraction)
        absolute_cutoff = threshold_config.get("absolute_cutoff", thresholds.absolute_cutoff)
    else:
        outlier_z = thresholds.outlier_z_score
        significance_z = thresholds.significance_z_score
        max_iterations = thresholds.max_outlier_iterations
        max_fraction = thresholds.max_outlier_fraction
        absolute_cutoff = thresholds.absolute_cutoff

    if absolute_cutoff is not None:
        threshold_info = _make_absolute_threshold(values_abs, float(absolute_cutoff))
    else:
        threshold_info = compute_threshold_with_outlier_removal(
            values_abs,
            outlier_z,
            significance_z,
            max_iterations,
            max_fraction,
        )

    cutoff = threshold_info.threshold

    for entry in per_residue:
        entry["significant"] = bool(entry["delta_abs"] >= cutoff)

    return AtomDeltaResult(
        atom=atom_upper,
        offset=best_offset,
        cutoff=cutoff,
        values_abs=values_abs,
        per_residue=per_residue,
    )


def _build_param_slug(*, h_min: float, h_max: float, h_step: float, n_min: float, n_max: float, n_step: float, cutoff: float) -> str:
    def fmt(v: float, decimals: int) -> str:
        return (f"{v:.{decimals}f}").rstrip('0').rstrip('.') if '.' in f"{v:.{decimals}f}" else f"{v:.{decimals}f}"
    return (
        f"H_{fmt(h_min, 3)}_{fmt(h_max, 3)}_{fmt(h_step, 3)}__"
        f"N_{fmt(n_min, 3)}_{fmt(n_max, 3)}_{fmt(n_step, 3)}__"
        f"C_{fmt(cutoff, 3)}"
    )


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_grid_heatmap_png(output_png: str, grid_result: Dict[str, object], title: str) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore
    except Exception:
        return
    h_values = np.array(grid_result["h_values"])  # type: ignore
    n_values = np.array(grid_result["n_values"])  # type: ignore
    matrix = grid_result["count_matrix"]  # type: ignore
    data = np.array(matrix, dtype=float)
    best_h = float(grid_result["best_h_offset"])  # type: ignore
    best_n = float(grid_result["best_n_offset"])  # type: ignore

    # Small figure: ~2 inch square heatmap, extra space for horizontal colorbar
    fig, ax = plt.subplots(figsize=(4, 3.5))
    im = ax.imshow(data, origin='lower', aspect='auto', cmap='viridis')
    # Horizontal colorbar below heatmap, pad adds buffer between x-axis label and colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.8, aspect=25)
    cbar.set_label(r'Number of aligned peaks ($\Delta\delta$ < 0.05 ppm)', fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlabel('¹⁵N Offset', fontsize=10)
    ax.set_ylabel('¹H Offset', fontsize=10)
    ax.set_title(title, fontsize=11)
    # Mark optimal offset with black star
    best_h_idx = int(np.argmin(np.abs(h_values - best_h)))
    best_n_idx = int(np.argmin(np.abs(n_values - best_n)))
    ax.scatter(best_n_idx, best_h_idx, marker='*', s=150, c='black', zorder=5, edgecolors='white', linewidths=0.5)
    # ~5 tick marks including 0 (linspace over symmetric range yields 0 at center)
    num_ticks = 5
    n_idx = np.linspace(0, len(n_values) - 1, num=num_ticks).astype(int)
    h_idx = np.linspace(0, len(h_values) - 1, num=num_ticks).astype(int)
    ax.set_xticks(n_idx)
    ax.set_xticklabels([f"{n_values[i]:.2f}" for i in n_idx], fontsize=8, rotation=45, ha='right')
    ax.set_yticks(h_idx)
    ax.set_yticklabels([f"{h_values[i]:.2f}" for i in h_idx], fontsize=8)
    # Add optimal offset result as text
    ax.text(0.02, 0.98, f"(¹⁵N offset, ¹H offset) = ({best_n:.2f}, {best_h:.2f})",
            transform=ax.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout(pad=1.5)
    _ensure_dir(os.path.dirname(output_png))
    fig.savefig(output_png, dpi=200, bbox_inches='tight')
    plt.close(fig)


def save_grid_csv(output_csv: str, grid_result: Dict[str, object], *,
                  h_min: float, h_max: float, h_step: float,
                  n_min: float, n_max: float, n_step: float,
                  cutoff: float) -> None:
    _ensure_dir(os.path.dirname(output_csv))
    h_values: List[float] = grid_result["h_values"]  # type: ignore
    n_values: List[float] = grid_result["n_values"]  # type: ignore
    matrix: List[List[int]] = grid_result["count_matrix"]  # type: ignore
    best_h: float = grid_result["best_h_offset"]  # type: ignore
    best_n: float = grid_result["best_n_offset"]  # type: ignore
    best_count: int = grid_result["best_count"]  # type: ignore

    with open(output_csv, 'w') as f:
        f.write(f"# h_min,{h_min}\n")
        f.write(f"# h_max,{h_max}\n")
        f.write(f"# h_step,{h_step}\n")
        f.write(f"# n_min,{n_min}\n")
        f.write(f"# n_max,{n_max}\n")
        f.write(f"# n_step,{n_step}\n")
        f.write(f"# cutoff,{cutoff}\n")
        f.write("h_values," + ",".join(str(v) for v in h_values) + "\n")
        f.write("n_values," + ",".join(str(v) for v in n_values) + "\n")
        # Matrix rows labelled by h value; columns correspond to n_values order
        for i, h in enumerate(h_values):
            row = matrix[i]
            f.write(str(h) + "," + ",".join(str(c) for c in row) + "\n")
        f.write(f"best_n_offset,{best_n}\n")
        f.write(f"best_h_offset,{best_h}\n")
        f.write(f"best_count,{best_count}\n")


def load_grid_csv(input_csv: str) -> Optional[Dict[str, object]]:
    if not os.path.exists(input_csv):
        return None
    try:
        h_values: List[float] = []
        n_values: List[float] = []
        matrix: List[List[int]] = []
        best_h = 0.0
        best_n = 0.0
        best_count = 0
        stage = 'meta'
        with open(input_csv, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                parts = line.split(',')
                if parts[0] == 'h_values':
                    h_values = [float(x) for x in parts[1:] if x]
                    stage = 'h_values'
                    continue
                if parts[0] == 'n_values':
                    n_values = [float(x) for x in parts[1:] if x]
                    stage = 'matrix'
                    continue
                if parts[0] == 'best_n_offset':
                    best_n = float(parts[1])
                    continue
                if parts[0] == 'best_h_offset':
                    best_h = float(parts[1])
                    continue
                if parts[0] == 'best_count':
                    best_count = int(float(parts[1]))
                    continue
                # Matrix rows: first value is h label
                if stage == 'matrix':
                    # ignore first column (h label)
                    row_counts = [int(float(x)) for x in parts[1:] if x]
                    matrix.append(row_counts)
        if not h_values or not n_values or not matrix:
            return None
        return {
            "h_values": h_values,
            "n_values": n_values,
            "count_matrix": matrix,
            "best_h_offset": best_h,
            "best_n_offset": best_n,
            "best_count": best_count,
        }
    except Exception:
        return None

def compute_shift_reference_offsets(
    results: List[CSPResult],
    threshold: float
) -> Tuple[float, float]:
    """
    Compute optimal offsets for chemical shift referencing using insignificant CSPs.
    
    Parameters
    ----------
    results : List[CSPResult]
        List of CSP results with significance already determined
    threshold : float
        Significance threshold for identifying reference subset
        
    Returns
    -------
    Tuple[float, float]
        (N_offset, H_offset) - optimal offsets to apply to holo shifts
    """
    # Identify insignificant CSPs (below threshold) for referencing
    reference_subset = []
    for r in results:
        if (r.csp_A is not None and r.csp_A < threshold and 
            r.H_apo is not None and r.N_apo is not None and 
            r.H_holo is not None and r.N_holo is not None):
            reference_subset.append(r)
    
    if not reference_subset:
        # No reference subset available, return zero offsets
        return 0.0, 0.0
    
    # Calculate optimal offsets using mean difference for reference subset
    # For shift type S: optimal_offset = mean(S_apo - S_holo)
    n_differences = []
    h_differences = []
    
    for r in reference_subset:
        n_differences.append(r.N_apo - r.N_holo)
        h_differences.append(r.H_apo - r.H_holo)
    
    n_offset = statistics.mean(n_differences) if n_differences else 0.0
    h_offset = statistics.mean(h_differences) if h_differences else 0.0
    
    return n_offset, h_offset


def apply_shift_offsets(
    H_holo: Dict[int, float],
    N_holo: Dict[int, float],
    h_offset: float,
    n_offset: float
) -> Tuple[Dict[int, float], Dict[int, float]]:
    """
    Apply referencing offsets to holo chemical shifts.
    
    Parameters
    ----------
    H_holo : Dict[int, float]
        Original holo H shifts
    N_holo : Dict[int, float]
        Original holo N shifts
    h_offset : float
        Offset to apply to H shifts
    n_offset : float
        Offset to apply to N shifts
        
    Returns
    -------
    Tuple[Dict[int, float], Dict[int, float]]
        (H_holo_referenced, N_holo_referenced) - referenced shift dictionaries
    """
    H_holo_referenced = {pos: shift + h_offset for pos, shift in H_holo.items()}
    N_holo_referenced = {pos: shift + n_offset for pos, shift in N_holo.items()}
    
    return H_holo_referenced, N_holo_referenced


def compute_csp_A(
    mapping: List[Tuple[int, int]],
    apo_seq: str,
    holo_seq: str,
    H_apo: Dict[int, float],
    N_apo: Dict[int, float],
    H_holo: Dict[int, float],
    N_holo: Dict[int, float],
    threshold_config: Optional[Dict] = None,
    enable_referencing: bool = True,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
) -> List[CSPResult]:
    results: List[CSPResult] = []
    values: List[float] = []
    # Debug verbosity check
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))
    
    # Store original shifts for referencing
    H_holo_original = H_holo.copy()
    N_holo_original = N_holo.copy()
    
    for i_apo, i_holo in mapping:
        aa_apo = apo_seq[i_apo - 1] if 0 < i_apo <= len(apo_seq) else "?"
        aa_holo = holo_seq[i_holo - 1] if 0 < i_holo <= len(holo_seq) else "?"
        h_a = H_apo.get(i_apo)
        n_a = N_apo.get(i_apo)
        h_h = H_holo.get(i_holo)
        n_h = N_holo.get(i_holo)
        dH = (h_h - h_a) if (h_h is not None and h_a is not None) else None
        dN = (n_h - n_a) if (n_h is not None and n_a is not None) else None
        csp_val: Optional[float] = None
        if dH is not None and dN is not None:
            csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
            values.append(csp_val)
        if _verbose and (dH is not None or dN is not None):
            _dH = f"{dH:.3f}" if dH is not None else "NA"
            _dN = f"{dN:.3f}" if dN is not None else "NA"
            _A = f"{csp_val:.3f}" if csp_val is not None else "NA"
            print(f"[CSP] apo {i_apo}{aa_apo} ↔ holo {i_holo}{aa_holo} : dH={_dH} dN={_dN} A={_A}")
        results.append(
            CSPResult(
                apo_index=i_apo,
                holo_index=i_holo,
                apo_aa=aa_apo,
                holo_aa=aa_holo,
                H_apo=h_a,
                N_apo=n_a,
                H_holo=h_h,
                N_holo=n_h,
                dH=dH,
                dN=dN,
                csp_A=csp_val,
                significant=None,
                significant_1sd=None,
                significant_2sd=None,
                H_holo_original=H_holo_original.get(i_holo),
                N_holo_original=N_holo_original.get(i_holo),
                H_offset=None,
                N_offset=None,
            )
        )

    threshold_info: Optional[ThresholdComputation] = None
    threshold_info_ref: Optional[ThresholdComputation] = None
    final_stats: Optional[ThresholdComputation] = None
    values_referenced: List[float] = []

    # determine significance
    if values:
        # Use dynamic threshold config if provided, otherwise use global config
        if threshold_config:
            outlier_z = threshold_config.get('outlier_z_score', thresholds.outlier_z_score)
            significance_z = threshold_config.get('significance_z_score', thresholds.significance_z_score)
            max_iterations = threshold_config.get('max_outlier_iterations', thresholds.max_outlier_iterations)
            max_fraction = threshold_config.get('max_outlier_fraction', thresholds.max_outlier_fraction)
            absolute_cutoff = threshold_config.get('absolute_cutoff', thresholds.absolute_cutoff)
        else:
            outlier_z = thresholds.outlier_z_score
            significance_z = thresholds.significance_z_score
            max_iterations = thresholds.max_outlier_iterations
            max_fraction = thresholds.max_outlier_fraction
            absolute_cutoff = thresholds.absolute_cutoff

        if absolute_cutoff is not None:
            threshold_info = _make_absolute_threshold(values, float(absolute_cutoff))
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP] Using absolute cutoff: {cutoff:.3f}")
        else:
            threshold_info = compute_threshold_with_outlier_removal(
                values,
                outlier_z,
                significance_z,
                max_iterations,
                max_fraction
            )
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP] values={len(values)} mean={threshold_info.mean:.3f} sd={threshold_info.sd:.3f}")
                print(f"[CSP] outlier removal: {threshold_info.iterations} iterations, {threshold_info.outliers_removed} outliers removed")
                print(f"[CSP] final threshold: {cutoff:.3f}")

        final_stats = threshold_info

        # Apply referencing if enabled
        if enable_referencing:
            if _verbose:
                print("[CSP] Computing referencing offsets")

            method = (referencing_method or getattr(_Referencing, 'method', 'mean'))
            n_offset = 0.0
            h_offset = 0.0

            if method == 'grid':
                # Resolve grid parameters (config defaults overridden by grid_params)
                cfg = _Referencing()
                h_min = (grid_params or {}).get('h_min', cfg.grid_h_min)
                h_max = (grid_params or {}).get('h_max', cfg.grid_h_max)
                h_step = (grid_params or {}).get('h_step', cfg.grid_h_step)
                n_min = (grid_params or {}).get('n_min', cfg.grid_n_min)
                n_max = (grid_params or {}).get('n_max', cfg.grid_n_max)
                n_step = (grid_params or {}).get('n_step', cfg.grid_n_step)
                gs_cutoff = (grid_params or {}).get('cutoff', cfg.grid_cutoff)

                if _verbose:
                    print(f"[CSP] Grid search over H[{h_min},{h_max},{h_step}] and N[{n_min},{n_max},{n_step}] with cutoff={gs_cutoff}")
                # Determine cache paths
                grid_result: Dict[str, object]
                slug = _build_param_slug(h_min=h_min, h_max=h_max, h_step=h_step, n_min=n_min, n_max=n_max, n_step=n_step, cutoff=float(gs_cutoff))
                out_dir = os.path.join(paths.outputs_dir, target_id) if target_id else None
                csv_path = os.path.join(out_dir, f"offset_grid_{slug}.csv") if out_dir else None
                png_path = os.path.join(out_dir, f"offset_grid_{slug}.png") if out_dir else None
                use_cache = bool(getattr(cfg, 'cache_results', True) and target_id and out_dir)

                loaded = False
                if use_cache and csv_path and os.path.exists(csv_path):
                    loaded_result = load_grid_csv(csv_path)
                    if loaded_result is not None:
                        grid_result = loaded_result
                        loaded = True
                        if _verbose:
                            print(f"[CSP] Loaded grid result from cache: {csv_path}")
                if not loaded:
                    grid_result = run_offset_grid_search(
                        results,
                        h_min=h_min,
                        h_max=h_max,
                        h_step=h_step,
                        n_min=n_min,
                        n_max=n_max,
                        n_step=n_step,
                        cutoff=float(gs_cutoff),
                    )
                    if use_cache and csv_path:
                        save_grid_csv(
                            csv_path,
                            grid_result,
                            h_min=h_min, h_max=h_max, h_step=h_step,
                            n_min=n_min, n_max=n_max, n_step=n_step,
                            cutoff=float(gs_cutoff),
                        )
                        if _verbose:
                            print(f"[CSP] Saved grid result to: {csv_path}")
                # Save heatmap if desired
                cfg_save = bool(getattr(cfg, 'save_heatmap', True))
                if cfg_save and png_path:
                    title = f"Grid Search for Optimal ¹H/¹⁵N Offsets\n(PDB {target_id.upper()})" if target_id else "Grid Search for Optimal ¹H/¹⁵N Offsets"
                    save_grid_heatmap_png(png_path, grid_result, title)
                h_offset = float(grid_result["best_h_offset"])  # type: ignore
                n_offset = float(grid_result["best_n_offset"])  # type: ignore
                if _verbose:
                    print(f"[CSP] Grid best offsets: N_offset={n_offset:.3f}, H_offset={h_offset:.3f} (count={grid_result['best_count']})")
            else:
                if _verbose:
                    print("[CSP] Mean-based referencing using insignificant CSPs")
                n_offset, h_offset = compute_shift_reference_offsets(results, cutoff)
            if n_offset != 0.0 or h_offset != 0.0:
                if _verbose:
                    print(f"[CSP] Applying offsets: N_offset={n_offset:.3f}, H_offset={h_offset:.3f}")
                
                # Apply offsets to holo shifts
                H_holo_referenced, N_holo_referenced = apply_shift_offsets(
                    H_holo_original, N_holo_original, h_offset, n_offset
                )
                
                # Recalculate CSPs with referenced shifts
                values_referenced = []
                for r in results:
                    if r.H_holo_original is not None and r.N_holo_original is not None:
                        h_h_ref = H_holo_referenced.get(r.holo_index)
                        n_h_ref = N_holo_referenced.get(r.holo_index)
                        
                        if h_h_ref is not None and n_h_ref is not None:
                            dH_ref = h_h_ref - r.H_apo if r.H_apo is not None else None
                            dN_ref = n_h_ref - r.N_apo if r.N_apo is not None else None
                            
                            if dH_ref is not None and dN_ref is not None:
                                csp_val_ref = math.sqrt(0.5 * (dH_ref * dH_ref + _csp_n_term(dN_ref)))
                                values_referenced.append(csp_val_ref)
                                
                                # Update result with referenced values
                                r.H_holo = h_h_ref
                                r.N_holo = n_h_ref
                                r.dH = dH_ref
                                r.dN = dN_ref
                                r.csp_A = csp_val_ref
                                r.H_offset = h_offset
                                r.N_offset = n_offset
                
                # Recalculate significance threshold with referenced CSPs
                if values_referenced:
                    if absolute_cutoff is not None:
                        threshold_info_ref = _make_absolute_threshold(values_referenced, float(absolute_cutoff))
                    else:
                        threshold_info_ref = compute_threshold_with_outlier_removal(
                            values_referenced,
                            outlier_z,
                            significance_z,
                            max_iterations,
                            max_fraction
                        )
                    if threshold_info_ref:
                        cutoff = threshold_info_ref.threshold
                        final_stats = threshold_info_ref
                        if _verbose:
                            print(f"[CSP] Referenced CSPs: {len(values_referenced)} mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}")
                            print(f"[CSP] Referenced threshold: {cutoff:.3f}")

    threshold_1sd: Optional[float] = None
    threshold_2sd: Optional[float] = None

    if final_stats and final_stats.cleaned_values:
        threshold_1sd = final_stats.mean + 1.0 * final_stats.sd
        threshold_2sd = final_stats.mean + 2.0 * final_stats.sd
        if _verbose:
            print(f"[CSP] Additional thresholds: 1SD={threshold_1sd:.3f}, 2SD={threshold_2sd:.3f}")

    for r in results:
        r.significant = (r.csp_A is not None and final_stats is not None and r.csp_A >= cutoff) if values else None
        if final_stats and r.csp_A is not None:
            if final_stats.sd > 0.0:
                r.z_score = (r.csp_A - final_stats.mean) / final_stats.sd
            else:
                r.z_score = 0.0
        else:
            r.z_score = None

        if threshold_1sd is not None and r.csp_A is not None:
            r.significant_1sd = (r.csp_A >= threshold_1sd)
        else:
            r.significant_1sd = None

        if threshold_2sd is not None and r.csp_A is not None:
            r.significant_2sd = (r.csp_A >= threshold_2sd)
        else:
            r.significant_2sd = None

    return results


def compute_csp_multiple_saveframes(
    apo_sequences: List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]],
    holo_sequences: List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]],
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    threshold_config: Optional[Dict] = None,
    *,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CSPs for all possible apo-holo sequence pairs from multiple saveframes.
    
    Only aligns sequences that have both H and N shifts.
    Returns CSPs for the alignment with the best alignment score.
    """
    from .align import align_global
    
    # Debug verbosity check
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))
    
    if _verbose:
        print(f"[CSP] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")
    
    # First pass: find the alignment with the best score
    best_alignment_score = float('-inf')
    best_alignment_info = None
    
    for apo_seq, H_apo, N_apo, _CA_apo, apo_saveframe in apo_sequences:
        for holo_seq, H_holo, N_holo, _CA_holo, holo_saveframe in holo_sequences:
            if _verbose:
                print(f"[CSP] Aligning apo({apo_saveframe}) vs holo({holo_saveframe})")
            
            # Perform global alignment
            aligned_apo, aligned_holo, mapping, alignment_score = align_global(apo_seq, holo_seq)
            
            if _verbose:
                print(f"[CSP] Alignment score: {alignment_score}, mapped pairs: {len(mapping)}")
            
            # Keep the alignment with the best score
            if alignment_score > best_alignment_score:
                best_alignment_score = alignment_score
                best_alignment_info = (aligned_apo, aligned_holo, H_apo, N_apo, H_holo, N_holo, apo_saveframe, holo_saveframe)
    
    # Second pass: compute CSPs only for the best alignment
    if best_alignment_info is None:
        return []
    
    aligned_apo, aligned_holo, H_apo, N_apo, H_holo, N_holo, apo_saveframe, holo_saveframe = best_alignment_info
    
    if _verbose:
        print(f"[CSP] Computing CSPs for best alignment: apo({apo_saveframe}) vs holo({holo_saveframe}) with score {best_alignment_score}")
    
    # Compute CSPs for the best alignment
    results = compute_csp_from_aligned_sequences(
        aligned_apo,
        aligned_holo,
        H_apo,
        N_apo,
        H_holo,
        N_holo,
        threshold_config,
        enable_referencing=True,
        referencing_method=referencing_method,
        grid_params=grid_params,
        target_id=(target_id or holo_pdb),
    )
    
    if _verbose:
        valid_csps = sum(1 for r in results if r.csp_A is not None)
        print(f"[CSP] Valid CSPs: {valid_csps}")
    
    return results


def compute_csp_from_aligned_sequences(
    aligned_apo: str,
    aligned_holo: str,
    H_apo: Dict[int, float],
    N_apo: Dict[int, float],
    H_holo: Dict[int, float],
    N_holo: Dict[int, float],
    threshold_config: Optional[Dict] = None,
    enable_referencing: bool = True,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CSPs from aligned sequences, properly handling gaps.
    
    This function takes the aligned sequences and maps chemical shifts
    according to the alignment, accounting for gaps.
    """
    results: List[CSPResult] = []
    values: List[float] = []
    
    # Debug verbosity check
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))
    
    # Store original shifts for referencing
    H_holo_original = H_holo.copy()
    N_holo_original = N_holo.copy()
    
    # Track original sequence positions
    apo_pos = 0
    holo_pos = 0
    
    for i, (aa_apo, aa_holo) in enumerate(zip(aligned_apo, aligned_holo)):
        # Update positions based on gaps
        if aa_apo != "-":
            apo_pos += 1
        if aa_holo != "-":
            holo_pos += 1
            
        # Only compute CSP for positions where both sequences have amino acids (no gaps)
        if aa_apo != "-" and aa_holo != "-":
            # Get chemical shifts for current positions
            h_a = H_apo.get(apo_pos)
            n_a = N_apo.get(apo_pos)
            h_h = H_holo.get(holo_pos)
            n_h = N_holo.get(holo_pos)
            
            # Calculate differences
            dH = (h_h - h_a) if (h_h is not None and h_a is not None) else None
            dN = (n_h - n_a) if (n_h is not None and n_a is not None) else None
            
            # Calculate CSP
            csp_val: Optional[float] = None
            if dH is not None and dN is not None:
                csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
                values.append(csp_val)
            
            if _verbose and (dH is not None or dN is not None):
                _dH = f"{dH:.3f}" if dH is not None else "NA"
                _dN = f"{dN:.3f}" if dN is not None else "NA"
                _A = f"{csp_val:.3f}" if csp_val is not None else "NA"
                print(f"[CSP] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : dH={_dH} dN={_dN} A={_A}")
            
            results.append(
                CSPResult(
                    apo_index=apo_pos,
                    holo_index=holo_pos,
                    apo_aa=aa_apo,
                    holo_aa=aa_holo,
                    H_apo=h_a,
                    N_apo=n_a,
                    H_holo=h_h,
                    N_holo=n_h,
                    dH=dH,
                    dN=dN,
                    csp_A=csp_val,
                    significant=None,
                    significant_1sd=None,
                    significant_2sd=None,
                    H_holo_original=H_holo_original.get(holo_pos),
                    N_holo_original=N_holo_original.get(holo_pos),
                    H_offset=None,
                    N_offset=None,
                )
            )

    threshold_info: Optional[ThresholdComputation] = None
    threshold_info_ref: Optional[ThresholdComputation] = None
    final_stats: Optional[ThresholdComputation] = None
    values_referenced: List[float] = []

    # Determine significance
    if values:
        # Use dynamic threshold config if provided, otherwise use global config
        if threshold_config:
            outlier_z = threshold_config.get('outlier_z_score', thresholds.outlier_z_score)
            significance_z = threshold_config.get('significance_z_score', thresholds.significance_z_score)
            max_iterations = threshold_config.get('max_outlier_iterations', thresholds.max_outlier_iterations)
            max_fraction = threshold_config.get('max_outlier_fraction', thresholds.max_outlier_fraction)
            absolute_cutoff = threshold_config.get('absolute_cutoff', thresholds.absolute_cutoff)
        else:
            outlier_z = thresholds.outlier_z_score
            significance_z = thresholds.significance_z_score
            max_iterations = thresholds.max_outlier_iterations
            max_fraction = thresholds.max_outlier_fraction
            absolute_cutoff = thresholds.absolute_cutoff

        if absolute_cutoff is not None:
            threshold_info = _make_absolute_threshold(values, float(absolute_cutoff))
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP] Using absolute cutoff: {cutoff:.3f}")
        else:
            threshold_info = compute_threshold_with_outlier_removal(
                values,
                outlier_z,
                significance_z,
                max_iterations,
                max_fraction
            )
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP] values={len(values)} mean={threshold_info.mean:.3f} sd={threshold_info.sd:.3f}")
                print(f"[CSP] outlier removal: {threshold_info.iterations} iterations, {threshold_info.outliers_removed} outliers removed")
                print(f"[CSP] final threshold: {cutoff:.3f}")

        final_stats = threshold_info

        # Apply referencing if enabled
        if enable_referencing:
            if _verbose:
                print("[CSP] Computing referencing offsets")

            method = (referencing_method or getattr(_Referencing, 'method', 'mean'))
            n_offset = 0.0
            h_offset = 0.0

            if method == 'grid':
                # Resolve grid parameters (config defaults overridden by grid_params)
                cfg = _Referencing()
                h_min = (grid_params or {}).get('h_min', cfg.grid_h_min)
                h_max = (grid_params or {}).get('h_max', cfg.grid_h_max)
                h_step = (grid_params or {}).get('h_step', cfg.grid_h_step)
                n_min = (grid_params or {}).get('n_min', cfg.grid_n_min)
                n_max = (grid_params or {}).get('n_max', cfg.grid_n_max)
                n_step = (grid_params or {}).get('n_step', cfg.grid_n_step)
                gs_cutoff = (grid_params or {}).get('cutoff', cfg.grid_cutoff)

                if _verbose:
                    print(f"[CSP] Grid search over H[{h_min},{h_max},{h_step}] and N[{n_min},{n_max},{n_step}] with cutoff={gs_cutoff}")
                # Determine cache paths
                grid_result: Dict[str, object]
                slug = _build_param_slug(h_min=h_min, h_max=h_max, h_step=h_step, n_min=n_min, n_max=n_max, n_step=n_step, cutoff=float(gs_cutoff))
                out_dir = os.path.join(paths.outputs_dir, target_id) if target_id else None
                csv_path = os.path.join(out_dir, f"offset_grid_{slug}.csv") if out_dir else None
                png_path = os.path.join(out_dir, f"offset_grid_{slug}.png") if out_dir else None
                use_cache = bool(getattr(cfg, 'cache_results', True) and target_id and out_dir)

                loaded = False
                if use_cache and csv_path and os.path.exists(csv_path):
                    loaded_result = load_grid_csv(csv_path)
                    if loaded_result is not None:
                        grid_result = loaded_result
                        loaded = True
                        if _verbose:
                            print(f"[CSP] Loaded grid result from cache: {csv_path}")
                if not loaded:
                    grid_result = run_offset_grid_search(
                        results,
                        h_min=h_min,
                        h_max=h_max,
                        h_step=h_step,
                        n_min=n_min,
                        n_max=n_max,
                        n_step=n_step,
                        cutoff=float(gs_cutoff),
                    )
                    if use_cache and csv_path:
                        save_grid_csv(
                            csv_path,
                            grid_result,
                            h_min=h_min, h_max=h_max, h_step=h_step,
                            n_min=n_min, n_max=n_max, n_step=n_step,
                            cutoff=float(gs_cutoff),
                        )
                        if _verbose:
                            print(f"[CSP] Saved grid result to: {csv_path}")
                # Save heatmap if desired
                cfg_save = bool(getattr(cfg, 'save_heatmap', True))
                if cfg_save and png_path:
                    title = f"Grid Search for Optimal ¹H/¹⁵N Offsets\n(PDB {target_id.upper()})" if target_id else "Grid Search for Optimal ¹H/¹⁵N Offsets"
                    save_grid_heatmap_png(png_path, grid_result, title)
                h_offset = float(grid_result["best_h_offset"])  # type: ignore
                n_offset = float(grid_result["best_n_offset"])  # type: ignore
                if _verbose:
                    print(f"[CSP] Grid best offsets: N_offset={n_offset:.3f}, H_offset={h_offset:.3f} (count={grid_result['best_count']})")
            else:
                if _verbose:
                    print("[CSP] Mean-based referencing using insignificant CSPs")
                n_offset, h_offset = compute_shift_reference_offsets(results, cutoff)
            if n_offset != 0.0 or h_offset != 0.0:
                if _verbose:
                    print(f"[CSP] Applying offsets: N_offset={n_offset:.3f}, H_offset={h_offset:.3f}")
                
                # Apply offsets to holo shifts
                H_holo_referenced, N_holo_referenced = apply_shift_offsets(
                    H_holo_original, N_holo_original, h_offset, n_offset
                )
                
                # Recalculate CSPs with referenced shifts
                values_referenced = []
                for r in results:
                    if r.H_holo_original is not None and r.N_holo_original is not None:
                        h_h_ref = H_holo_referenced.get(r.holo_index)
                        n_h_ref = N_holo_referenced.get(r.holo_index)
                        
                        if h_h_ref is not None and n_h_ref is not None:
                            dH_ref = h_h_ref - r.H_apo if r.H_apo is not None else None
                            dN_ref = n_h_ref - r.N_apo if r.N_apo is not None else None
                            
                            if dH_ref is not None and dN_ref is not None:
                                csp_val_ref = math.sqrt(0.5 * (dH_ref * dH_ref + _csp_n_term(dN_ref)))
                                values_referenced.append(csp_val_ref)
                                
                                # Update result with referenced values
                                r.H_holo = h_h_ref
                                r.N_holo = n_h_ref
                                r.dH = dH_ref
                                r.dN = dN_ref
                                r.csp_A = csp_val_ref
                                r.H_offset = h_offset
                                r.N_offset = n_offset
                
                # Recalculate significance threshold with referenced CSPs
                if values_referenced:
                    if absolute_cutoff is not None:
                        threshold_info_ref = _make_absolute_threshold(values_referenced, float(absolute_cutoff))
                    else:
                        threshold_info_ref = compute_threshold_with_outlier_removal(
                            values_referenced,
                            outlier_z,
                            significance_z,
                            max_iterations,
                            max_fraction
                        )
                    if threshold_info_ref:
                        cutoff = threshold_info_ref.threshold
                        final_stats = threshold_info_ref
                        if _verbose:
                            print(f"[CSP] Referenced CSPs: {len(values_referenced)} mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}")
                            print(f"[CSP] Referenced threshold: {cutoff:.3f}")

    threshold_1sd: Optional[float] = None
    threshold_2sd: Optional[float] = None

    if final_stats and final_stats.cleaned_values:
        threshold_1sd = final_stats.mean + 1.0 * final_stats.sd
        threshold_2sd = final_stats.mean + 2.0 * final_stats.sd
        if _verbose:
            print(f"[CSP] Additional thresholds: 1SD={threshold_1sd:.3f}, 2SD={threshold_2sd:.3f}")

    for r in results:
        r.significant = (r.csp_A is not None and final_stats is not None and r.csp_A >= cutoff) if values else None
        if final_stats and r.csp_A is not None:
            if final_stats.sd > 0.0:
                r.z_score = (r.csp_A - final_stats.mean) / final_stats.sd
            else:
                r.z_score = 0.0
        else:
            r.z_score = None

        if threshold_1sd is not None and r.csp_A is not None:
            r.significant_1sd = (r.csp_A >= threshold_1sd)
        else:
            r.significant_1sd = None

        if threshold_2sd is not None and r.csp_A is not None:
            r.significant_2sd = (r.csp_A >= threshold_2sd)
        else:
            r.significant_2sd = None
    
    return results


def compute_csp_from_aligned_sequences_ca(
    aligned_apo: str,
    aligned_holo: str,
    H_apo: Dict[int, float],
    N_apo: Dict[int, float],
    CA_apo: Dict[int, float],
    H_holo: Dict[int, float],
    N_holo: Dict[int, float],
    CA_holo: Dict[int, float],
    threshold_config: Optional[Dict] = None,
    enable_referencing: bool = True,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CA-inclusive CSPs from aligned sequences, handling gaps.

    Uses the CA-inclusive CSP formula:
        CSP = sqrt(1/3*(dH^2+(wN*dN)^2+(wCA*dCA)^2)), wN/wCA from config Compute
    """
    results: List[CSPResult] = []
    values: List[float] = []

    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))

    # Store original shifts for referencing
    H_holo_original = H_holo.copy()
    N_holo_original = N_holo.copy()
    CA_holo_original = CA_holo.copy()

    # Track original sequence positions
    apo_pos = 0
    holo_pos = 0

    # Points used for 3D grid search:
    # (H_apo, N_apo, CA_apo, H_holo_orig, N_holo_orig, CA_holo_orig)
    grid_points: List[Tuple[float, float, float, float, float, float]] = []

    for aa_apo, aa_holo in zip(aligned_apo, aligned_holo):
        if aa_apo != "-":
            apo_pos += 1
        if aa_holo != "-":
            holo_pos += 1

        if aa_apo != "-" and aa_holo != "-":
            h_a = H_apo.get(apo_pos)
            n_a = N_apo.get(apo_pos)
            ca_a = CA_apo.get(apo_pos)
            h_h = H_holo.get(holo_pos)
            n_h = N_holo.get(holo_pos)
            ca_h = CA_holo.get(holo_pos)

            dH = (h_h - h_a) if (h_h is not None and h_a is not None) else None
            dN = (n_h - n_a) if (n_h is not None and n_a is not None) else None
            dCA = (ca_h - ca_a) if (ca_h is not None and ca_a is not None) else None

            csp_val: Optional[float] = None
            if dH is not None and dN is not None and dCA is not None:
                csp_val = math.sqrt(
                    (1.0 / 3.0) * (dH * dH + _csp_n_term(dN) + _csp_ca_term(dCA))
                )
                values.append(csp_val)
                grid_points.append(
                    (
                        float(h_a),
                        float(n_a),
                        float(ca_a),
                        float(H_holo_original.get(holo_pos, h_h)),
                        float(N_holo_original.get(holo_pos, n_h)),
                        float(CA_holo_original.get(holo_pos, ca_h)),
                    )
                )

            if _verbose and (dH is not None or dN is not None or dCA is not None):
                _dH = f"{dH:.3f}" if dH is not None else "NA"
                _dN = f"{dN:.3f}" if dN is not None else "NA"
                _dCA = f"{dCA:.3f}" if dCA is not None else "NA"
                _C = f"{csp_val:.3f}" if csp_val is not None else "NA"
                print(f"[CSP-CA] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : dH={_dH} dN={_dN} dCA={_dCA} C={_C}")

            results.append(
                CSPResult(
                    apo_index=apo_pos,
                    holo_index=holo_pos,
                    apo_aa=aa_apo,
                    holo_aa=aa_holo,
                    H_apo=h_a,
                    N_apo=n_a,
                    H_holo=h_h,
                    N_holo=n_h,
                    dH=dH,
                    dN=dN,
                    csp_A=csp_val,
                    significant=None,
                    significant_1sd=None,
                    significant_2sd=None,
                    H_holo_original=H_holo_original.get(holo_pos),
                    N_holo_original=N_holo_original.get(holo_pos),
                    H_offset=None,
                    N_offset=None,
                    CA_apo=ca_a,
                    CA_holo=ca_h,
                    CA_holo_original=CA_holo_original.get(holo_pos),
                    CA_offset=None,
                    dCA=dCA,
                )
            )

    threshold_info: Optional[ThresholdComputation] = None
    threshold_info_ref: Optional[ThresholdComputation] = None
    final_stats: Optional[ThresholdComputation] = None
    values_referenced: List[float] = []

    if values:
        if threshold_config:
            outlier_z = threshold_config.get('outlier_z_score', thresholds.outlier_z_score)
            significance_z = threshold_config.get('significance_z_score', thresholds.significance_z_score)
            max_iterations = threshold_config.get('max_outlier_iterations', thresholds.max_outlier_iterations)
            max_fraction = threshold_config.get('max_outlier_fraction', thresholds.max_outlier_fraction)
            absolute_cutoff = threshold_config.get('absolute_cutoff', thresholds.absolute_cutoff)
        else:
            outlier_z = thresholds.outlier_z_score
            significance_z = thresholds.significance_z_score
            max_iterations = thresholds.max_outlier_iterations
            max_fraction = thresholds.max_outlier_fraction
            absolute_cutoff = thresholds.absolute_cutoff

        if absolute_cutoff is not None:
            threshold_info = _make_absolute_threshold(values, float(absolute_cutoff))
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP-CA] Using absolute cutoff: {cutoff:.3f}")
        else:
            threshold_info = compute_threshold_with_outlier_removal(
                values,
                outlier_z,
                significance_z,
                max_iterations,
                max_fraction,
            )
            cutoff = threshold_info.threshold
            if _verbose:
                print(f"[CSP-CA] values={len(values)} mean={threshold_info.mean:.3f} sd={threshold_info.sd:.3f}")
                print(f"[CSP-CA] outlier removal: {threshold_info.iterations} iterations, {threshold_info.outliers_removed} outliers removed")
                print(f"[CSP-CA] final threshold: {cutoff:.3f}")

        final_stats = threshold_info

        if enable_referencing and grid_points:
            if _verbose:
                print("[CSP-CA] Computing referencing offsets (3D grid)")

            method = (referencing_method or getattr(_Referencing, 'method', 'grid'))
            n_offset = 0.0
            h_offset = 0.0
            ca_offset = 0.0

            if method == 'grid':
                cfg = _Referencing()
                h_min = (grid_params or {}).get('h_min', cfg.grid_h_min)
                h_max = (grid_params or {}).get('h_max', cfg.grid_h_max)
                h_step = (grid_params or {}).get('h_step', cfg.grid_h_step)
                n_min = (grid_params or {}).get('n_min', cfg.grid_n_min)
                n_max = (grid_params or {}).get('n_max', cfg.grid_n_max)
                n_step = (grid_params or {}).get('n_step', cfg.grid_n_step)
                ca_min = (grid_params or {}).get('ca_min', cfg.grid_ca_min)
                ca_max = (grid_params or {}).get('ca_max', cfg.grid_ca_max)
                ca_step = (grid_params or {}).get('ca_step', cfg.grid_ca_step)
                gs_cutoff = (grid_params or {}).get('cutoff', cfg.grid_cutoff)

                if _verbose:
                    print(
                        f"[CSP-CA] 3D grid over H[{h_min},{h_max},{h_step}] "
                        f"N[{n_min},{n_max},{n_step}] CA[{ca_min},{ca_max},{ca_step}] "
                        f"with cutoff={gs_cutoff}"
                    )

                grid_result = run_offset_grid_search_3d(
                    grid_points,
                    h_min=h_min,
                    h_max=h_max,
                    h_step=h_step,
                    n_min=n_min,
                    n_max=n_max,
                    n_step=n_step,
                    ca_min=ca_min,
                    ca_max=ca_max,
                    ca_step=ca_step,
                    cutoff=float(gs_cutoff),
                )
                h_offset = float(grid_result["best_h_offset"])  # type: ignore
                n_offset = float(grid_result["best_n_offset"])  # type: ignore
                ca_offset = float(grid_result["best_ca_offset"])  # type: ignore
                if _verbose:
                    print(
                        f"[CSP-CA] Grid best offsets: "
                        f"N_offset={n_offset:.3f}, H_offset={h_offset:.3f}, CA_offset={ca_offset:.3f} "
                        f"(count={grid_result['best_count']})"
                    )
            else:
                if _verbose:
                    print("[CSP-CA] Mean-based referencing not implemented for CA; skipping referencing")

            if n_offset != 0.0 or h_offset != 0.0 or ca_offset != 0.0:
                if _verbose:
                    print(
                        f"[CSP-CA] Applying offsets: "
                        f"N_offset={n_offset:.3f}, H_offset={h_offset:.3f}, CA_offset={ca_offset:.3f}"
                    )

                H_holo_referenced, N_holo_referenced = apply_shift_offsets(
                    H_holo_original, N_holo_original, h_offset, n_offset
                )
                CA_holo_referenced = {pos: shift + ca_offset for pos, shift in CA_holo_original.items()}

                values_referenced = []
                for r in results:
                    if (
                        r.H_holo_original is not None
                        and r.N_holo_original is not None
                        and r.CA_holo_original is not None
                    ):
                        h_h_ref = H_holo_referenced.get(r.holo_index)
                        n_h_ref = N_holo_referenced.get(r.holo_index)
                        ca_h_ref = CA_holo_referenced.get(r.holo_index)

                        if (
                            h_h_ref is not None
                            and n_h_ref is not None
                            and ca_h_ref is not None
                            and r.H_apo is not None
                            and r.N_apo is not None
                            and r.CA_apo is not None
                        ):
                            dH_ref = h_h_ref - r.H_apo
                            dN_ref = n_h_ref - r.N_apo
                            dCA_ref = ca_h_ref - r.CA_apo

                            csp_val_ref = math.sqrt(
                                (1.0 / 3.0)
                                * (dH_ref * dH_ref + _csp_n_term(dN_ref) + _csp_ca_term(dCA_ref))
                            )
                            values_referenced.append(csp_val_ref)

                            r.H_holo = h_h_ref
                            r.N_holo = n_h_ref
                            r.dH = dH_ref
                            r.dN = dN_ref
                            r.CA_holo = ca_h_ref
                            r.dCA = dCA_ref
                            r.csp_A = csp_val_ref
                            r.H_offset = h_offset
                            r.N_offset = n_offset
                            r.CA_offset = ca_offset

                if values_referenced:
                    if absolute_cutoff is not None:
                        threshold_info_ref = _make_absolute_threshold(values_referenced, float(absolute_cutoff))
                    else:
                        threshold_info_ref = compute_threshold_with_outlier_removal(
                            values_referenced,
                            outlier_z,
                            significance_z,
                            max_iterations,
                            max_fraction,
                        )
                    if threshold_info_ref:
                        cutoff = threshold_info_ref.threshold
                        final_stats = threshold_info_ref
                        if _verbose:
                            print(
                                f"[CSP-CA] Referenced CSPs: {len(values_referenced)} "
                                f"mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}"
                            )
                            print(f"[CSP-CA] Referenced threshold: {cutoff:.3f}")

    threshold_1sd: Optional[float] = None
    threshold_2sd: Optional[float] = None

    if final_stats and final_stats.cleaned_values:
        threshold_1sd = final_stats.mean + 1.0 * final_stats.sd
        threshold_2sd = final_stats.mean + 2.0 * final_stats.sd
        if _verbose:
            print(f"[CSP-CA] Additional thresholds: 1SD={threshold_1sd:.3f}, 2SD={threshold_2sd:.3f}")

    for r in results:
        r.significant = (r.csp_A is not None and final_stats is not None and r.csp_A >= cutoff) if values else None
        if final_stats and r.csp_A is not None:
            if final_stats.sd > 0.0:
                r.z_score = (r.csp_A - final_stats.mean) / final_stats.sd
            else:
                r.z_score = 0.0
        else:
            r.z_score = None

        if threshold_1sd is not None and r.csp_A is not None:
            r.significant_1sd = (r.csp_A >= threshold_1sd)
        else:
            r.significant_1sd = None

        if threshold_2sd is not None and r.csp_A is not None:
            r.significant_2sd = (r.csp_A >= threshold_2sd)
        else:
            r.significant_2sd = None

    return results


def compute_csp_multiple_saveframes_ca(
    apo_sequences: List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]],
    holo_sequences: List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]],
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    threshold_config: Optional[Dict] = None,
    *,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CA-inclusive CSPs for all possible apo-holo sequence pairs from multiple saveframes.

    Only aligns sequences that have H, N, and CA shifts.
    Returns CSPs for the alignment with the best alignment score.
    """
    from .align import align_global

    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))

    if _verbose:
        print(f"[CSP-CA] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")

    best_alignment_score = float('-inf')
    best_alignment_info = None

    for apo_seq, H_apo, N_apo, CA_apo, apo_saveframe in apo_sequences:
        for holo_seq, H_holo, N_holo, CA_holo, holo_saveframe in holo_sequences:
            if _verbose:
                print(f"[CSP-CA] Aligning apo({apo_saveframe}) vs holo({holo_saveframe})")

            aligned_apo, aligned_holo, mapping, alignment_score = align_global(apo_seq, holo_seq)

            if _verbose:
                print(f"[CSP-CA] Alignment score: {alignment_score}, mapped pairs: {len(mapping)}")

            if alignment_score > best_alignment_score:
                best_alignment_score = alignment_score
                best_alignment_info = (
                    aligned_apo,
                    aligned_holo,
                    H_apo,
                    N_apo,
                    CA_apo,
                    H_holo,
                    N_holo,
                    CA_holo,
                    apo_saveframe,
                    holo_saveframe,
                )

    if best_alignment_info is None:
        if _verbose:
            print("[CSP-CA] No valid alignment found for CA-inclusive CSPs")
        return []

    (
        aligned_apo,
        aligned_holo,
        H_apo,
        N_apo,
        CA_apo,
        H_holo,
        N_holo,
        CA_holo,
        apo_saveframe,
        holo_saveframe,
    ) = best_alignment_info

    if _verbose:
        print(
            f"[CSP-CA] Computing CA-inclusive CSPs for best alignment: "
            f"apo({apo_saveframe}) vs holo({holo_saveframe}) with score {best_alignment_score}"
        )

    results = compute_csp_from_aligned_sequences_ca(
        aligned_apo,
        aligned_holo,
        H_apo,
        N_apo,
        CA_apo,
        H_holo,
        N_holo,
        CA_holo,
        threshold_config,
        enable_referencing=True,
        referencing_method=referencing_method,
        grid_params=grid_params,
        target_id=(target_id or holo_pdb),
    )

    if _verbose:
        valid_csps = sum(1 for r in results if r.csp_A is not None)
        print(f"[CSP-CA] Valid CA-inclusive CSPs: {valid_csps}")

    return results


