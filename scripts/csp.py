"""
Compute chemical shift perturbations (CSP) using equation A and summarize.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import os
import sys
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


def _csp_ha_term(d_ha: float) -> float:
    w = getattr(compute, "csp_delta_ha_scale", 1.0)
    return (d_ha * w) ** 2


def _csp_ca_term(d_ca: float) -> float:
    w = compute.csp_delta_ca_scale
    return (d_ca * w) ** 2


def _max_abs_delta_n_ppm() -> float:
    return float(getattr(compute, "max_abs_delta_n_ppm", 15.0))


def _exceeds_max_abs_delta_n(
    n_apo: Optional[float],
    n_holo_raw: Optional[float],
    *,
    threshold_ppm: Optional[float] = None,
) -> bool:
    """True when |N_holo_raw − N_apo| exceeds the configured ppm cutoff."""
    if n_apo is None or n_holo_raw is None:
        return False
    limit = _max_abs_delta_n_ppm() if threshold_ppm is None else float(threshold_ppm)
    return abs(float(n_holo_raw) - float(n_apo)) > limit


def _residues_identity_matched(aa_apo: str, aa_holo: str) -> bool:
    """True when apo/holo one-letter codes match (same residue; not gap/unknown)."""
    a = (aa_apo or "").strip().upper()
    b = (aa_holo or "").strip().upper()
    return bool(a) and a == b and a not in {"-", "?", "X"}


def collect_large_dn_exclusions(
    results: List["CSPResult"],
    *,
    threshold_ppm: Optional[float] = None,
) -> List[Dict[str, object]]:
    """Build manifest rows for residues excluded due to large |ΔN_raw|."""
    limit = _max_abs_delta_n_ppm() if threshold_ppm is None else float(threshold_ppm)
    rows: List[Dict[str, object]] = []
    for r in results:
        if not getattr(r, "excluded_large_dn", False):
            continue
        n_apo = r.N_apo
        n_holo_raw = r.N_holo_original if r.N_holo_original is not None else r.N_holo
        if n_apo is None or n_holo_raw is None:
            continue
        dN_raw = float(n_holo_raw) - float(n_apo)
        rows.append(
            {
                "apo_resi": r.apo_index,
                "apo_aa": r.apo_aa,
                "holo_resi": r.holo_index,
                "holo_aa": r.holo_aa,
                "N_apo": float(n_apo),
                "N_holo_raw": float(n_holo_raw),
                "dN_raw": dN_raw,
                "abs_dN_raw": abs(dN_raw),
                "threshold_ppm": limit,
            }
        )
    return rows


LARGE_DN_EXCLUSION_FIELDS = [
    "apo_bmrb",
    "holo_bmrb",
    "holo_pdb",
    "target_dir",
    "apo_resi",
    "apo_aa",
    "holo_resi",
    "holo_aa",
    "N_apo",
    "N_holo_raw",
    "dN_raw",
    "abs_dN_raw",
    "threshold_ppm",
]


def write_large_dn_exclusions_csv(
    path: str,
    rows: List[Dict[str, object]],
    *,
    apo_bmrb: str = "",
    holo_bmrb: str = "",
    holo_pdb: str = "",
    target_dir: str = "",
    write_if_empty: bool = False,
) -> None:
    """Write exclusion manifest CSV (creates parent dirs).

    By default skips writing when ``rows`` is empty (per-target files).
    Pass ``write_if_empty=True`` for aggregate manifests that should be refreshed.
    """
    import csv as _csv

    if not rows and not write_if_empty:
        return
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=LARGE_DN_EXCLUSION_FIELDS)
        w.writeheader()
        for row in rows:
            out = {
                "apo_bmrb": apo_bmrb or row.get("apo_bmrb", ""),
                "holo_bmrb": holo_bmrb or row.get("holo_bmrb", ""),
                "holo_pdb": holo_pdb or row.get("holo_pdb", ""),
                "target_dir": target_dir or row.get("target_dir", ""),
                "apo_resi": row.get("apo_resi", ""),
                "apo_aa": row.get("apo_aa", ""),
                "holo_resi": row.get("holo_resi", ""),
                "holo_aa": row.get("holo_aa", ""),
                "N_apo": f"{float(row['N_apo']):.4f}" if row.get("N_apo") is not None else "",
                "N_holo_raw": f"{float(row['N_holo_raw']):.4f}" if row.get("N_holo_raw") is not None else "",
                "dN_raw": f"{float(row['dN_raw']):.4f}" if row.get("dN_raw") is not None else "",
                "abs_dN_raw": f"{float(row['abs_dN_raw']):.4f}" if row.get("abs_dN_raw") is not None else "",
                "threshold_ppm": f"{float(row.get('threshold_ppm', _max_abs_delta_n_ppm())):.1f}",
            }
            w.writerow(out)


def scan_csp_tables_for_large_dn_exclusions(
    outputs_dir: str,
    *,
    threshold_ppm: Optional[float] = None,
) -> List[Dict[str, object]]:
    """Scan outputs/*/csp_table.csv for |ΔN_raw| > threshold; return manifest rows."""
    import csv as _csv
    import glob as _glob

    limit = _max_abs_delta_n_ppm() if threshold_ppm is None else float(threshold_ppm)
    rows: List[Dict[str, object]] = []
    pattern = os.path.join(outputs_dir, "*", "csp_table.csv")
    for table_path in sorted(_glob.glob(pattern)):
        target_dir = os.path.basename(os.path.dirname(table_path))
        with open(table_path, newline="") as f:
            reader = _csv.DictReader(f)
            for row in reader:
                try:
                    n_apo_s = (row.get("N_apo") or "").strip()
                    n_holo_s = (row.get("N_holo_original") or row.get("N_holo") or "").strip()
                    if not n_apo_s or not n_holo_s:
                        continue
                    n_apo = float(n_apo_s)
                    n_holo_raw = float(n_holo_s)
                except (TypeError, ValueError):
                    continue
                dN_raw = n_holo_raw - n_apo
                if abs(dN_raw) <= limit:
                    continue
                rows.append(
                    {
                        "apo_bmrb": (row.get("apo_bmrb") or "").strip(),
                        "holo_bmrb": (row.get("holo_bmrb") or "").strip(),
                        "holo_pdb": (row.get("holo_pdb") or "").strip(),
                        "target_dir": target_dir,
                        "apo_resi": (row.get("apo_resi") or "").strip(),
                        "apo_aa": (row.get("apo_aa") or "").strip(),
                        "holo_resi": (row.get("holo_resi") or "").strip(),
                        "holo_aa": (row.get("holo_aa") or "").strip(),
                        "N_apo": n_apo,
                        "N_holo_raw": n_holo_raw,
                        "dN_raw": dN_raw,
                        "abs_dN_raw": abs(dN_raw),
                        "threshold_ppm": limit,
                    }
                )
    rows.sort(key=lambda r: (str(r.get("target_dir", "")), abs(float(r.get("dN_raw", 0))), str(r.get("apo_resi", ""))))
    return rows


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


def _floor_primary_hn_cutoff(cutoff: float) -> float:
    """Primary HN significance cutoff: max(cleaned-mean path, cutoff_05_ppm)."""
    return max(float(cutoff), float(thresholds.cutoff_05_ppm))


def percentile_linear(values: List[float], percentile: float) -> float:
    """
    Linear-interpolation percentile (numpy default style).

    Parameters
    ----------
    values : List[float]
        Sample values (need not be sorted).
    percentile : float
        Percentile in [0, 100].
    """
    if not values:
        raise ValueError("percentile_linear requires a non-empty values list")
    if percentile < 0.0 or percentile > 100.0:
        raise ValueError(f"percentile must be in [0, 100], got {percentile}")
    xs = sorted(values)
    if len(xs) == 1:
        return float(xs[0])
    rank = (percentile / 100.0) * (len(xs) - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return float(xs[lo])
    weight = rank - lo
    return float(xs[lo] * (1.0 - weight) + xs[hi] * weight)


def _apply_hn_significance_masks(
    results: List["CSPResult"],
    *,
    cutoff: Optional[float],
    final_stats: Optional[ThresholdComputation],
    csp_values: List[float],
    verbose: bool = False,
) -> None:
    """
    Assign legacy and alternate HN significance boolean masks on CSPResult rows.

    Sigma masks use cleaned mean/SD from iterative outlier removal (same stats as
    significant_1sd / significant_2sd). Raw-mean masks use mean/SD of the full HN
    CSP value list with no outlier removal. Percentile and fixed-ppm masks use the
    provided HN CSP value list (typically the values used for final thresholding).
    """
    threshold_1sd: Optional[float] = None
    threshold_2sd: Optional[float] = None
    threshold_sigma_0: Optional[float] = None
    threshold_max_05_cleaned: Optional[float] = None

    if final_stats and final_stats.cleaned_values:
        threshold_sigma_0 = final_stats.mean
        threshold_1sd = final_stats.mean + 1.0 * final_stats.sd
        threshold_2sd = final_stats.mean + 2.0 * final_stats.sd
        if verbose:
            print(
                f"[CSP] Additional thresholds: "
                f"sigma_0={threshold_sigma_0:.3f}, 1SD={threshold_1sd:.3f}, 2SD={threshold_2sd:.3f}"
            )

    top_10_frac = float(getattr(thresholds, "top_10_percentile_fraction", 0.10))
    top_5_frac = float(getattr(thresholds, "top_5_percentile_fraction", 0.05))
    ppm_03 = float(getattr(thresholds, "cutoff_03_ppm", 0.03))
    ppm_05 = float(getattr(thresholds, "cutoff_05_ppm", 0.05))
    ppm_10 = float(getattr(thresholds, "cutoff_10_ppm", 0.10))

    if final_stats and final_stats.cleaned_values:
        threshold_max_05_cleaned = max(ppm_05, float(final_stats.mean))
        if verbose:
            print(
                f"[CSP] max(0.05, cleaned mean) threshold: {threshold_max_05_cleaned:.3f}"
            )

    threshold_top_10: Optional[float] = None
    threshold_top_5: Optional[float] = None
    threshold_raw_mean: Optional[float] = None
    threshold_raw_1sd: Optional[float] = None
    threshold_raw_2sd: Optional[float] = None
    if csp_values:
        threshold_top_10 = percentile_linear(csp_values, (1.0 - top_10_frac) * 100.0)
        threshold_top_5 = percentile_linear(csp_values, (1.0 - top_5_frac) * 100.0)
        raw_mean, raw_sd = _compute_basic_stats(csp_values)
        threshold_raw_mean = raw_mean
        threshold_raw_1sd = raw_mean + 1.0 * raw_sd
        threshold_raw_2sd = raw_mean + 2.0 * raw_sd
        if verbose:
            print(
                f"[CSP] Percentile thresholds: "
                f"top_10>={threshold_top_10:.3f}, top_5>={threshold_top_5:.3f}"
            )
            print(
                f"[CSP] Raw (no outlier removal) thresholds: "
                f"mean={threshold_raw_mean:.3f}, 1SD={threshold_raw_1sd:.3f}, "
                f"2SD={threshold_raw_2sd:.3f}"
            )

    has_values = bool(csp_values)

    for r in results:
        if has_values and r.csp_A is not None and final_stats is not None and cutoff is not None:
            r.significant = r.csp_A >= cutoff
        else:
            r.significant = None

        if final_stats and r.csp_A is not None:
            if final_stats.sd > 0.0:
                r.z_score = (r.csp_A - final_stats.mean) / final_stats.sd
            else:
                r.z_score = 0.0
        else:
            r.z_score = None

        if threshold_1sd is not None and r.csp_A is not None:
            r.significant_1sd = r.csp_A >= threshold_1sd
        else:
            r.significant_1sd = None

        if threshold_2sd is not None and r.csp_A is not None:
            r.significant_2sd = r.csp_A >= threshold_2sd
        else:
            r.significant_2sd = None

        if threshold_sigma_0 is not None and r.csp_A is not None:
            r.significant_sigma_0 = r.csp_A >= threshold_sigma_0
            r.significant_sigma_1 = r.csp_A >= threshold_1sd if threshold_1sd is not None else None
            r.significant_sigma_2 = r.csp_A >= threshold_2sd if threshold_2sd is not None else None
        else:
            r.significant_sigma_0 = None
            r.significant_sigma_1 = None
            r.significant_sigma_2 = None

        if threshold_top_10 is not None and r.csp_A is not None:
            r.significant_top_10_percentile = r.csp_A >= threshold_top_10
        else:
            r.significant_top_10_percentile = None

        if threshold_top_5 is not None and r.csp_A is not None:
            r.significant_top_5_percentile = r.csp_A >= threshold_top_5
        else:
            r.significant_top_5_percentile = None

        if r.csp_A is not None:
            r.significant_03_ppm = r.csp_A >= ppm_03
            r.significant_05_ppm = r.csp_A >= ppm_05
            r.significant_10_ppm = r.csp_A >= ppm_10
        else:
            r.significant_03_ppm = None
            r.significant_05_ppm = None
            r.significant_10_ppm = None

        if threshold_raw_mean is not None and r.csp_A is not None:
            r.significant_raw_mean = r.csp_A >= threshold_raw_mean
            r.significant_raw_1sd = (
                r.csp_A >= threshold_raw_1sd if threshold_raw_1sd is not None else None
            )
            r.significant_raw_2sd = (
                r.csp_A >= threshold_raw_2sd if threshold_raw_2sd is not None else None
            )
        else:
            r.significant_raw_mean = None
            r.significant_raw_1sd = None
            r.significant_raw_2sd = None

        if threshold_max_05_cleaned is not None and r.csp_A is not None:
            r.significant_max_05_cleaned_mean = r.csp_A >= threshold_max_05_cleaned
        else:
            r.significant_max_05_cleaned_mean = None


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
    # Alternate HN significance masks (sigma / percentile / fixed ppm)
    significant_sigma_0: Optional[bool] = None
    significant_sigma_1: Optional[bool] = None
    significant_sigma_2: Optional[bool] = None
    significant_top_10_percentile: Optional[bool] = None
    significant_top_5_percentile: Optional[bool] = None
    significant_03_ppm: Optional[bool] = None
    significant_05_ppm: Optional[bool] = None
    significant_10_ppm: Optional[bool] = None
    # Raw mean/SD (no outlier removal) and floor-on-cleaned-mean masks
    significant_raw_mean: Optional[bool] = None
    significant_raw_1sd: Optional[bool] = None
    significant_raw_2sd: Optional[bool] = None
    significant_max_05_cleaned_mean: Optional[bool] = None
    # Referencing fields
    H_holo_original: Optional[float] = None
    N_holo_original: Optional[float] = None
    H_offset: Optional[float] = None
    N_offset: Optional[float] = None
    z_score: Optional[float] = None
    # Optional HA-related fields (used for HA-inclusive CSP analysis)
    HA_apo: Optional[float] = None
    HA_holo: Optional[float] = None
    HA_holo_original: Optional[float] = None
    HA_offset: Optional[float] = None
    dHA: Optional[float] = None
    # Optional CA-related fields (used for CA-inclusive CSP analysis)
    CA_apo: Optional[float] = None
    CA_holo: Optional[float] = None
    CA_holo_original: Optional[float] = None
    CA_offset: Optional[float] = None
    dCA: Optional[float] = None
    # True when |ΔN_raw| exceeded max_abs_delta_n_ppm and CSP was not recorded
    excluded_large_dn: bool = False
    # True when apo/holo amino acids differ (misaligned pair); CSP was not recorded
    excluded_aa_mismatch: bool = False


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


def run_offset_grid_search_ha_ca(
    points: List[Tuple[float, float, float, float]],
    *,
    ha_min: float,
    ha_max: float,
    ha_step: float,
    ca_min: float,
    ca_max: float,
    ca_step: float,
    cutoff: float,
) -> Dict[str, object]:
    """
    Grid search for HA/CA offsets that maximize the number of HA/CA CSPs < cutoff.

    Each point is a tuple:
      (HA_apo, CA_apo, HA_holo_orig, CA_holo_orig)
    using the HA/CA CSP formula:
      CSP = sqrt(1/2*((wHA*dHA)^2 + (wCA*dCA)^2))
    """
    ha_decimals = max(0, min(6, len(str(ha_step).split(".")[-1]) if "." in str(ha_step) else 0))
    ca_decimals = max(0, min(6, len(str(ca_step).split(".")[-1]) if "." in str(ca_step) else 0))

    ha_values = _frange(ha_min, ha_max, ha_step, ha_decimals)
    ca_values = _frange(ca_min, ca_max, ca_step, ca_decimals)

    count_matrix: List[List[int]] = []
    best_count = -1
    best_pairs: List[Tuple[float, float]] = []

    cutoff_sq = cutoff * cutoff
    for ha_off in ha_values:
        row: List[int] = []
        for ca_off in ca_values:
            count_lt = 0
            for ha_a, ca_a, ha_h_orig, ca_h_orig in points:
                dHA = (ha_h_orig + ha_off) - ha_a
                dCA = (ca_h_orig + ca_off) - ca_a
                inner = _csp_ha_term(dHA) + _csp_ca_term(dCA)
                if inner < 2.0 * cutoff_sq:
                    count_lt += 1
            row.append(count_lt)

            if count_lt > best_count:
                best_count = count_lt
                best_pairs = [(ca_off, ha_off)]
            elif count_lt == best_count:
                best_pairs.append((ca_off, ha_off))
        count_matrix.append(row)

    if best_pairs:
        best_ca_offset, best_ha_offset = min(
            best_pairs,
            key=lambda p: (abs(p[0]) + abs(p[1]), abs(p[0]), abs(p[1])),
        )
    else:
        best_ca_offset, best_ha_offset = 0.0, 0.0

    return {
        "ha_values": ha_values,
        "ca_values": ca_values,
        "count_matrix": count_matrix,
        "best_ha_offset": best_ha_offset,
        "best_ca_offset": best_ca_offset,
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
    if atom_upper not in ("H", "N", "CA", "HA"):
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
        elif atom_upper == "CA":
            apo = getattr(r, "CA_apo", None)
            holo_orig = getattr(r, "CA_holo_original", None)
        else:  # "HA"
            apo = getattr(r, "HA_apo", None)
            holo_orig = getattr(r, "HA_holo_original", None)
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
    elif atom_upper == "CA":
        min_offset = (grid_params or {}).get("ca_min", cfg.grid_ca_min)
        max_offset = (grid_params or {}).get("ca_max", cfg.grid_ca_max)
        step = (grid_params or {}).get("ca_step", cfg.grid_ca_step)
    else:  # "HA"
        min_offset = (grid_params or {}).get("ha_min", cfg.grid_ha_min)
        max_offset = (grid_params or {}).get("ha_max", cfg.grid_ha_max)
        step = (grid_params or {}).get("ha_step", cfg.grid_ha_step)

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
        elif atom_upper == "CA":
            apo = getattr(r, "CA_apo", None)
            holo_orig = getattr(r, "CA_holo_original", None)
        else:
            apo = getattr(r, "HA_apo", None)
            holo_orig = getattr(r, "HA_holo_original", None)

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


def _build_param_slug_ha_ca(*, ha_min: float, ha_max: float, ha_step: float, ca_min: float, ca_max: float, ca_step: float, cutoff: float) -> str:
    def fmt(v: float, decimals: int) -> str:
        return (f"{v:.{decimals}f}").rstrip('0').rstrip('.') if '.' in f"{v:.{decimals}f}" else f"{v:.{decimals}f}"
    return (
        f"HA_{fmt(ha_min, 3)}_{fmt(ha_max, 3)}_{fmt(ha_step, 3)}__"
        f"CA_{fmt(ca_min, 3)}_{fmt(ca_max, 3)}_{fmt(ca_step, 3)}__"
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
    # count_matrix is (H rows × N cols); transpose so ¹H is x, ¹⁵N is y
    data = np.array(matrix, dtype=float).T
    best_h = float(grid_result["best_h_offset"])  # type: ignore
    best_n = float(grid_result["best_n_offset"])  # type: ignore

    # Small figure: ~2 inch square heatmap, extra space for horizontal colorbar
    fig, ax = plt.subplots(figsize=(4, 3.5))
    im = ax.imshow(data, origin='lower', aspect='auto', cmap='viridis')
    # Horizontal colorbar below heatmap, pad adds buffer between x-axis label and colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.8, aspect=25)
    cbar.set_label(r'Number of aligned peaks ($\Delta\delta$ < 0.05 ppm)', fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlabel('¹H Offset (ppm)', fontsize=10)
    ax.set_ylabel('¹⁵N Offset (ppm)', fontsize=10)
    ax.set_title(title, fontsize=11)
    # Mark optimal offset with black star (x=H, y=N after transpose)
    best_h_idx = int(np.argmin(np.abs(h_values - best_h)))
    best_n_idx = int(np.argmin(np.abs(n_values - best_n)))
    ax.scatter(best_h_idx, best_n_idx, marker='*', s=150, c='black', zorder=5, edgecolors='white', linewidths=0.5)
    # ~5 tick marks including 0 (linspace over symmetric range yields 0 at center)
    num_ticks = 5
    h_idx = np.linspace(0, len(h_values) - 1, num=num_ticks).astype(int)
    n_idx = np.linspace(0, len(n_values) - 1, num=num_ticks).astype(int)
    ax.set_xticks(h_idx)
    _tick_fs = 6
    ax.set_xticklabels(
        [f"{h_values[i]:.2f}" for i in h_idx],
        fontsize=_tick_fs,
        rotation=45,
        ha='right',
    )
    ax.set_yticks(n_idx)
    ax.set_yticklabels([f"{n_values[i]:.2f}" for i in n_idx], fontsize=_tick_fs)
    # Add optimal offset result as text (H, N to match axis order)
    ax.text(
        0.02,
        0.98,
        f"(¹H offset, ¹⁵N offset) = ({best_h:.2f}, {best_n:.2f}) ppm",
        transform=ax.transAxes,
        fontsize=9,
        va='top',
        ha='left',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
    )
    plt.tight_layout(pad=1.5)
    _ensure_dir(os.path.dirname(output_png))
    fig.savefig(output_png, dpi=200, bbox_inches='tight')
    plt.close(fig)


def save_grid_heatmap_png_ha_ca(output_png: str, grid_result: Dict[str, object], title: str) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore
    except Exception:
        return
    ha_values = np.array(grid_result["ha_values"])  # type: ignore
    ca_values = np.array(grid_result["ca_values"])  # type: ignore
    matrix = grid_result["count_matrix"]  # type: ignore
    data = np.array(matrix, dtype=float)
    best_ha = float(grid_result["best_ha_offset"])  # type: ignore
    best_ca = float(grid_result["best_ca_offset"])  # type: ignore

    fig, ax = plt.subplots(figsize=(4, 3.5))
    im = ax.imshow(data, origin='lower', aspect='auto', cmap='viridis')
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.25, shrink=0.8, aspect=25)
    cbar.set_label(r'Number of aligned peaks ($\Delta\delta$ < 0.05 ppm)', fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlabel('Cα Offset (ppm)', fontsize=10)
    ax.set_ylabel('Hα Offset (ppm)', fontsize=10)
    ax.set_title(title, fontsize=11)

    best_ha_idx = int(np.argmin(np.abs(ha_values - best_ha)))
    best_ca_idx = int(np.argmin(np.abs(ca_values - best_ca)))
    ax.scatter(best_ca_idx, best_ha_idx, marker='*', s=150, c='black', zorder=5, edgecolors='white', linewidths=0.5)

    num_ticks = 5
    ca_idx = np.linspace(0, len(ca_values) - 1, num=num_ticks).astype(int)
    ha_idx = np.linspace(0, len(ha_values) - 1, num=num_ticks).astype(int)
    ax.set_xticks(ca_idx)
    _tick_fs = 6
    ax.set_xticklabels(
        [f"{ca_values[i]:.2f}" for i in ca_idx],
        fontsize=_tick_fs,
        rotation=45,
        ha='right',
    )
    ax.set_yticks(ha_idx)
    ax.set_yticklabels([f"{ha_values[i]:.2f}" for i in ha_idx], fontsize=_tick_fs)
    ax.text(
        0.02,
        0.98,
        f"(Cα offset, Hα offset) = ({best_ca:.2f}, {best_ha:.2f}) ppm",
        transform=ax.transAxes,
        fontsize=9,
        va='top',
        ha='left',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
    )
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


def save_grid_csv_ha_ca(output_csv: str, grid_result: Dict[str, object], *,
                        ha_min: float, ha_max: float, ha_step: float,
                        ca_min: float, ca_max: float, ca_step: float,
                        cutoff: float) -> None:
    _ensure_dir(os.path.dirname(output_csv))
    ha_values: List[float] = grid_result["ha_values"]  # type: ignore
    ca_values: List[float] = grid_result["ca_values"]  # type: ignore
    matrix: List[List[int]] = grid_result["count_matrix"]  # type: ignore
    best_ha: float = grid_result["best_ha_offset"]  # type: ignore
    best_ca: float = grid_result["best_ca_offset"]  # type: ignore
    best_count: int = grid_result["best_count"]  # type: ignore

    with open(output_csv, 'w') as f:
        f.write(f"# ha_min,{ha_min}\n")
        f.write(f"# ha_max,{ha_max}\n")
        f.write(f"# ha_step,{ha_step}\n")
        f.write(f"# ca_min,{ca_min}\n")
        f.write(f"# ca_max,{ca_max}\n")
        f.write(f"# ca_step,{ca_step}\n")
        f.write(f"# cutoff,{cutoff}\n")
        f.write("ha_values," + ",".join(str(v) for v in ha_values) + "\n")
        f.write("ca_values," + ",".join(str(v) for v in ca_values) + "\n")
        for i, ha in enumerate(ha_values):
            row = matrix[i]
            f.write(str(ha) + "," + ",".join(str(c) for c in row) + "\n")
        f.write(f"best_ca_offset,{best_ca}\n")
        f.write(f"best_ha_offset,{best_ha}\n")
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


def load_grid_csv_ha_ca(input_csv: str) -> Optional[Dict[str, object]]:
    if not os.path.exists(input_csv):
        return None
    try:
        ha_values: List[float] = []
        ca_values: List[float] = []
        matrix: List[List[int]] = []
        best_ha = 0.0
        best_ca = 0.0
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
                if parts[0] == 'ha_values':
                    ha_values = [float(x) for x in parts[1:] if x]
                    continue
                if parts[0] == 'ca_values':
                    ca_values = [float(x) for x in parts[1:] if x]
                    stage = 'matrix'
                    continue
                if parts[0] == 'best_ca_offset':
                    best_ca = float(parts[1])
                    continue
                if parts[0] == 'best_ha_offset':
                    best_ha = float(parts[1])
                    continue
                if parts[0] == 'best_count':
                    best_count = int(float(parts[1]))
                    continue
                if stage == 'matrix':
                    row_counts = [int(float(x)) for x in parts[1:] if x]
                    matrix.append(row_counts)
        if not ha_values or not ca_values or not matrix:
            return None
        return {
            "ha_values": ha_values,
            "ca_values": ca_values,
            "count_matrix": matrix,
            "best_ha_offset": best_ha,
            "best_ca_offset": best_ca,
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
    output_root: Optional[str] = None,
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
        excluded_aa_mismatch = not _residues_identity_matched(aa_apo, aa_holo)
        excluded_large_dn = _exceeds_max_abs_delta_n(n_a, n_h)
        csp_val: Optional[float] = None
        if (
            dH is not None
            and dN is not None
            and not excluded_large_dn
            and not excluded_aa_mismatch
        ):
            csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
            values.append(csp_val)
        if _verbose and excluded_aa_mismatch:
            print(
                f"[CSP] apo {i_apo}{aa_apo} ↔ holo {i_holo}{aa_holo} : "
                f"excluded AA mismatch (require identical residues)"
            )
        elif _verbose and excluded_large_dn:
            _dN = f"{dN:.3f}" if dN is not None else "NA"
            print(
                f"[CSP] apo {i_apo}{aa_apo} ↔ holo {i_holo}{aa_holo} : "
                f"excluded |ΔN|={_dN} > {_max_abs_delta_n_ppm():.1f} ppm"
            )
        elif _verbose and (dH is not None or dN is not None):
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
                excluded_large_dn=excluded_large_dn,
                excluded_aa_mismatch=excluded_aa_mismatch,
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
            cutoff = _floor_primary_hn_cutoff(threshold_info.threshold)
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
                output_base_dir = output_root or paths.outputs_dir
                out_dir = os.path.join(output_base_dir, target_id) if target_id else None
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

                            if getattr(r, "excluded_large_dn", False) or getattr(
                                r, "excluded_aa_mismatch", False
                            ):
                                r.H_holo = h_h_ref
                                r.N_holo = n_h_ref
                                r.dH = dH_ref
                                r.dN = dN_ref
                                r.csp_A = None
                                r.H_offset = h_offset
                                r.N_offset = n_offset
                                continue
                            
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
                        if absolute_cutoff is not None:
                            cutoff = threshold_info_ref.threshold
                        else:
                            cutoff = _floor_primary_hn_cutoff(threshold_info_ref.threshold)
                        final_stats = threshold_info_ref
                        if _verbose:
                            print(f"[CSP] Referenced CSPs: {len(values_referenced)} mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}")
                            print(f"[CSP] Referenced threshold: {cutoff:.3f}")

    mask_values = values_referenced if values_referenced else values
    _apply_hn_significance_masks(
        results,
        cutoff=cutoff if values else None,
        final_stats=final_stats if values else None,
        csp_values=mask_values,
        verbose=_verbose,
    )

    return results


def _unpack_chem_shift_entry(entry):
    """Support 6-tuples (legacy) and 7-tuples (seq_id_min, saveframe)."""
    if len(entry) == 7:
        seq, H, N, CA, HA, seq_id_min, name = entry
        return seq, H, N, CA, HA, int(seq_id_min), name
    seq, H, N, CA, HA, name = entry
    keys = set(H) | set(N) | set(CA) | set(HA)
    seq_id_min = min(keys) if keys else 1
    return seq, H, N, CA, HA, int(seq_id_min), name


def _best_seqid_alignment(apo_sequences, holo_sequences):
    """Pick apo/holo polymers with the most exact Seq_ID-offset AA matches."""
    from .align import align_by_seqid_offset

    best_score = float("-inf")
    best = None
    for apo_entry in apo_sequences:
        apo_seq, H_apo, N_apo, CA_apo, HA_apo, apo_sid_min, apo_sf = _unpack_chem_shift_entry(
            apo_entry
        )
        for holo_entry in holo_sequences:
            holo_seq, H_holo, N_holo, CA_holo, HA_holo, holo_sid_min, holo_sf = (
                _unpack_chem_shift_entry(holo_entry)
            )
            aligned_apo, aligned_holo, mapping, score, offset = align_by_seqid_offset(
                apo_seq,
                holo_seq,
                apo_seq_id_min=apo_sid_min,
                holo_seq_id_min=holo_sid_min,
                H_apo=H_apo,
                N_apo=N_apo,
                H_holo=H_holo,
                N_holo=N_holo,
            )
            if score > best_score:
                best_score = score
                best = {
                    "aligned_apo": aligned_apo,
                    "aligned_holo": aligned_holo,
                    "mapping": mapping,
                    "score": score,
                    "offset": offset,
                    "H_apo": H_apo,
                    "N_apo": N_apo,
                    "CA_apo": CA_apo,
                    "HA_apo": HA_apo,
                    "H_holo": H_holo,
                    "N_holo": N_holo,
                    "CA_holo": CA_holo,
                    "HA_holo": HA_holo,
                    "apo_saveframe": apo_sf,
                    "holo_saveframe": holo_sf,
                    "apo_seq": apo_seq,
                    "holo_seq": holo_seq,
                    "apo_seq_id_min": apo_sid_min,
                    "holo_seq_id_min": holo_sid_min,
                }
    return best


def _csp_from_seqid_mapping(
    mapping: List[Tuple[int, int]],
    H_apo: Dict[int, float],
    N_apo: Dict[int, float],
    H_holo: Dict[int, float],
    N_holo: Dict[int, float],
    apo_aa_by_sid: Dict[int, str],
    holo_aa_by_sid: Dict[int, str],
    threshold_config: Optional[Dict],
    *,
    referencing_method: Optional[str],
    grid_params: Optional[Dict],
    target_id: Optional[str],
    output_root: Optional[str],
) -> List[CSPResult]:
    """Compute HN CSPs for exact Seq_ID pairs; remap indices back to Seq_IDs."""
    aligned_apo = "".join(apo_aa_by_sid[a] for a, _h in mapping)
    aligned_holo = "".join(holo_aa_by_sid[h] for _a, h in mapping)
    H_apo_d = {i: H_apo[a] for i, (a, _h) in enumerate(mapping, 1) if a in H_apo}
    N_apo_d = {i: N_apo[a] for i, (a, _h) in enumerate(mapping, 1) if a in N_apo}
    H_holo_d = {i: H_holo[h] for i, (_a, h) in enumerate(mapping, 1) if h in H_holo}
    N_holo_d = {i: N_holo[h] for i, (_a, h) in enumerate(mapping, 1) if h in N_holo}
    results = compute_csp_from_aligned_sequences(
        aligned_apo,
        aligned_holo,
        H_apo_d,
        N_apo_d,
        H_holo_d,
        N_holo_d,
        threshold_config,
        enable_referencing=True,
        referencing_method=referencing_method,
        grid_params=grid_params,
        target_id=target_id,
        output_root=output_root,
    )
    for r in results:
        asid, hsid = mapping[r.apo_index - 1]
        r.apo_index = asid
        r.holo_index = hsid
    return results


def compute_csp_multiple_saveframes(
    apo_sequences: List[Tuple],
    holo_sequences: List[Tuple],
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    threshold_config: Optional[Dict] = None,
    *,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
    output_root: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CSPs for all possible apo-holo sequence pairs from multiple saveframes.
    
    Pairs residues by BMRB Seq_ID offset (exact AA/Comp_ID matches only).
    Returns CSPs for the polymer pair with the most exact matches.
    """
    # Debug verbosity check
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))
    
    if _verbose:
        print(f"[CSP] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")
    
    best = _best_seqid_alignment(apo_sequences, holo_sequences)
    if best is None or not best["mapping"]:
        return []

    if _verbose:
        print(
            f"[CSP] Seq_ID-offset alignment: apo({best['apo_saveframe']}) vs "
            f"holo({best['holo_saveframe']}) offset={best['offset']} "
            f"exact_matches={len(best['mapping'])} score={best['score']}"
        )

    apo_aa = {
        best["apo_seq_id_min"] + i: aa for i, aa in enumerate(best["apo_seq"])
    }
    holo_aa = {
        best["holo_seq_id_min"] + i: aa for i, aa in enumerate(best["holo_seq"])
    }
    results = _csp_from_seqid_mapping(
        best["mapping"],
        best["H_apo"],
        best["N_apo"],
        best["H_holo"],
        best["N_holo"],
        apo_aa,
        holo_aa,
        threshold_config,
        referencing_method=referencing_method,
        grid_params=grid_params,
        target_id=(target_id or holo_pdb),
        output_root=output_root,
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
    output_root: Optional[str] = None,
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
            excluded_aa_mismatch = not _residues_identity_matched(aa_apo, aa_holo)
            excluded_large_dn = _exceeds_max_abs_delta_n(n_a, n_h)
            
            # Calculate CSP
            csp_val: Optional[float] = None
            if (
                dH is not None
                and dN is not None
                and not excluded_large_dn
                and not excluded_aa_mismatch
            ):
                csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
                values.append(csp_val)
            
            if _verbose and excluded_aa_mismatch:
                print(
                    f"[CSP] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : "
                    f"excluded AA mismatch (require identical residues)"
                )
            elif _verbose and excluded_large_dn:
                _dN = f"{dN:.3f}" if dN is not None else "NA"
                print(
                    f"[CSP] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : "
                    f"excluded |ΔN|={_dN} > {_max_abs_delta_n_ppm():.1f} ppm"
                )
            elif _verbose and (dH is not None or dN is not None):
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
                    excluded_large_dn=excluded_large_dn,
                    excluded_aa_mismatch=excluded_aa_mismatch,
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
            cutoff = _floor_primary_hn_cutoff(threshold_info.threshold)
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
                output_base_dir = output_root or paths.outputs_dir
                out_dir = os.path.join(output_base_dir, target_id) if target_id else None
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

                            if getattr(r, "excluded_large_dn", False) or getattr(
                                r, "excluded_aa_mismatch", False
                            ):
                                r.H_holo = h_h_ref
                                r.N_holo = n_h_ref
                                r.dH = dH_ref
                                r.dN = dN_ref
                                r.csp_A = None
                                r.H_offset = h_offset
                                r.N_offset = n_offset
                                continue
                            
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
                        if absolute_cutoff is not None:
                            cutoff = threshold_info_ref.threshold
                        else:
                            cutoff = _floor_primary_hn_cutoff(threshold_info_ref.threshold)
                        final_stats = threshold_info_ref
                        if _verbose:
                            print(f"[CSP] Referenced CSPs: {len(values_referenced)} mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}")
                            print(f"[CSP] Referenced threshold: {cutoff:.3f}")

    mask_values = values_referenced if values_referenced else values
    _apply_hn_significance_masks(
        results,
        cutoff=cutoff if values else None,
        final_stats=final_stats if values else None,
        csp_values=mask_values,
        verbose=_verbose,
    )
    
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
    output_root: Optional[str] = None,
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
            excluded_aa_mismatch = not _residues_identity_matched(aa_apo, aa_holo)
            excluded_large_dn = _exceeds_max_abs_delta_n(n_a, n_h)

            csp_val: Optional[float] = None
            if (
                dH is not None
                and dN is not None
                and dCA is not None
                and not excluded_large_dn
                and not excluded_aa_mismatch
            ):
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

            if _verbose and excluded_aa_mismatch:
                print(
                    f"[CSP-CA] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : "
                    f"excluded AA mismatch (require identical residues)"
                )
            elif _verbose and excluded_large_dn:
                _dN = f"{dN:.3f}" if dN is not None else "NA"
                print(
                    f"[CSP-CA] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : "
                    f"excluded |ΔN|={_dN} > {_max_abs_delta_n_ppm():.1f} ppm"
                )
            elif _verbose and (dH is not None or dN is not None or dCA is not None):
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
                    excluded_large_dn=excluded_large_dn,
                    excluded_aa_mismatch=excluded_aa_mismatch,
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

                            if getattr(r, "excluded_large_dn", False) or getattr(
                                r, "excluded_aa_mismatch", False
                            ):
                                r.H_holo = h_h_ref
                                r.N_holo = n_h_ref
                                r.dH = dH_ref
                                r.dN = dN_ref
                                r.CA_holo = ca_h_ref
                                r.dCA = dCA_ref
                                r.csp_A = None
                                r.H_offset = h_offset
                                r.N_offset = n_offset
                                r.CA_offset = ca_offset
                                continue

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
    apo_sequences: List[Tuple],
    holo_sequences: List[Tuple],
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    threshold_config: Optional[Dict] = None,
    *,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
    output_root: Optional[str] = None,
) -> List[CSPResult]:
    """Compute CA-inclusive CSPs using Seq_ID-offset exact AA matches."""
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))

    if _verbose:
        print(f"[CSP-CA] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")

    best = _best_seqid_alignment(apo_sequences, holo_sequences)
    if best is None or not best["mapping"]:
        if _verbose:
            print("[CSP-CA] No valid Seq_ID alignment found for CA-inclusive CSPs")
        return []

    mapping = best["mapping"]
    apo_aa = {best["apo_seq_id_min"] + i: aa for i, aa in enumerate(best["apo_seq"])}
    holo_aa = {best["holo_seq_id_min"] + i: aa for i, aa in enumerate(best["holo_seq"])}
    aligned_apo = "".join(apo_aa[a] for a, _h in mapping)
    aligned_holo = "".join(holo_aa[h] for _a, h in mapping)
    H_apo = {i: best["H_apo"][a] for i, (a, _h) in enumerate(mapping, 1) if a in best["H_apo"]}
    N_apo = {i: best["N_apo"][a] for i, (a, _h) in enumerate(mapping, 1) if a in best["N_apo"]}
    CA_apo = {i: best["CA_apo"][a] for i, (a, _h) in enumerate(mapping, 1) if a in best["CA_apo"]}
    H_holo = {i: best["H_holo"][h] for i, (_a, h) in enumerate(mapping, 1) if h in best["H_holo"]}
    N_holo = {i: best["N_holo"][h] for i, (_a, h) in enumerate(mapping, 1) if h in best["N_holo"]}
    CA_holo = {i: best["CA_holo"][h] for i, (_a, h) in enumerate(mapping, 1) if h in best["CA_holo"]}

    if _verbose:
        print(
            f"[CSP-CA] Seq_ID-offset: apo({best['apo_saveframe']}) vs "
            f"holo({best['holo_saveframe']}) offset={best['offset']} "
            f"exact_matches={len(mapping)}"
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
        output_root=output_root,
    )
    for r in results:
        asid, hsid = mapping[r.apo_index - 1]
        r.apo_index = asid
        r.holo_index = hsid

    if _verbose:
        valid_csps = sum(1 for r in results if r.csp_A is not None)
        print(f"[CSP-CA] Valid CA-inclusive CSPs: {valid_csps}")

    return results


def compute_csp_from_aligned_sequences_ha_ca(
    aligned_apo: str,
    aligned_holo: str,
    HA_apo: Dict[int, float],
    CA_apo: Dict[int, float],
    HA_holo: Dict[int, float],
    CA_holo: Dict[int, float],
    threshold_config: Optional[Dict] = None,
    enable_referencing: bool = True,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
    output_root: Optional[str] = None,
) -> List[CSPResult]:
    """Compute HA/CA CSPs from aligned sequences, handling gaps."""
    results: List[CSPResult] = []
    values: List[float] = []

    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))

    HA_holo_original = HA_holo.copy()
    CA_holo_original = CA_holo.copy()

    apo_pos = 0
    holo_pos = 0
    grid_points: List[Tuple[float, float, float, float]] = []

    for aa_apo, aa_holo in zip(aligned_apo, aligned_holo):
        if aa_apo != "-":
            apo_pos += 1
        if aa_holo != "-":
            holo_pos += 1

        if aa_apo != "-" and aa_holo != "-":
            ha_a = HA_apo.get(apo_pos)
            ca_a = CA_apo.get(apo_pos)
            ha_h = HA_holo.get(holo_pos)
            ca_h = CA_holo.get(holo_pos)

            dHA = (ha_h - ha_a) if (ha_h is not None and ha_a is not None) else None
            dCA = (ca_h - ca_a) if (ca_h is not None and ca_a is not None) else None
            excluded_aa_mismatch = not _residues_identity_matched(aa_apo, aa_holo)

            csp_val: Optional[float] = None
            if dHA is not None and dCA is not None and not excluded_aa_mismatch:
                csp_val = math.sqrt(0.5 * (_csp_ha_term(dHA) + _csp_ca_term(dCA)))
                values.append(csp_val)
                grid_points.append(
                    (
                        float(ha_a),
                        float(ca_a),
                        float(HA_holo_original.get(holo_pos, ha_h)),
                        float(CA_holo_original.get(holo_pos, ca_h)),
                    )
                )

            if _verbose and excluded_aa_mismatch:
                print(
                    f"[CSP-HA-CA] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : "
                    f"excluded AA mismatch (require identical residues)"
                )
            elif _verbose and (dHA is not None or dCA is not None):
                _dHA = f"{dHA:.3f}" if dHA is not None else "NA"
                _dCA = f"{dCA:.3f}" if dCA is not None else "NA"
                _C = f"{csp_val:.3f}" if csp_val is not None else "NA"
                print(f"[CSP-HA-CA] apo {apo_pos}{aa_apo} ↔ holo {holo_pos}{aa_holo} : dHA={_dHA} dCA={_dCA} C={_C}")

            results.append(
                CSPResult(
                    apo_index=apo_pos,
                    holo_index=holo_pos,
                    apo_aa=aa_apo,
                    holo_aa=aa_holo,
                    H_apo=None,
                    N_apo=None,
                    H_holo=None,
                    N_holo=None,
                    dH=None,
                    dN=None,
                    csp_A=csp_val,
                    significant=None,
                    significant_1sd=None,
                    significant_2sd=None,
                    HA_apo=ha_a,
                    HA_holo=ha_h,
                    HA_holo_original=HA_holo_original.get(holo_pos),
                    HA_offset=None,
                    dHA=dHA,
                    CA_apo=ca_a,
                    CA_holo=ca_h,
                    CA_holo_original=CA_holo_original.get(holo_pos),
                    CA_offset=None,
                    dCA=dCA,
                    excluded_aa_mismatch=excluded_aa_mismatch,
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
                print(f"[CSP-HA-CA] Using absolute cutoff: {cutoff:.3f}")
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
                print(f"[CSP-HA-CA] values={len(values)} mean={threshold_info.mean:.3f} sd={threshold_info.sd:.3f}")
                print(f"[CSP-HA-CA] outlier removal: {threshold_info.iterations} iterations, {threshold_info.outliers_removed} outliers removed")
                print(f"[CSP-HA-CA] final threshold: {cutoff:.3f}")

        final_stats = threshold_info

        if enable_referencing and grid_points:
            if _verbose:
                print("[CSP-HA-CA] Computing referencing offsets (2D grid)")

            method = (referencing_method or getattr(_Referencing, 'method', 'grid'))
            ha_offset = 0.0
            ca_offset = 0.0

            if method == 'grid':
                cfg = _Referencing()
                ha_min = (grid_params or {}).get('ha_min', cfg.grid_ha_min)
                ha_max = (grid_params or {}).get('ha_max', cfg.grid_ha_max)
                ha_step = (grid_params or {}).get('ha_step', cfg.grid_ha_step)
                ca_min = (grid_params or {}).get('ca_min', cfg.grid_ca_min)
                ca_max = (grid_params or {}).get('ca_max', cfg.grid_ca_max)
                ca_step = (grid_params or {}).get('ca_step', cfg.grid_ca_step)
                gs_cutoff = (grid_params or {}).get('cutoff', cfg.grid_cutoff)

                if _verbose:
                    print(
                        f"[CSP-HA-CA] Grid search over HA[{ha_min},{ha_max},{ha_step}] "
                        f"and CA[{ca_min},{ca_max},{ca_step}] with cutoff={gs_cutoff}"
                    )

                slug = _build_param_slug_ha_ca(
                    ha_min=ha_min,
                    ha_max=ha_max,
                    ha_step=ha_step,
                    ca_min=ca_min,
                    ca_max=ca_max,
                    ca_step=ca_step,
                    cutoff=float(gs_cutoff),
                )
                output_base_dir = output_root or paths.outputs_dir
                out_dir = os.path.join(output_base_dir, target_id) if target_id else None
                csv_path = os.path.join(out_dir, f"offset_grid_{slug}.csv") if out_dir else None
                png_path = os.path.join(out_dir, f"offset_grid_{slug}.png") if out_dir else None
                use_cache = bool(getattr(cfg, 'cache_results', True) and target_id and out_dir)

                loaded = False
                if use_cache and csv_path and os.path.exists(csv_path):
                    loaded_result = load_grid_csv_ha_ca(csv_path)
                    if loaded_result is not None:
                        grid_result = loaded_result
                        loaded = True
                        if _verbose:
                            print(f"[CSP-HA-CA] Loaded grid result from cache: {csv_path}")
                if not loaded:
                    grid_result = run_offset_grid_search_ha_ca(
                        grid_points,
                        ha_min=ha_min,
                        ha_max=ha_max,
                        ha_step=ha_step,
                        ca_min=ca_min,
                        ca_max=ca_max,
                        ca_step=ca_step,
                        cutoff=float(gs_cutoff),
                    )
                    if use_cache and csv_path:
                        save_grid_csv_ha_ca(
                            csv_path,
                            grid_result,
                            ha_min=ha_min, ha_max=ha_max, ha_step=ha_step,
                            ca_min=ca_min, ca_max=ca_max, ca_step=ca_step,
                            cutoff=float(gs_cutoff),
                        )
                        if _verbose:
                            print(f"[CSP-HA-CA] Saved grid result to: {csv_path}")

                if bool(getattr(cfg, 'save_heatmap', True)) and png_path:
                    title = f"Grid Search for Optimal Hα/Cα Offsets\n(PDB {target_id.upper()})" if target_id else "Grid Search for Optimal Hα/Cα Offsets"
                    save_grid_heatmap_png_ha_ca(png_path, grid_result, title)

                ha_offset = float(grid_result["best_ha_offset"])  # type: ignore
                ca_offset = float(grid_result["best_ca_offset"])  # type: ignore
                if _verbose:
                    print(
                        f"[CSP-HA-CA] Grid best offsets: "
                        f"CA_offset={ca_offset:.3f}, HA_offset={ha_offset:.3f} "
                        f"(count={grid_result['best_count']})"
                    )
            else:
                if _verbose:
                    print("[CSP-HA-CA] Mean-based referencing not implemented for HA/CA; skipping referencing")

            if ha_offset != 0.0 or ca_offset != 0.0:
                if _verbose:
                    print(
                        f"[CSP-HA-CA] Applying offsets: "
                        f"CA_offset={ca_offset:.3f}, HA_offset={ha_offset:.3f}"
                    )

                HA_holo_referenced = {pos: shift + ha_offset for pos, shift in HA_holo_original.items()}
                CA_holo_referenced = {pos: shift + ca_offset for pos, shift in CA_holo_original.items()}

                values_referenced = []
                for r in results:
                    if r.HA_holo_original is None or r.CA_holo_original is None:
                        continue

                    ha_h_ref = HA_holo_referenced.get(r.holo_index)
                    ca_h_ref = CA_holo_referenced.get(r.holo_index)
                    if ha_h_ref is None or ca_h_ref is None or r.HA_apo is None or r.CA_apo is None:
                        continue

                    dHA_ref = ha_h_ref - r.HA_apo
                    dCA_ref = ca_h_ref - r.CA_apo

                    if getattr(r, "excluded_aa_mismatch", False):
                        r.HA_holo = ha_h_ref
                        r.CA_holo = ca_h_ref
                        r.dHA = dHA_ref
                        r.dCA = dCA_ref
                        r.csp_A = None
                        r.HA_offset = ha_offset
                        r.CA_offset = ca_offset
                        continue

                    csp_val_ref = math.sqrt(0.5 * (_csp_ha_term(dHA_ref) + _csp_ca_term(dCA_ref)))
                    values_referenced.append(csp_val_ref)

                    r.HA_holo = ha_h_ref
                    r.CA_holo = ca_h_ref
                    r.dHA = dHA_ref
                    r.dCA = dCA_ref
                    r.csp_A = csp_val_ref
                    r.HA_offset = ha_offset
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
                                f"[CSP-HA-CA] Referenced CSPs: {len(values_referenced)} "
                                f"mean={threshold_info_ref.mean:.3f} sd={threshold_info_ref.sd:.3f}"
                            )
                            print(f"[CSP-HA-CA] Referenced threshold: {cutoff:.3f}")

    threshold_1sd: Optional[float] = None
    threshold_2sd: Optional[float] = None

    if final_stats and final_stats.cleaned_values:
        threshold_1sd = final_stats.mean + 1.0 * final_stats.sd
        threshold_2sd = final_stats.mean + 2.0 * final_stats.sd
        if _verbose:
            print(f"[CSP-HA-CA] Additional thresholds: 1SD={threshold_1sd:.3f}, 2SD={threshold_2sd:.3f}")

    for r in results:
        r.significant = (r.csp_A is not None and final_stats is not None and r.csp_A >= cutoff) if values else None
        if final_stats and r.csp_A is not None:
            r.z_score = (r.csp_A - final_stats.mean) / final_stats.sd if final_stats.sd > 0.0 else 0.0
        else:
            r.z_score = None

        r.significant_1sd = (r.csp_A >= threshold_1sd) if (threshold_1sd is not None and r.csp_A is not None) else None
        r.significant_2sd = (r.csp_A >= threshold_2sd) if (threshold_2sd is not None and r.csp_A is not None) else None

    return results


def _ha_ca_shift_pair_coverage_for_alignment(
    aligned_apo: str,
    aligned_holo: str,
    HA_apo: Dict[int, float],
    CA_apo: Dict[int, float],
    HA_holo: Dict[int, float],
    CA_holo: Dict[int, float],
) -> Tuple[float, int, int]:
    """
    Fraction of matched (non-gap) alignment columns where apo and holo each have
    both HA and CA chemical shifts — the positions that can produce a full HA-CA CSP.

    Returns (coverage, n_columns_with_full_ha_ca, n_matched_alignment_columns).
    Mirrors the position walk in ``compute_csp_from_aligned_sequences_ha_ca``.
    """
    apo_pos = 0
    holo_pos = 0
    n_matched = 0
    n_full_pair = 0
    for aa_apo, aa_holo in zip(aligned_apo, aligned_holo):
        if aa_apo != "-":
            apo_pos += 1
        if aa_holo != "-":
            holo_pos += 1
        if aa_apo != "-" and aa_holo != "-":
            n_matched += 1
            ha_a = HA_apo.get(apo_pos)
            ca_a = CA_apo.get(apo_pos)
            ha_h = HA_holo.get(holo_pos)
            ca_h = CA_holo.get(holo_pos)
            if (
                ha_a is not None
                and ca_a is not None
                and ha_h is not None
                and ca_h is not None
            ):
                n_full_pair += 1
    cov = (n_full_pair / n_matched) if n_matched else 0.0
    return cov, n_full_pair, n_matched


def compute_csp_multiple_saveframes_ha_ca(
    apo_sequences: List[Tuple],
    holo_sequences: List[Tuple],
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    threshold_config: Optional[Dict] = None,
    *,
    referencing_method: Optional[str] = None,
    grid_params: Optional[Dict] = None,
    target_id: Optional[str] = None,
    output_root: Optional[str] = None,
) -> List[CSPResult]:
    """Compute HA/CA CSPs using Seq_ID-offset exact AA matches."""
    import os as _os
    _verbose = (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"))

    if _verbose:
        print(f"[CSP-HA-CA] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")

    best = _best_seqid_alignment(apo_sequences, holo_sequences)
    if best is None or not best["mapping"]:
        if _verbose:
            print("[CSP-HA-CA] No valid Seq_ID alignment found for HA/CA CSPs")
        return []

    mapping = best["mapping"]
    apo_aa = {best["apo_seq_id_min"] + i: aa for i, aa in enumerate(best["apo_seq"])}
    holo_aa = {best["holo_seq_id_min"] + i: aa for i, aa in enumerate(best["holo_seq"])}
    aligned_apo = "".join(apo_aa[a] for a, _h in mapping)
    aligned_holo = "".join(holo_aa[h] for _a, h in mapping)
    HA_apo = {i: best["HA_apo"][a] for i, (a, _h) in enumerate(mapping, 1) if a in best["HA_apo"]}
    CA_apo = {i: best["CA_apo"][a] for i, (a, _h) in enumerate(mapping, 1) if a in best["CA_apo"]}
    HA_holo = {i: best["HA_holo"][h] for i, (_a, h) in enumerate(mapping, 1) if h in best["HA_holo"]}
    CA_holo = {i: best["CA_holo"][h] for i, (_a, h) in enumerate(mapping, 1) if h in best["CA_holo"]}

    min_cov = float(compute.ha_ca_min_shift_coverage)
    env_mc = os.environ.get("CSP_HA_CA_MIN_COVERAGE")
    if env_mc is not None:
        try:
            min_cov = float(env_mc)
        except ValueError:
            pass

    cov, n_pair, n_match = _ha_ca_shift_pair_coverage_for_alignment(
        aligned_apo,
        aligned_holo,
        HA_apo,
        CA_apo,
        HA_holo,
        CA_holo,
    )
    if n_match > 0 and cov < min_cov:
        tid = target_id or holo_pdb
        print(
            f"[CSP-HA-CA] Skipping HA/CA CSP for {tid}: HA+CA shift coverage on aligned "
            f"sequence is {cov:.1%} ({n_pair}/{n_match} columns), minimum {min_cov:.1%}",
            file=sys.stderr,
        )
        return []
    if _verbose and n_match > 0:
        print(
            f"[CSP-HA-CA] HA+CA shift coverage {cov:.1%} ({n_pair}/{n_match} columns), "
            f"minimum {min_cov:.1%}"
        )

    if _verbose:
        print(
            f"[CSP-HA-CA] Seq_ID-offset: apo({best['apo_saveframe']}) vs "
            f"holo({best['holo_saveframe']}) offset={best['offset']} "
            f"exact_matches={len(mapping)}"
        )

    results = compute_csp_from_aligned_sequences_ha_ca(
        aligned_apo,
        aligned_holo,
        HA_apo,
        CA_apo,
        HA_holo,
        CA_holo,
        threshold_config,
        enable_referencing=True,
        referencing_method=referencing_method,
        grid_params=grid_params,
        target_id=(target_id or holo_pdb),
        output_root=output_root,
    )
    for r in results:
        asid, hsid = mapping[r.apo_index - 1]
        r.apo_index = asid
        r.holo_index = hsid

    if _verbose:
        valid_csps = sum(1 for r in results if r.csp_A is not None)
        print(f"[CSP-HA-CA] Valid HA/CA CSPs: {valid_csps}")

    return results
