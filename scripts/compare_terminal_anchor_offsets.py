#!/usr/bin/env python3
"""Compare terminal-anchor RMSD H/N offsets vs global CSP-count grid offsets.

For selected targets, terminal residues with high RCI in both apo and holo
(RCI > 0.3, gap-of-1) are used as anchors. Analytic least-squares offsets
minimize per-nucleus RMSD on those anchors and are compared to the stored
global grid-search offsets (maximize count of CSP < 0.05 ppm).
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    _HAS_PLT = True
except Exception:
    _HAS_PLT = False

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.backfill_terminal_rci import load_state_rci, terminal_high_residues
from scripts.csp import (
    _csp_n_term,
    _floor_primary_hn_cutoff,
    compute_threshold_with_outlier_removal,
    run_offset_grid_search,
)
from scripts.config import Referencing, classification_colors, thresholds as csp_thresholds
from scripts.merge_csv import exceeds_max_classification_csp_z

DEFAULT_TARGETS = (
    "2KZU_17018",
    "1YWI_6558",
    "2RR4_11362",
    "2K17_15670",
)

DEFAULT_GRID_NAME = "offset_grid_H_-0.2_0.2_0.01__N_-1.5_1.5_0.05__C_0.05.csv"

_APO_COLOR = "#7f7f7f"
_HOLO_COLOR = "#e74c3c"
_ANCHOR_EDGE = "limegreen"
_LINE_COLOR = "#555555"


def _parse_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def _holo_hn(row: Dict[str, str]) -> Tuple[Optional[float], Optional[float]]:
    h = _parse_float(row.get("H_holo_original"))
    if h is None:
        h = _parse_float(row.get("H_holo"))
    n = _parse_float(row.get("N_holo_original"))
    if n is None:
        n = _parse_float(row.get("N_holo"))
    return h, n


def load_global_offsets(
    tgt_dir: Path,
    rows: Sequence[Dict[str, str]],
) -> Tuple[float, float, int]:
    """Return (H_offset, N_offset, best_count) from stored grid or recompute."""
    grid_path = tgt_dir / DEFAULT_GRID_NAME
    if grid_path.is_file():
        best_h = best_n = None
        best_count = 0
        for line in grid_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not line.startswith("best_"):
                continue
            parts = line.split(",", 1)
            if len(parts) != 2:
                continue
            key, raw = parts[0].strip(), parts[1].strip()
            if key == "best_h_offset":
                best_h = float(raw)
            elif key == "best_n_offset":
                best_n = float(raw)
            elif key == "best_count":
                best_count = int(float(raw))
        if best_h is not None and best_n is not None:
            return best_h, best_n, best_count

    # Fallback: rebuild minimal CSPResult-like objects for grid search
    from scripts.csp import CSPResult

    results: List[CSPResult] = []
    for row in rows:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        apo_resi = int(float(row["apo_resi"]))
        holo_resi = int(float(row["holo_resi"]))
        results.append(
            CSPResult(
                apo_index=apo_resi,
                holo_index=holo_resi,
                apo_aa=(row.get("apo_aa") or "").strip() or "X",
                holo_aa=(row.get("holo_aa") or "").strip() or "X",
                H_apo=h_apo,
                N_apo=n_apo,
                H_holo=h_holo,
                N_holo=n_holo,
                dH=None,
                dN=None,
                csp_A=None,
                significant=False,
                H_holo_original=h_holo,
                N_holo_original=n_holo,
            )
        )
    ref = Referencing()
    grid = run_offset_grid_search(
        results,
        h_min=ref.grid_h_min,
        h_max=ref.grid_h_max,
        h_step=ref.grid_h_step,
        n_min=ref.grid_n_min,
        n_max=ref.grid_n_max,
        n_step=ref.grid_n_step,
        cutoff=ref.grid_cutoff,
    )
    return (
        float(grid["best_h_offset"]),
        float(grid["best_n_offset"]),
        int(grid["best_count"]),
    )


def both_terminal_anchor_rows(
    tgt_dir: Path,
    rows: Sequence[Dict[str, str]],
    *,
    threshold: float = 0.3,
    max_gap: int = 1,
) -> List[Dict[str, str]]:
    """Master/csp rows that are terminal-high in both apo and holo with usable HN."""
    # Prefer csp_table for shift completeness; fall back to master rows.
    csp_path = tgt_dir / "csp_table.csv"
    shift_rows = _read_csv_rows(csp_path) if csp_path.is_file() else list(rows)
    master_rows = list(rows)

    apo_flags = terminal_high_residues(
        load_state_rci(tgt_dir, "apo", master_rows),
        threshold=threshold,
        max_gap=max_gap,
    )
    holo_flags = terminal_high_residues(
        load_state_rci(tgt_dir, "holo", master_rows),
        threshold=threshold,
        max_gap=max_gap,
    )

    anchors: List[Dict[str, str]] = []
    for row in shift_rows:
        try:
            apo_resi = int(float(row["apo_resi"]))
            holo_resi = int(float(row["holo_resi"]))
        except (KeyError, TypeError, ValueError):
            continue
        if apo_resi not in apo_flags or holo_resi not in holo_flags:
            continue
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        anchors.append(row)
    return anchors


def terminal_anchor_offsets(
    anchors: Sequence[Dict[str, str]],
) -> Tuple[float, float, float]:
    """Return (H_offset, N_offset, post-fit CSP RMSD) for anchor residues."""
    if not anchors:
        return 0.0, 0.0, float("nan")

    h_diffs: List[float] = []
    n_diffs: List[float] = []
    for row in anchors:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        assert h_apo is not None and n_apo is not None
        assert h_holo is not None and n_holo is not None
        h_diffs.append(h_apo - h_holo)
        n_diffs.append(n_apo - n_holo)

    h_off = statistics.mean(h_diffs)
    n_off = statistics.mean(n_diffs)

    csp_sq: List[float] = []
    for row in anchors:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        assert h_apo is not None and n_apo is not None
        assert h_holo is not None and n_holo is not None
        dH = (h_holo + h_off) - h_apo
        dN = (n_holo + n_off) - n_apo
        csp_sq.append(0.5 * (dH * dH + _csp_n_term(dN)))
    rmsd = math.sqrt(statistics.mean(csp_sq)) if csp_sq else float("nan")
    return h_off, n_off, rmsd


def count_aligned_hn(
    rows: Sequence[Dict[str, str]],
    h_offset: float,
    n_offset: float,
    *,
    cutoff: float = 0.05,
) -> Tuple[int, int]:
    """Count HN pairs with CSP < cutoff after applying holo offsets.

    Same definition as ``run_offset_grid_search``: complete H/N in both states,
    CSP = sqrt(0.5*(dH^2 + (dN/5)^2)). Returns (aligned_count, n_evaluable).
    """
    aligned = 0
    evaluable = 0
    for row in rows:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        evaluable += 1
        dH = (h_holo + h_offset) - h_apo
        dN = (n_holo + n_offset) - n_apo
        csp_val = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
        if csp_val < cutoff:
            aligned += 1
    return aligned, evaluable


def process_target(tgt_dir: Path) -> Dict[str, object]:
    master_path = tgt_dir / "master_alignment.csv"
    if not master_path.is_file():
        raise FileNotFoundError(f"missing {master_path}")
    master_rows = _read_csv_rows(master_path)
    anchors = both_terminal_anchor_rows(tgt_dir, master_rows)
    h_anchor, n_anchor, anchor_rmsd = terminal_anchor_offsets(anchors)

    csp_path = tgt_dir / "csp_table.csv"
    shift_rows = _read_csv_rows(csp_path) if csp_path.is_file() else master_rows
    h_global, n_global, best_count = load_global_offsets(tgt_dir, shift_rows)

    cutoff = float(Referencing().grid_cutoff)
    aligned_anchor, n_hn = count_aligned_hn(shift_rows, h_anchor, n_anchor, cutoff=cutoff)
    aligned_global, _ = count_aligned_hn(shift_rows, h_global, n_global, cutoff=cutoff)

    return {
        "target": tgt_dir.name,
        "tgt_dir": tgt_dir,
        "shift_rows": shift_rows,
        "anchors": anchors,
        "n_anchors": len(anchors),
        "n_hn_pairs": n_hn,
        "H_offset_anchor": h_anchor,
        "N_offset_anchor": n_anchor,
        "anchor_csp_rmsd": anchor_rmsd,
        "aligned_count_anchor": aligned_anchor,
        "H_offset_global": h_global,
        "N_offset_global": n_global,
        "global_best_count": best_count,
        "aligned_count_global": aligned_global,
        "delta_H": h_anchor - h_global,
        "delta_N": n_anchor - n_global,
        "cutoff": cutoff,
    }


def _hn_pairs(
    rows: Sequence[Dict[str, str]],
) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Return [(apo_HN, holo_HN_orig), ...] with H on x and N on y."""
    pairs: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    for row in rows:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        pairs.append(((h_apo, n_apo), (h_holo, n_holo)))
    return pairs


def _parse_bool01(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def _axis_limits(
    points: Sequence[Tuple[float, float]],
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Descending ppm limits (high left / high bottom) with a small margin."""
    if not points:
        return (1.0, 0.0), (1.0, 0.0)
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_margin = 0.05 * (max(xs) - min(xs) or 1.0)
    y_margin = 0.05 * (max(ys) - min(ys) or 1.0)
    return (
        (max(xs) + x_margin, min(xs) - x_margin),
        (max(ys) + y_margin, min(ys) - y_margin),
    )


def compute_csp_series(
    rows: Sequence[Dict[str, str]],
    h_offset: float,
    n_offset: float,
) -> Tuple[List[Dict[str, object]], float]:
    """Recompute HN CSPs and significance for a given holo offset scheme.

    Returns (entries, threshold) where each entry has holo_resi, holo_aa, csp,
    significant, occluded, and classification (TP/FP/TN/FN).
    """
    raw: List[Dict[str, object]] = []
    values: List[float] = []
    for row in rows:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        try:
            holo_resi = int(float(row["holo_resi"]))
        except (KeyError, TypeError, ValueError):
            continue
        dH = (h_holo + h_offset) - h_apo
        dN = (n_holo + n_offset) - n_apo
        csp = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
        values.append(csp)
        raw.append(
            {
                "holo_resi": holo_resi,
                "holo_aa": (row.get("holo_aa") or "").strip() or "X",
                "csp": csp,
                "occluded": _parse_bool01(row.get("occluded")),
            }
        )

    thr_info = compute_threshold_with_outlier_removal(
        values,
        float(csp_thresholds.outlier_z_score),
        float(csp_thresholds.significance_z_score),
        int(csp_thresholds.max_outlier_iterations),
        float(csp_thresholds.max_outlier_fraction),
    )
    # Match primary pipeline HN cutoff: max(cleaned mean (+ z·SD), 0.05 ppm).
    threshold = _floor_primary_hn_cutoff(float(thr_info.threshold))
    cleaned_mean = float(thr_info.mean)
    cleaned_sd = float(thr_info.sd)

    entries: List[Dict[str, object]] = []
    for item in raw:
        csp_val = float(item["csp"])
        csp_z = (
            (csp_val - cleaned_mean) / cleaned_sd
            if cleaned_sd > 0.0
            else 0.0
        )
        # Match pipeline: extreme z-scores are left unclassified (omit bar).
        if exceeds_max_classification_csp_z(csp_z):
            continue
        is_sig = csp_val >= threshold
        is_bind = bool(item["occluded"])
        if is_sig and is_bind:
            cls = "TP"
        elif is_sig and not is_bind:
            cls = "FP"
        elif (not is_sig) and (not is_bind):
            cls = "TN"
        else:
            cls = "FN"
        entries.append(
            {
                **item,
                "significant": is_sig,
                "classification": cls,
            }
        )
    return entries, threshold


def _anchor_holo_resi(anchors: Sequence[Dict[str, str]]) -> set:
    out = set()
    for row in anchors:
        try:
            out.add(int(float(row["holo_resi"])))
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _draw_overlay_scheme(
    ax,
    *,
    pairs: Sequence[Tuple[Tuple[float, float], Tuple[float, float]]],
    h_offset: float,
    n_offset: float,
    aligned_count: int,
    n_hn: int,
    scheme_title: str,
    highlight_pairs: Optional[Sequence[Tuple[Tuple[float, float], Tuple[float, float]]]] = None,
    x_limits: Tuple[float, float],
    y_limits: Tuple[float, float],
    legend_fontsize: float = 8,
    anchor_label: str = "Terminal anchors",
) -> None:
    """Draw apo/holo overlay with applied holo offsets on one axes."""
    apo_pts = [apo for apo, _ in pairs]
    holo_pts = [(h[0] + h_offset, h[1] + n_offset) for _, h in pairs]
    segments = list(zip(apo_pts, holo_pts))

    if segments:
        ax.add_collection(
            LineCollection(
                segments,
                colors=_LINE_COLOR,
                linewidths=0.7,
                alpha=0.45,
                zorder=1,
            )
        )
    if apo_pts:
        ax.scatter(
            [p[0] for p in apo_pts],
            [p[1] for p in apo_pts],
            s=18,
            c=_APO_COLOR,
            alpha=0.85,
            edgecolors="none",
            label="Apo",
            zorder=2,
        )
    if holo_pts:
        ax.scatter(
            [p[0] for p in holo_pts],
            [p[1] for p in holo_pts],
            s=18,
            c=_HOLO_COLOR,
            alpha=0.85,
            edgecolors="none",
            label="Holo (offset)",
            zorder=2,
        )

    if highlight_pairs:
        h_apo = [apo for apo, _ in highlight_pairs]
        h_holo = [(h[0] + h_offset, h[1] + n_offset) for _, h in highlight_pairs]
        for pts in (h_apo, h_holo):
            if not pts:
                continue
            ax.scatter(
                [p[0] for p in pts],
                [p[1] for p in pts],
                s=90,
                facecolors="none",
                edgecolors=_ANCHOR_EDGE,
                linewidths=1.6,
                zorder=3,
            )
        ax.scatter(
            [],
            [],
            s=90,
            facecolors="none",
            edgecolors=_ANCHOR_EDGE,
            linewidths=1.6,
            label=f"{anchor_label} (n={len(highlight_pairs)})",
        )

    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xlabel(r"$^{1}$H $\delta$ (ppm)", fontsize=11)
    ax.set_ylabel(r"$^{15}$N $\delta$ (ppm)", fontsize=11)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
    ax.set_title(scheme_title, fontsize=11)

    legend = ax.legend(
        title=(
            f"ΔH={h_offset:+.3f} ppm, ΔN={n_offset:+.3f} ppm\n"
            f"Aligned CSP<0.05: {aligned_count}/{n_hn}"
        ),
        fontsize=legend_fontsize,
        title_fontsize=legend_fontsize,
        loc="best",
        framealpha=0.9,
    )
    if legend is not None:
        legend.set_title(legend.get_title().get_text())


def _draw_csp_classification_bars(
    ax,
    *,
    entries: Sequence[Dict[str, object]],
    threshold: float,
    h_offset: float,
    n_offset: float,
    aligned_count: int,
    n_hn: int,
    scheme_title: str,
    mark_anchor_resi: Optional[set] = None,
    legend_fontsize: float = 7,
    anchor_label: str = "Terminal anchors",
) -> None:
    """Per-residue CSP bars colored by TP/FP/TN/FN; optional green stars on anchors."""
    if not entries:
        ax.text(0.5, 0.5, "No CSP values", transform=ax.transAxes, ha="center", va="center")
        ax.set_title(scheme_title, fontsize=11)
        return

    residue_numbers = [int(e["holo_resi"]) for e in entries]
    csp_values = [float(e["csp"]) for e in entries]
    colors = {
        "TP": classification_colors.TP,
        "FP": classification_colors.FP,
        "TN": classification_colors.TN,
        "FN": classification_colors.FN,
    }
    bar_colors = [colors[str(e["classification"])] for e in entries]
    ax.bar(
        residue_numbers,
        csp_values,
        color=bar_colors,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.4,
        width=0.85,
    )
    ax.axhline(
        y=threshold,
        color="black",
        linestyle="--",
        linewidth=1.2,
        alpha=0.85,
        label=f"Threshold = {threshold:.3f}",
    )

    if mark_anchor_resi:
        star_x: List[float] = []
        star_y: List[float] = []
        y_pad = 0.04 * (max(csp_values) if csp_values else 1.0)
        for resi, csp in zip(residue_numbers, csp_values):
            if resi in mark_anchor_resi:
                star_x.append(resi)
                star_y.append(csp + y_pad)
        if star_x:
            ax.scatter(
                star_x,
                star_y,
                marker="*",
                s=70,
                c=_ANCHOR_EDGE,
                zorder=4,
                label=f"{anchor_label} (n={len(star_x)})",
            )

    counts = {k: sum(1 for e in entries if e["classification"] == k) for k in ("TP", "FP", "TN", "FN")}
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=colors["TP"], alpha=0.85, label=f"TP ({counts['TP']})"),
        Rectangle((0, 0), 1, 1, facecolor=colors["FP"], alpha=0.85, label=f"FP ({counts['FP']})"),
        Rectangle((0, 0), 1, 1, facecolor=colors["TN"], alpha=0.85, label=f"TN ({counts['TN']})"),
        Rectangle((0, 0), 1, 1, facecolor=colors["FN"], alpha=0.85, label=f"FN ({counts['FN']})"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1.2, label=f"Thr={threshold:.3f}"),
    ]
    if mark_anchor_resi:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="*",
                color="none",
                markerfacecolor=_ANCHOR_EDGE,
                markersize=10,
                label=anchor_label,
            )
        )

    ax.legend(
        handles=handles,
        title=(
            f"ΔH={h_offset:+.3f} ppm, ΔN={n_offset:+.3f} ppm\n"
            f"Aligned CSP<0.05: {aligned_count}/{n_hn}"
        ),
        fontsize=legend_fontsize,
        title_fontsize=legend_fontsize,
        loc="best",
        framealpha=0.9,
    )
    ax.set_xlabel("Holo residue", fontsize=10)
    ax.set_ylabel("CSP (ppm)", fontsize=10)
    ax.set_title(scheme_title, fontsize=11)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    if csp_values:
        ymax = max(max(csp_values), threshold) * 1.18
        ax.set_ylim(0.0, ymax if ymax > 0 else 0.1)


def _scheme_plot_data(result: Dict[str, object]) -> Dict[str, object]:
    """Shared per-target arrays for HSQC + CSP panels."""
    shift_rows = result["shift_rows"]  # type: ignore[assignment]
    anchors = result["anchors"]  # type: ignore[assignment]
    pairs = _hn_pairs(shift_rows)
    highlight = _hn_pairs(anchors)
    h_a = float(result["H_offset_anchor"])
    n_a = float(result["N_offset_anchor"])
    h_g = float(result["H_offset_global"])
    n_g = float(result["N_offset_global"])

    axis_pts: List[Tuple[float, float]] = []
    for apo, holo in pairs:
        axis_pts.append(apo)
        axis_pts.append((holo[0] + h_a, holo[1] + n_a))
        axis_pts.append((holo[0] + h_g, holo[1] + n_g))
    x_limits, y_limits = _axis_limits(axis_pts)

    entries_a, thr_a = compute_csp_series(shift_rows, h_a, n_a)
    entries_g, thr_g = compute_csp_series(shift_rows, h_g, n_g)
    return {
        "pairs": pairs,
        "highlight": highlight,
        "h_a": h_a,
        "n_a": n_a,
        "h_g": h_g,
        "n_g": n_g,
        "x_limits": x_limits,
        "y_limits": y_limits,
        "entries_a": entries_a,
        "thr_a": thr_a,
        "entries_g": entries_g,
        "thr_g": thr_g,
        "anchor_resi": _anchor_holo_resi(anchors),
    }


def plot_target_overlay(
    result: Dict[str, object],
    out_path: Path,
) -> None:
    """Write a 2x2 figure: HSQC overlays (top) + CSP classification bars (bottom)."""
    if not _HAS_PLT:
        raise RuntimeError("matplotlib is required to generate HSQC plots")

    data = _scheme_plot_data(result)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    _draw_overlay_scheme(
        axes[0, 0],
        pairs=data["pairs"],
        h_offset=float(data["h_a"]),
        n_offset=float(data["n_a"]),
        aligned_count=int(result["aligned_count_anchor"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Terminal-anchor RMSD — HSQC",
        highlight_pairs=data["highlight"],
        x_limits=data["x_limits"],
        y_limits=data["y_limits"],
    )
    _draw_overlay_scheme(
        axes[0, 1],
        pairs=data["pairs"],
        h_offset=float(data["h_g"]),
        n_offset=float(data["n_g"]),
        aligned_count=int(result["aligned_count_global"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Global CSP-count grid — HSQC",
        highlight_pairs=None,
        x_limits=data["x_limits"],
        y_limits=data["y_limits"],
    )
    _draw_csp_classification_bars(
        axes[1, 0],
        entries=data["entries_a"],
        threshold=float(data["thr_a"]),
        h_offset=float(data["h_a"]),
        n_offset=float(data["n_a"]),
        aligned_count=int(result["aligned_count_anchor"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Terminal-anchor RMSD — CSP bars",
        mark_anchor_resi=data["anchor_resi"],
    )
    _draw_csp_classification_bars(
        axes[1, 1],
        entries=data["entries_g"],
        threshold=float(data["thr_g"]),
        h_offset=float(data["h_g"]),
        n_offset=float(data["n_g"]),
        aligned_count=int(result["aligned_count_global"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Global CSP-count grid — CSP bars",
        mark_anchor_resi=None,
    )

    fig.suptitle(str(result["target"]), fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_aggregate_overlay(
    results: Sequence[Dict[str, object]],
    out_path: Path,
) -> None:
    """Write a rows×4 summary: HSQC_anchor | HSQC_global | CSP_anchor | CSP_global."""
    if not _HAS_PLT:
        raise RuntimeError("matplotlib is required to generate HSQC plots")
    if not results:
        return

    n_rows = len(results)
    fig, axes = plt.subplots(n_rows, 4, figsize=(22, 4.8 * n_rows), squeeze=False)

    for row_idx, result in enumerate(results):
        data = _scheme_plot_data(result)
        target = str(result["target"])

        _draw_overlay_scheme(
            axes[row_idx, 0],
            pairs=data["pairs"],
            h_offset=float(data["h_a"]),
            n_offset=float(data["n_a"]),
            aligned_count=int(result["aligned_count_anchor"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Anchor HSQC",
            highlight_pairs=data["highlight"],
            x_limits=data["x_limits"],
            y_limits=data["y_limits"],
            legend_fontsize=6,
        )
        _draw_overlay_scheme(
            axes[row_idx, 1],
            pairs=data["pairs"],
            h_offset=float(data["h_g"]),
            n_offset=float(data["n_g"]),
            aligned_count=int(result["aligned_count_global"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Global HSQC",
            highlight_pairs=None,
            x_limits=data["x_limits"],
            y_limits=data["y_limits"],
            legend_fontsize=6,
        )
        _draw_csp_classification_bars(
            axes[row_idx, 2],
            entries=data["entries_a"],
            threshold=float(data["thr_a"]),
            h_offset=float(data["h_a"]),
            n_offset=float(data["n_a"]),
            aligned_count=int(result["aligned_count_anchor"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Anchor CSP bars",
            mark_anchor_resi=data["anchor_resi"],
            legend_fontsize=6,
        )
        _draw_csp_classification_bars(
            axes[row_idx, 3],
            entries=data["entries_g"],
            threshold=float(data["thr_g"]),
            h_offset=float(data["h_g"]),
            n_offset=float(data["n_g"]),
            aligned_count=int(result["aligned_count_global"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Global CSP bars",
            mark_anchor_resi=None,
            legend_fontsize=6,
        )

        # Row label (a., b., c., …) at top-left of each row; leave titles unchanged.
        if row_idx < 26:
            axes[row_idx, 0].text(
                -0.12,
                1.08,
                f"{chr(ord('a') + row_idx)}.",
                transform=axes[row_idx, 0].transAxes,
                ha="left",
                va="bottom",
                fontsize=16,
                fontweight="bold",
                color="black",
                clip_on=False,
            )

    fig.suptitle(
        "Terminal-anchor vs global offsets: HSQC overlays and CSP classification bars",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_ROOT / "outputs",
        help="Per-target output root",
    )
    parser.add_argument(
        "--ids",
        nargs="*",
        default=list(DEFAULT_TARGETS),
        help="Target directory names (default: the four ≥5 both-terminal targets)",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Output CSV path (default: <outputs>/terminal_anchor_offset_compare.csv)",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip HSQC overlay figure generation",
    )
    args = parser.parse_args(argv)
    out_csv = args.out_csv or (args.outputs / "terminal_anchor_offset_compare.csv")

    results: List[Dict[str, object]] = []
    failures: List[Tuple[str, str]] = []
    for name in args.ids:
        tgt = args.outputs / name
        try:
            results.append(process_target(tgt))
        except Exception as e:
            failures.append((name, str(e)))
            print(f"FAIL {name}: {e}", flush=True)

    fieldnames = [
        "target",
        "n_anchors",
        "n_hn_pairs",
        "H_offset_anchor",
        "N_offset_anchor",
        "anchor_csp_rmsd",
        "aligned_count_anchor",
        "H_offset_global",
        "N_offset_global",
        "global_best_count",
        "aligned_count_global",
        "delta_H",
        "delta_N",
    ]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(
                {
                    "target": row["target"],
                    "n_anchors": row["n_anchors"],
                    "n_hn_pairs": row["n_hn_pairs"],
                    "H_offset_anchor": f"{float(row['H_offset_anchor']):.6f}",
                    "N_offset_anchor": f"{float(row['N_offset_anchor']):.6f}",
                    "anchor_csp_rmsd": f"{float(row['anchor_csp_rmsd']):.6f}",
                    "aligned_count_anchor": row["aligned_count_anchor"],
                    "H_offset_global": f"{float(row['H_offset_global']):.6f}",
                    "N_offset_global": f"{float(row['N_offset_global']):.6f}",
                    "global_best_count": row["global_best_count"],
                    "aligned_count_global": row["aligned_count_global"],
                    "delta_H": f"{float(row['delta_H']):.6f}",
                    "delta_N": f"{float(row['delta_N']):.6f}",
                }
            )

    # Console table
    hdr = (
        f"{'target':<14} {'nA':>3} {'nHN':>4} "
        f"{'H_anc':>8} {'N_anc':>8} {'rmsd':>8} {'alnA':>5} "
        f"{'H_glb':>8} {'N_glb':>8} {'alnG':>5} {'dH':>8} {'dN':>8}"
    )
    print(hdr, flush=True)
    print("-" * len(hdr), flush=True)
    for row in results:
        print(
            f"{row['target']:<14} {row['n_anchors']:>3} {row['n_hn_pairs']:>4} "
            f"{float(row['H_offset_anchor']):>8.4f} {float(row['N_offset_anchor']):>8.4f} "
            f"{float(row['anchor_csp_rmsd']):>8.4f} {row['aligned_count_anchor']:>5} "
            f"{float(row['H_offset_global']):>8.4f} {float(row['N_offset_global']):>8.4f} "
            f"{row['aligned_count_global']:>5} "
            f"{float(row['delta_H']):>8.4f} {float(row['delta_N']):>8.4f}",
            flush=True,
        )
    tot_a = sum(int(r["aligned_count_anchor"]) for r in results)
    tot_g = sum(int(r["aligned_count_global"]) for r in results)
    tot_hn = sum(int(r["n_hn_pairs"]) for r in results)
    print(
        f"\nTotals (CSP < {float(Referencing().grid_cutoff):g}): "
        f"anchor={tot_a}/{tot_hn}  global={tot_g}/{tot_hn}",
        flush=True,
    )
    print(f"Wrote {out_csv}", flush=True)

    if not args.no_plots:
        if not _HAS_PLT:
            print("WARNING: matplotlib unavailable; skipping plots", flush=True)
        else:
            for row in results:
                out_png = Path(row["tgt_dir"]) / "hsqc_overlay_anchor_vs_global.png"  # type: ignore[arg-type]
                plot_target_overlay(row, out_png)
                print(f"Wrote {out_png}", flush=True)
            agg_png = args.outputs / "hsqc_overlay_anchor_vs_global_all.png"
            plot_aggregate_overlay(results, agg_png)
            print(f"Wrote {agg_png}", flush=True)

    if failures:
        print(f"Failures: {len(failures)}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
