#!/usr/bin/env python3
"""
Rebuild selected apo–apo ``aggregate_csp_hsqc_grid.png`` figures: offset-grid
heatmap (top left), Query/Match BMRB conditions (top middle/right),
offset-applied HSQC (middle left), CSP along the query sequence (middle right),
and a sequence-alignment strip (bottom; same ``align_global`` pairing as
apo–apo CSP). Existing offset-grid PNGs are reused. HSQC and CSP panels are
redrawn from ``csp_table.csv`` (HSQC legend = BMRB IDs; CSP threshold 0.05 ppm).

Default pairs:
  52080_52079, 28070_28071, 34000_34001, 34394_6354, 17769_51725
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import biotite.sequence.graphics as graphics
import numpy as np

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

try:
    from .apo_holo_exp_conditions import (
        ensure_star_path,
        extract_conditions_from_star,
        extract_sample_components_from_star,
    )
    from .align import align_global
    from .bmrb_io import parse_sequence_and_shifts_from_saveframes
    from .build_apo_holo_cs_list_alignment_deck import (
        biotite_alignment_from_gapped,
        entity_names_joined,
    )
    from .case_study_3 import (
        _conditions_block_lines,
        _draw_conditions_panel,
        _trim_whitespace,
    )
    from .config import Paths
    from .HSQC_visualize import _draw_overlay_panel, format_offset_annotation
    from .visualize import plot_csp_distribution
except Exception:
    from scripts.apo_holo_exp_conditions import (  # type: ignore
        ensure_star_path,
        extract_conditions_from_star,
        extract_sample_components_from_star,
    )
    from scripts.align import align_global  # type: ignore
    from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # type: ignore
    from scripts.build_apo_holo_cs_list_alignment_deck import (  # type: ignore
        biotite_alignment_from_gapped,
        entity_names_joined,
    )
    from scripts.case_study_3 import (  # type: ignore
        _conditions_block_lines,
        _draw_conditions_panel,
        _trim_whitespace,
    )
    from scripts.config import Paths  # type: ignore
    from scripts.HSQC_visualize import (  # type: ignore
        _draw_overlay_panel,
        format_offset_annotation,
    )
    from scripts.visualize import plot_csp_distribution  # type: ignore


DEFAULT_PAIRS_ROOT = _REPO / "outputs" / "apo_apo_matches"
DEFAULT_IDS = (
    "52080_52079",
    "28070_28071",
    "34000_34001",
    "34394_6354",
    "17769_51725",
)


def _parse_pair_id(pair_id: str) -> tuple[str, str]:
    text = pair_id.strip()
    if "_" not in text:
        raise ValueError(f"Pair id must be query_match (got {pair_id!r})")
    query, match = text.split("_", 1)
    query, match = query.strip(), match.strip()
    if not query or not match:
        raise ValueError(f"Pair id must be query_match (got {pair_id!r})")
    return query, match


def _find_offset_grid_png(pair_dir: Path) -> Path:
    preferred = pair_dir / "offset_grid_H_-0.12_0.12_0.01__N_-1.2_1.2_0.05__C_0.05.png"
    if preferred.is_file():
        return preferred
    matches = sorted(pair_dir.glob("offset_grid_*.png"))
    if not matches:
        raise FileNotFoundError(f"No offset_grid_*.png in {pair_dir}")
    return matches[0]


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def _resolve_star(bmrb_id: str, cs_dir: Path, extra_dirs: Sequence[Path], do_fetch: bool) -> Path:
    for directory in (cs_dir, *extra_dirs):
        found = ensure_star_path(bmrb_id, directory, do_fetch=False)
        if found is not None:
            return found
    found = ensure_star_path(bmrb_id, cs_dir, do_fetch=do_fetch)
    if found is None:
        raise FileNotFoundError(f"NMR-STAR not found for BMRB {bmrb_id}")
    return found


def _best_global_alignment(query_sequences: Sequence, match_sequences: Sequence):
    """Same polymer pairing as ``run_apo_apo_csp.write_sequence_alignment``."""
    best_score = float("-inf")
    best: Optional[tuple] = None
    for query_entry in query_sequences:
        query_seq = query_entry[0]
        for match_entry in match_sequences:
            match_seq = match_entry[0]
            aligned_q, aligned_m, _mapping, score = align_global(query_seq, match_seq)
            if score > best_score:
                best_score = score
                best = (aligned_q, aligned_m, score)
    if best is None:
        raise RuntimeError("no global alignment between query/match CS lists")
    aligned_q, aligned_m, score = best
    return biotite_alignment_from_gapped(aligned_q, aligned_m, float(score))


def _one_side_conditions(label: str, bmrb_id: str, star_path: Path) -> List[str]:
    return _conditions_block_lines(
        label,
        bmrb_id,
        "",
        extract_conditions_from_star(str(star_path)),
        entity_names_joined(star_path),
        extract_sample_components_from_star(str(star_path)),
        wrap_width=52,
    )


def _crop_baked_title(image, *, min_major_px: int = 80, min_major_frac: float = 0.12):
    """Drop the short top title band(s); keep the first tall plot block."""
    arr = np.asarray(_trim_whitespace(image))
    if arr.ndim == 3:
        rgb = arr[..., :3]
        if arr.shape[-1] == 4:
            content = (arr[..., 3] > 0.01) & np.any(rgb < 0.985, axis=-1)
        else:
            content = np.any(rgb < 0.985, axis=-1)
    else:
        content = arr < 0.985
    is_content = content.mean(axis=1) > 0.002
    h = int(arr.shape[0])
    start = 0
    i = 0
    while i < h:
        is_c = bool(is_content[i])
        j = i + 1
        while j < h and bool(is_content[j]) == is_c:
            j += 1
        if is_c:
            run = j - i
            if run >= min_major_px or run / max(h, 1) >= min_major_frac:
                start = max(0, i - 2)
                break
        i = j
    return arr[start:]


CSP_DISPLAY_THRESHOLD_PPM = 0.05


def _csp_rows_from_table(csp_table: Path) -> List[SimpleNamespace]:
    rows: List[SimpleNamespace] = []
    with csp_table.open(newline="", encoding="utf-8") as handle:
        for rec in csv.DictReader(handle):
            raw = (rec.get("csp_A") or "").strip()
            if not raw:
                continue
            try:
                csp_val = float(raw)
                pos = int(float((rec.get("apo_resi") or "").strip()))
            except ValueError:
                continue
            aa = (rec.get("apo_aa") or "").strip() or "X"
            rows.append(SimpleNamespace(apo_index=pos, apo_aa=aa, csp_A=csp_val, significant=False))
    if not rows:
        raise FileNotFoundError(f"No recorded CSPs in {csp_table}")
    return rows


def _render_csp_panel(pair_dir: Path):
    table = _require_file(pair_dir / "csp_table.csv", "CSP table")
    rows = _csp_rows_from_table(table)
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        plot_csp_distribution(
            rows,  # type: ignore[arg-type]
            tmp_path,
            title="",
            threshold=CSP_DISPLAY_THRESHOLD_PPM,
            ylabel=r"$\Delta\delta_{NH}$ (ppm)",
        )
        return _trim_whitespace(mpimg.imread(tmp_path))
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _optional_float(raw: Optional[str]) -> Optional[float]:
    text = (raw or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _render_hsqc_offset_panel(pair_dir: Path, query_id: str, match_id: str):
    """Offset-applied HSQC overlay; legend lists BMRB IDs and applied ΔH/ΔN."""
    table = _require_file(pair_dir / "csp_table.csv", "CSP table")
    apo_pts: List[tuple] = []
    match_pts: List[tuple] = []
    segments: List[tuple] = []
    h_offset = None
    n_offset = None
    with table.open(newline="", encoding="utf-8") as handle:
        for rec in csv.DictReader(handle):
            h_a = _optional_float(rec.get("H_apo"))
            n_a = _optional_float(rec.get("N_apo"))
            h_m = _optional_float(rec.get("H_match"))
            n_m = _optional_float(rec.get("N_match"))
            if h_offset is None:
                h_offset = _optional_float(rec.get("H_offset"))
            if n_offset is None:
                n_offset = _optional_float(rec.get("N_offset"))
            if h_a is not None and n_a is not None:
                apo_pts.append((h_a, n_a))
            if h_m is not None and n_m is not None:
                match_pts.append((h_m, n_m))
            if h_a is not None and n_a is not None and h_m is not None and n_m is not None:
                segments.append(((h_a, n_a), (h_m, n_m)))
    if not apo_pts and not match_pts:
        raise FileNotFoundError(f"No H/N shifts in {table}")
    legend_title = format_offset_annotation(
        "H",
        0.0 if h_offset is None else h_offset,
        "N",
        0.0 if n_offset is None else n_offset,
    )

    axis_pts = apo_pts + match_pts
    xs = [p[0] for p in axis_pts]
    ys = [p[1] for p in axis_pts]
    x_margin = 0.05 * (max(xs) - min(xs) or 1.0)
    y_margin = 0.05 * (max(ys) - min(ys) or 1.0)
    x_limits = (max(xs) + x_margin, min(xs) - x_margin)
    y_limits = (max(ys) + y_margin, min(ys) - y_margin)

    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    ax.set_xlabel(r"$^{1}$H $\delta$ (ppm)", fontsize=12)
    ax.set_ylabel(r"$^{15}$N $\delta$ (ppm)", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
    _draw_overlay_panel(
        ax,
        apo_points=apo_pts,
        holo_points=match_pts,
        segments=segments,
        segment_colors=["#7f7f7f"] * len(segments),
        apo_label=query_id,
        holo_label=match_id,
        legend_title=legend_title,
    )
    fig.tight_layout()
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        fig.savefig(tmp_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return _trim_whitespace(mpimg.imread(tmp_path))
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def compose_aggregate(
    pair_dir: Path,
    *,
    query_id: str,
    match_id: str,
    query_conditions: Sequence[str],
    match_conditions: Sequence[str],
    alignment,
    alignment_matrix,
    output_path: Path,
) -> None:
    hsqc = _render_hsqc_offset_panel(pair_dir, query_id, match_id)
    offset = _crop_baked_title(mpimg.imread(_find_offset_grid_png(pair_dir)))
    csp = _render_csp_panel(pair_dir)

    # Narrower than the old full-width strip so the alignment stays under CSP.
    align_symbols = 50
    aln_len = len(alignment)
    n_blocks = max(1, math.ceil(aln_len / align_symbols))
    align_h_in = max(1.70, 0.865 * n_blocks + 0.45)
    n_cond = max(len(query_conditions), len(match_conditions), 1)
    grid_h_in = 6.4
    hsqc_h_in = 8.8
    fig_w = 22.0
    left_ratio, right_ratio = 1.55, 2.50
    right_w_in = fig_w * 0.97 * right_ratio / (left_ratio + right_ratio)
    # Cell height matches the CSP image so it does not letterbox (no size drop).
    csp_h_in = max(6.8, right_w_in * (csp.shape[0] / max(csp.shape[1], 1)))
    # Room for 2.5× line spacing in the BMRB conditions blocks.
    cond_h_in = max(3.1, 0.50 * n_cond + 0.45)
    left_h = grid_h_in + hsqc_h_in
    right_h = cond_h_in + csp_h_in + align_h_in
    col_h = max(left_h, right_h)
    left_spacer = max(col_h - left_h, 0.01)
    right_spacer = max(col_h - right_h, 0.01)
    fig_h = col_h + 0.85

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.02, right=0.99, top=0.94, bottom=0.025)
    outer = GridSpec(
        1,
        2,
        figure=fig,
        width_ratios=[left_ratio, right_ratio],
        wspace=0.06,
    )
    left = outer[0].subgridspec(
        3,
        1,
        height_ratios=[grid_h_in, hsqc_h_in, left_spacer],
        hspace=0.08,
    )
    right = outer[1].subgridspec(
        4,
        1,
        height_ratios=[cond_h_in, csp_h_in, align_h_in, right_spacer],
        hspace=0.06,
    )
    conds = right[0].subgridspec(1, 2, wspace=0.06)
    # Inset the alignment so Seq_ID labels do not spill into the HSQC.
    aln_slot = right[2].subgridspec(1, 2, width_ratios=[0.08, 0.92], wspace=0.0)
    ax_grid = fig.add_subplot(left[0])
    ax_hsqc = fig.add_subplot(left[1])
    ax_cond_query = fig.add_subplot(conds[0, 0])
    ax_cond_match = fig.add_subplot(conds[0, 1])
    ax_csp = fig.add_subplot(right[1])
    ax_aln = fig.add_subplot(aln_slot[0, 1])

    ax_grid.imshow(offset)
    ax_grid.axis("off")
    ax_grid.set_title(
        "Grid Search for optimal ¹H/¹⁵N Offsets",
        fontsize=13,
        fontweight="bold",
        pad=4,
    )
    _draw_conditions_panel(
        ax_cond_query,
        query_conditions,
        max_dy=0.040,
        fontsize=15,
        line_spacing=2.5,
    )
    _draw_conditions_panel(
        ax_cond_match,
        match_conditions,
        max_dy=0.040,
        fontsize=15,
        line_spacing=2.5,
    )
    ax_hsqc.imshow(hsqc)
    ax_hsqc.axis("off")
    ax_hsqc.set_title(
        f"{query_id}/{match_id} HSQC Offset",
        fontsize=13,
        fontweight="bold",
        pad=4,
    )
    ax_csp.imshow(csp)
    ax_csp.set_aspect("equal")
    ax_csp.axis("off")
    ax_csp.set_title("CSP along query apo sequence", fontsize=13, fontweight="bold", pad=4)

    graphics.plot_alignment_similarity_based(
        ax_aln,
        alignment,
        matrix=alignment_matrix,
        labels=[query_id, match_id],
        show_numbers=True,
        show_line_position=True,
        symbols_per_line=align_symbols,
    )
    ax_aln.set_title(
        "Query / Match CS-list sequence alignment",
        fontsize=13,
        fontweight="bold",
        pad=4,
    )

    fig.suptitle(
        f"{query_id} vs {match_id} — apo–apo CSP / HSQC / offset grid",
        fontsize=16,
        fontweight="bold",
        y=0.975,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def compose_pair(
    pair_id: str,
    pairs_root: Path,
    *,
    cs_dir: Path,
    do_fetch: bool,
) -> Path:
    query_id, match_id = _parse_pair_id(pair_id)
    pair_dir = pairs_root / f"{query_id}_{match_id}"
    if not pair_dir.is_dir():
        raise FileNotFoundError(f"Missing pair directory: {pair_dir}")

    extra_star_dirs = [pairs_root / "shifts"]
    query_star = _resolve_star(query_id, cs_dir, extra_star_dirs, do_fetch)
    match_star = _resolve_star(match_id, cs_dir, extra_star_dirs, do_fetch)
    query_cond = _one_side_conditions("Query", query_id, query_star)
    match_cond = _one_side_conditions("Match", match_id, match_star)

    query_sequences = parse_sequence_and_shifts_from_saveframes(str(query_star))
    match_sequences = parse_sequence_and_shifts_from_saveframes(str(match_star))
    alignment, alignment_matrix = _best_global_alignment(query_sequences, match_sequences)

    out_path = pair_dir / "aggregate_csp_hsqc_grid.png"
    compose_aggregate(
        pair_dir,
        query_id=query_id,
        match_id=match_id,
        query_conditions=query_cond,
        match_conditions=match_cond,
        alignment=alignment,
        alignment_matrix=alignment_matrix,
        output_path=out_path,
    )
    return out_path


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairs-root",
        type=Path,
        default=DEFAULT_PAIRS_ROOT,
        help="Root containing {query}_{match}/ directories (default: %(default)s).",
    )
    parser.add_argument(
        "--ids",
        type=str,
        default=",".join(DEFAULT_IDS),
        help="Comma-separated query_match pair ids (default: the five requested pairs).",
    )
    parser.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(Paths().cs_cache_dir),
        help="NMR-STAR cache directory (default: CS_Lists).",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Fetch missing NMR-STAR files from BMRB if not cached.",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(argv)
    pairs_root = args.pairs_root if args.pairs_root.is_absolute() else _REPO / args.pairs_root
    cs_dir = args.cs_dir if args.cs_dir.is_absolute() else _REPO / args.cs_dir
    pair_ids = [p.strip() for p in str(args.ids).split(",") if p.strip()]
    if not pair_ids:
        print("No pair ids provided.", file=sys.stderr)
        return 1

    failed = 0
    for pair_id in pair_ids:
        try:
            out = compose_pair(pair_id, pairs_root, cs_dir=cs_dir, do_fetch=args.fetch)
            print(f"[OK] {pair_id} -> {out}")
        except Exception as exc:
            failed += 1
            print(f"[FAIL] {pair_id}: {exc}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
