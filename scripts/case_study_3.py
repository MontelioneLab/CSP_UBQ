#!/usr/bin/env python3
"""
Case-study v3 layout: case_study (v1) panels plus top-left apo/holo conditions
and a bottom Seq_ID-offset sequence alignment strip.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import biotite.sequence.graphics as graphics
from pptx import Presentation
from pptx.util import Inches

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

try:
    from .apo_holo_exp_conditions import (
        ExtractedConditions,
        ensure_star_path,
        extract_conditions_from_star,
        extract_sample_components_from_star,
    )
    from .bmrb_io import parse_sequence_and_shifts_from_saveframes
    from .build_apo_holo_cs_list_alignment_deck import (
        best_seqid_pair_alignment,
        entity_names_joined,
        add_image_slide,
    )
    from .case_study import (
        format_case_study_metadata_header,
        format_hsqc_offset_panel_title,
        extract_hsqc_bottom_right_panel,
        render_pymol_panel_with_view,
        capture_user_view_interactive,
        _case_study_view_load_candidates,
        _case_study_view_save_path,
        _load_view,
    )
    from .config import Paths, paths
except Exception:
    from scripts.apo_holo_exp_conditions import (  # type: ignore
        ExtractedConditions,
        ensure_star_path,
        extract_conditions_from_star,
        extract_sample_components_from_star,
    )
    from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # type: ignore
    from scripts.build_apo_holo_cs_list_alignment_deck import (  # type: ignore
        best_seqid_pair_alignment,
        entity_names_joined,
        add_image_slide,
    )
    from scripts.case_study import (  # type: ignore
        format_case_study_metadata_header,
        format_hsqc_offset_panel_title,
        extract_hsqc_bottom_right_panel,
        render_pymol_panel_with_view,
        capture_user_view_interactive,
        _case_study_view_load_candidates,
        _case_study_view_save_path,
        _load_view,
    )
    from scripts.config import Paths, paths  # type: ignore

_SYMBOLS_PER_LINE = 60


def _trim_whitespace(image):
    arr = np.asarray(image)
    if arr.ndim == 2:
        content_mask = arr < 0.985
    else:
        rgb = arr[..., :3]
        if arr.shape[-1] == 4:
            alpha = arr[..., 3]
            content_mask = (alpha > 0.01) & (np.any(rgb < 0.985, axis=-1))
        else:
            content_mask = np.any(rgb < 0.985, axis=-1)

    rows = np.where(content_mask.any(axis=1))[0]
    cols = np.where(content_mask.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return arr

    pad = 4
    r0 = max(0, int(rows[0]) - pad)
    r1 = min(arr.shape[0], int(rows[-1]) + pad + 1)
    c0 = max(0, int(cols[0]) - pad)
    c1 = min(arr.shape[1], int(cols[-1]) + pad + 1)
    return arr[r0:r1, c0:c1]


def _fmt_num(val: Optional[float], digits: int) -> str:
    if val is None:
        return ""
    return f"{val:.{digits}f}".rstrip("0").rstrip(".")


def _conditions_block_lines(
    label: str,
    bmrb: str,
    pdb: str,
    ec: ExtractedConditions,
    entities: str,
    components: str,
    wrap_width: int = 44,
) -> List[str]:
    press = ""
    if ec.pressure is not None:
        units = ec.pressure_units or ""
        press = f"  pressure={_fmt_num(ec.pressure, 3)}{(' ' + units) if units else ''}"
    ph = _fmt_num(ec.pH, 3) or "n/a"
    temp = _fmt_num(ec.temperature_C, 2) or "n/a"
    ent = entities if entities else "(none)"
    comp = components if components else "(none)"
    lines = [
        f"{label} BMRB {bmrb} | PDB {(pdb or '').upper() or 'N/A'}",
        f"  pH={ph}  T={temp} °C{press}",
        f"  Entities: {ent}",
        f"  Components: {comp}",
    ]
    wrapped: List[str] = []
    for ln in lines:
        if ln.startswith("  Components:") or ln.startswith("  Entities:"):
            wrapped.extend(textwrap.wrap(ln, width=wrap_width, subsequent_indent="    ") or [ln])
        else:
            wrapped.append(ln)
    return wrapped


def build_conditions_text(
    *,
    apo_bmrb: str,
    holo_bmrb: str,
    apo_pdb: str,
    holo_pdb: str,
    apo_ec: ExtractedConditions,
    holo_ec: ExtractedConditions,
    apo_entities: str,
    holo_entities: str,
    apo_components: str,
    holo_components: str,
    pair_meta: Optional[str] = None,
) -> List[str]:
    lines: List[str] = []
    if pair_meta:
        lines.append(pair_meta)
        lines.append("")
    lines.extend(
        _conditions_block_lines("Apo", apo_bmrb, apo_pdb, apo_ec, apo_entities, apo_components)
    )
    lines.append("")
    lines.extend(
        _conditions_block_lines("Holo", holo_bmrb, holo_pdb, holo_ec, holo_entities, holo_components)
    )
    return lines


def _draw_conditions_panel(
    ax,
    conditions_lines: Sequence[str],
    *,
    max_dy: float = 0.10,
    fontsize: float = 13,
    line_spacing: float = 1.0,
) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    n_lines = max(len(conditions_lines), 1)
    spacing = max(float(line_spacing), 1.0)
    slots = 1.0 + (n_lines - 1) * spacing
    dy = min(max_dy, 0.95 / slots)
    y = 0.98
    for i, ln in enumerate(conditions_lines):
        weight = (
            "bold"
            if ln.startswith(("Apo ", "Holo ", "Query ", "Match "))
            else "normal"
        )
        ax.text(
            0.02,
            y - i * dy * spacing,
            ln,
            fontsize=fontsize,
            fontfamily="monospace",
            fontweight=weight,
            va="top",
            ha="left",
            transform=ax.transAxes,
            linespacing=1.05,
        )


def _add_panel_letter(axis, label: str) -> None:
    axis.text(
        0.015,
        0.985,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="black",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 2},
    )


def compose_case_study_3_figure(
    hsqc_bottom_right_panel,
    csp_classification_plot_path: str,
    rendered_panels: Dict[str, str],
    output_path: str,
    *,
    conditions_lines: Sequence[str],
    alignment,
    alignment_matrix,
    case_study_header: Optional[str] = None,
    hsqc_panel_title: Optional[str] = None,
) -> None:
    """
    Build case-study-3 figure based on v1 layout:
      - top: conditions | HSQC | CSP classification bars
      - middle: three PyMOL renders
      - bottom: Seq_ID apo/holo alignment
    """
    csp_plot = _trim_whitespace(mpimg.imread(csp_classification_plot_path))
    panel_left = _trim_whitespace(mpimg.imread(rendered_panels["left"]))
    panel_middle = _trim_whitespace(mpimg.imread(rendered_panels["middle"]))
    panel_right = _trim_whitespace(mpimg.imread(rendered_panels["right"]))
    hsqc_panel = _trim_whitespace(hsqc_bottom_right_panel)

    aln_len = len(alignment)
    n_blocks = max(1, math.ceil(aln_len / _SYMBOLS_PER_LINE))
    align_h_in = max(2.0, 1.05 * n_blocks + 0.45)
    main_h_in = 11.0
    fig_h = main_h_in + align_h_in + 0.8

    fig = plt.figure(figsize=(24, fig_h))
    outer = GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[main_h_in, align_h_in],
        hspace=0.10,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.03)
    if case_study_header:
        fig.suptitle(case_study_header, fontsize=18, fontweight="bold", x=0.5, y=0.985)

    main = outer[0].subgridspec(
        2,
        1,
        height_ratios=[1.0, 1.15],
        hspace=0.08,
    )
    top = main[0].subgridspec(
        1,
        3,
        wspace=0.02,
        width_ratios=[1.6, 2.2, 3.2],
    )
    mid = main[1].subgridspec(1, 3, wspace=0.02)

    ax_cond = fig.add_subplot(top[0, 0])
    ax_hsqc = fig.add_subplot(top[0, 1])
    ax_csp = fig.add_subplot(top[0, 2])
    ax_bottom_left = fig.add_subplot(mid[0, 0])
    ax_bottom_mid = fig.add_subplot(mid[0, 1])
    ax_bottom_right = fig.add_subplot(mid[0, 2])
    ax_aln = fig.add_subplot(outer[1])

    _draw_conditions_panel(ax_cond, conditions_lines)

    ax_hsqc.set_box_aspect(hsqc_panel.shape[0] / hsqc_panel.shape[1])
    ax_csp.set_box_aspect(csp_plot.shape[0] / csp_plot.shape[1])
    ax_bottom_left.set_box_aspect(panel_left.shape[0] / panel_left.shape[1])
    ax_bottom_mid.set_box_aspect(panel_middle.shape[0] / panel_middle.shape[1])
    ax_bottom_right.set_box_aspect(panel_right.shape[0] / panel_right.shape[1])

    ax_hsqc.imshow(hsqc_panel)
    ax_hsqc.set_title(
        hsqc_panel_title or format_hsqc_offset_panel_title(),
        fontsize=14,
        fontweight="bold",
        pad=4,
    )
    ax_csp.imshow(csp_plot)
    ax_csp.set_title("CSP Classification", fontsize=14, fontweight="bold", pad=4)
    ax_csp.set_xlabel("Sequence", fontsize=14)
    ax_csp.set_ylabel(r"$\Delta\delta_{NH}$ (ppm)", fontsize=14)
    ax_csp.set_xticks([])
    ax_csp.set_yticks([])
    for spine in ax_csp.spines.values():
        spine.set_visible(False)

    ax_bottom_left.imshow(panel_left)
    ax_bottom_left.set_title("CSP Mask", fontsize=14, fontweight="bold", pad=4)
    ax_bottom_mid.imshow(panel_middle)
    ax_bottom_mid.set_title("Binding Site Mask", fontsize=14, fontweight="bold", pad=4)
    ax_bottom_right.imshow(panel_right)
    ax_bottom_right.set_title("Confusion Matrix Classification", fontsize=14, fontweight="bold", pad=4)

    for axis, label in (
        (ax_hsqc, "A."),
        (ax_csp, "B."),
        (ax_bottom_left, "C."),
        (ax_bottom_mid, "D."),
        (ax_bottom_right, "E."),
    ):
        _add_panel_letter(axis, label)

    for axis in (ax_hsqc, ax_bottom_left, ax_bottom_mid, ax_bottom_right):
        axis.axis("off")

    graphics.plot_alignment_similarity_based(
        ax_aln,
        alignment,
        matrix=alignment_matrix,
        labels=["Apo", "Holo"],
        show_numbers=True,
        show_line_position=True,
        symbols_per_line=_SYMBOLS_PER_LINE,
    )
    ax_aln.set_title(
        "Apo / Holo CS-list Seq_ID alignment",
        fontsize=14,
        fontweight="bold",
        pad=6,
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _default_panel_paths(target_dir: str, *, for_render: bool = False) -> Dict[str, str]:
    """
    When re-rendering, write case_study3_* to avoid overwriting v1 assets.
    Reuse looks at case_study3_* first, then v1 case_study_*.
    """
    assets = os.path.join(target_dir, "case_study_assets")
    if for_render:
        return {
            "left": os.path.join(assets, "case_study3_color_csp_mask.png"),
            "middle": os.path.join(assets, "case_study3_color_occlusion.png"),
            "right": os.path.join(assets, "case_study3_csp_classification_original.png"),
        }
    return {
        "left": os.path.join(assets, "case_study_color_csp_mask.png"),
        "middle": os.path.join(assets, "case_study_color_occlusion.png"),
        "right": os.path.join(assets, "case_study_csp_classification_original.png"),
    }


def _find_reusable_v1_panels(
    target_dir: str,
    panel_search_dirs: Optional[Sequence[str]] = None,
) -> Optional[Dict[str, str]]:
    """Return existing PyMOL panel paths (case_study3_* first, then v1 case_study_*)."""
    candidates: List[str] = [target_dir]
    if panel_search_dirs:
        candidates.extend(str(p) for p in panel_search_dirs)
    seen = set()
    for base in candidates:
        base = os.path.abspath(base)
        if base in seen:
            continue
        seen.add(base)
        for for_render in (True, False):
            paths = _default_panel_paths(base, for_render=for_render)
            if all(os.path.exists(p) for p in paths.values()):
                return paths
    return None


def generate_case_study_3_figure(
    target_dir: str,
    pdb_id: str,
    apo_bmrb: Optional[str] = None,
    holo_bmrb: Optional[str] = None,
    apo_pdb: Optional[str] = None,
    force_view_reset: bool = False,
    view_key: Optional[str] = None,
    reuse_existing_panels: bool = True,
    cs_dir: Optional[str] = None,
    do_fetch: bool = False,
    panel_search_dirs: Optional[Sequence[str]] = None,
    allow_interactive_view: bool = True,
) -> str:
    """
    Build `{pdb_id}_case_study_3.png` in target_dir.

    When reuse_existing_panels is True and case_study3_* or v1 case_study_*
    panel PNGs exist (in target_dir or panel_search_dirs), skip PyMOL.
    """
    if not apo_bmrb or not holo_bmrb:
        raise ValueError("apo_bmrb and holo_bmrb are required for case_study_3")

    hsqc_scatter_path = os.path.join(target_dir, "hsqc_scatter.png")
    csp_bars_path = os.path.join(target_dir, "csp_classification_bars_original.png")
    for path in (hsqc_scatter_path, csp_bars_path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing case-study-3 input: {path}")

    reuse_paths = (
        _find_reusable_v1_panels(target_dir, panel_search_dirs)
        if reuse_existing_panels
        else None
    )
    if reuse_paths is not None:
        rendered = reuse_paths
        print(
            f"[CASE_STUDY_3] Reusing existing PyMOL panels from "
            f"{os.path.dirname(reuse_paths['left'])}"
        )
    else:
        color_csp_mask_pml = os.path.join(target_dir, "color_csp_mask.pml")
        color_occlusion_pml = os.path.join(target_dir, "color_occlusion.pml")
        classification_pml = os.path.join(target_dir, "csp_classification_original.pml")
        required = [color_csp_mask_pml, color_occlusion_pml, classification_pml]
        missing = [p for p in required if not os.path.exists(p)]
        if missing:
            raise FileNotFoundError(
                "Missing case-study-3 PyMOL inputs (and no reusable v1 panels):\n"
                + "\n".join(missing)
            )

        pymol_views_dir = paths.pymol_views_dir
        os.makedirs(pymol_views_dir, exist_ok=True)
        view_id = (view_key or os.path.basename(os.path.abspath(target_dir)) or pdb_id).strip()
        view_path_save = _case_study_view_save_path(pymol_views_dir, pdb_id, view_key)
        view: Optional[List[float]] = None
        if force_view_reset:
            print(f"[CASE_STUDY_3] Forcing view recapture for {view_id}; ignoring saved view.")
        else:
            for cand in _case_study_view_load_candidates(pymol_views_dir, pdb_id, view_key):
                if not os.path.exists(cand):
                    continue
                try:
                    view = _load_view(cand)
                    print(f"[CASE_STUDY_3] Reusing saved view: {cand}")
                    break
                except Exception as exc:
                    print(f"[CASE_STUDY_3] WARNING: Saved view is invalid ({cand}): {exc}")

        if view is None:
            if not allow_interactive_view:
                raise RuntimeError(
                    f"No saved PyMOL view for {view_id} and interactive capture is disabled."
                )
            capture_user_view_interactive(
                color_csp_mask_pml_path=color_csp_mask_pml,
                view_output_path=view_path_save,
                pdb_id=pdb_id,
            )
            view = _load_view(view_path_save)

        assets_dir = os.path.join(target_dir, "case_study_assets")
        os.makedirs(assets_dir, exist_ok=True)
        rendered = _default_panel_paths(target_dir, for_render=True)
        render_pymol_panel_with_view(color_csp_mask_pml, view, rendered["left"])
        render_pymol_panel_with_view(color_occlusion_pml, view, rendered["middle"])
        render_pymol_panel_with_view(classification_pml, view, rendered["right"])

    cs_path = Path(cs_dir) if cs_dir else Path(Paths().cs_cache_dir)
    apo_star = ensure_star_path(str(apo_bmrb), cs_path, do_fetch)
    holo_star = ensure_star_path(str(holo_bmrb), cs_path, do_fetch)
    if apo_star is None or holo_star is None:
        raise FileNotFoundError(
            f"Missing NMR-STAR for apo={apo_bmrb} and/or holo={holo_bmrb} under {cs_path}"
        )

    apo_sequences = parse_sequence_and_shifts_from_saveframes(str(apo_star))
    holo_sequences = parse_sequence_and_shifts_from_saveframes(str(holo_star))
    if not apo_sequences or not holo_sequences:
        raise RuntimeError("no H+N CS-list saveframes for apo/holo STAR files")

    pair = best_seqid_pair_alignment(apo_sequences, holo_sequences)
    apo_ec = extract_conditions_from_star(str(apo_star))
    holo_ec = extract_conditions_from_star(str(holo_star))
    apo_entities = entity_names_joined(apo_star)
    holo_entities = entity_names_joined(holo_star)
    apo_components = extract_sample_components_from_star(str(apo_star))
    holo_components = extract_sample_components_from_star(str(holo_star))

    conditions_lines = build_conditions_text(
        apo_bmrb=str(apo_bmrb),
        holo_bmrb=str(holo_bmrb),
        apo_pdb=apo_pdb or "",
        holo_pdb=pdb_id,
        apo_ec=apo_ec,
        holo_ec=holo_ec,
        apo_entities=apo_entities,
        holo_entities=holo_entities,
        apo_components=apo_components,
        holo_components=holo_components,
    )

    hsqc_panel = extract_hsqc_bottom_right_panel(hsqc_scatter_path)
    out_path = os.path.join(target_dir, f"{pdb_id.upper()}_case_study_3.png")
    header = format_case_study_metadata_header(pdb_id, apo_bmrb, holo_bmrb, apo_pdb)
    compose_case_study_3_figure(
        hsqc_panel,
        csp_bars_path,
        rendered,
        out_path,
        conditions_lines=conditions_lines,
        alignment=pair.alignment,
        alignment_matrix=pair.matrix,
        case_study_header=header,
        hsqc_panel_title=format_hsqc_offset_panel_title(),
    )
    return out_path


def write_one_slide_pptx(png_path: str, pptx_path: str) -> str:
    prs = Presentation()
    prs.slide_width = Inches(13.333333)
    prs.slide_height = Inches(7.5)
    add_image_slide(prs, Path(png_path))
    Path(pptx_path).parent.mkdir(parents=True, exist_ok=True)
    prs.save(pptx_path)
    return pptx_path


def _batch_generate_from_csv(
    csv_path: Path,
    *,
    outputs_dir: Path,
    panel_roots: Sequence[Path],
    cs_dir: Optional[str],
    do_fetch: bool,
    write_pptx: bool,
) -> int:
    import csv as _csv

    from scripts.target_resolution import (
        build_resolution_caches,
        canonical_output_dir_name,
        resolve_output_dir_from_csv_row,
    )

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(_csv.DictReader(handle))

    caches = build_resolution_caches(outputs_dir)
    n_ok = 0
    n_fail = 0
    failures: List[str] = []

    for i, row in enumerate(rows, start=1):
        apo_bmrb = (row.get("apo_bmrb") or "").strip()
        holo_bmrb = (row.get("holo_bmrb") or "").strip()
        apo_pdb = (row.get("apo_pdb") or "").strip()
        holo_pdb = (row.get("holo_pdb") or "").strip()
        label = canonical_output_dir_name(holo_pdb, apo_bmrb)
        print(f"\n[CASE_STUDY_3] ({i}/{len(rows)}) {label}")

        target = resolve_output_dir_from_csv_row(row, caches=caches)
        if target is None:
            target = outputs_dir / label
        if not target.is_dir():
            n_fail += 1
            failures.append(f"{label}: missing outputs dir {target}")
            print(f"[CASE_STUDY_3] FAIL: missing outputs dir {target}")
            continue

        panel_search = []
        for root in panel_roots:
            cand = root / label
            if cand.is_dir():
                panel_search.append(str(cand))

        try:
            out_png = generate_case_study_3_figure(
                target_dir=str(target),
                pdb_id=holo_pdb,
                apo_bmrb=apo_bmrb,
                holo_bmrb=holo_bmrb,
                apo_pdb=apo_pdb or None,
                view_key=label,
                reuse_existing_panels=True,
                cs_dir=cs_dir,
                do_fetch=do_fetch,
                panel_search_dirs=panel_search,
                allow_interactive_view=False,
            )
            print(f"[CASE_STUDY_3] Wrote {out_png}")
            if write_pptx:
                pptx_path = write_one_slide_pptx(out_png, str(Path(out_png).with_suffix(".pptx")))
                print(f"[CASE_STUDY_3] Wrote {pptx_path}")
            n_ok += 1
        except Exception as exc:
            n_fail += 1
            failures.append(f"{label}: {exc}")
            print(f"[CASE_STUDY_3] FAIL: {exc}")

    print(f"\n[CASE_STUDY_3] Batch done: ok={n_ok} fail={n_fail} total={len(rows)}")
    if failures:
        print("[CASE_STUDY_3] Failures:")
        for line in failures:
            print(f"  {line}")
    return 0 if n_fail == 0 else 1


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate case_study_3 figure (conditions + case_study v1 panels + alignment)."
    )
    parser.add_argument(
        "--csv",
        default="",
        help="Batch mode: generate for every row in this CSP CSV (e.g. data/CSP_UBQ.csv)",
    )
    parser.add_argument(
        "--outputs-dir",
        default="",
        help="Batch mode outputs root (default: Paths().outputs_dir)",
    )
    parser.add_argument(
        "--panel-roots",
        default="outputs",
        help="Comma-separated roots to search for case_study_assets (batch mode)",
    )
    parser.add_argument(
        "--target-dir",
        default="",
        help="Single-target output directory containing HSQC / CSP bars / case_study assets",
    )
    parser.add_argument("--holo-pdb", default="")
    parser.add_argument("--apo-bmrb", default="")
    parser.add_argument("--holo-bmrb", default="")
    parser.add_argument("--apo-pdb", default="")
    parser.add_argument(
        "--cs-dir",
        default="",
        help="NMR-STAR cache directory (default: Paths().cs_cache_dir)",
    )
    parser.add_argument("--fetch-bmrb", action="store_true")
    parser.add_argument(
        "--force-view-reset",
        action="store_true",
        help="Recapture PyMOL view (ignored when reusing existing panels)",
    )
    parser.add_argument(
        "--no-reuse-panels",
        action="store_true",
        help="Re-render PyMOL panels even if case_study (v1) assets exist",
    )
    parser.add_argument(
        "--pptx",
        default="",
        help="Optional path for a one-slide PPTX review deck (single-target mode)",
    )
    parser.add_argument(
        "--no-pptx",
        action="store_true",
        help="Do not write companion PPTX (default in batch mode)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.csv:
        csv_path = Path(args.csv)
        if not csv_path.is_file():
            raise SystemExit(f"Missing CSV: {csv_path}")
        outputs_dir = Path(args.outputs_dir) if args.outputs_dir else Path(Paths().outputs_dir)
        panel_roots = [
            Path(p.strip()) for p in args.panel_roots.split(",") if p.strip()
        ]
        raise SystemExit(
            _batch_generate_from_csv(
                csv_path,
                outputs_dir=outputs_dir,
                panel_roots=panel_roots,
                cs_dir=args.cs_dir or None,
                do_fetch=args.fetch_bmrb,
                write_pptx=bool(args.pptx) and not args.no_pptx,
            )
        )

    if not args.target_dir or not args.holo_pdb or not args.apo_bmrb or not args.holo_bmrb:
        raise SystemExit(
            "Single-target mode requires --target-dir --holo-pdb --apo-bmrb --holo-bmrb "
            "(or use --csv for batch mode)."
        )

    out_png = generate_case_study_3_figure(
        target_dir=args.target_dir,
        pdb_id=args.holo_pdb,
        apo_bmrb=args.apo_bmrb,
        holo_bmrb=args.holo_bmrb,
        apo_pdb=args.apo_pdb or None,
        force_view_reset=args.force_view_reset,
        reuse_existing_panels=not args.no_reuse_panels,
        cs_dir=args.cs_dir or None,
        do_fetch=args.fetch_bmrb,
    )
    print(f"[CASE_STUDY_3] Wrote {out_png}")

    if args.no_pptx:
        return
    if args.pptx:
        pptx_path = write_one_slide_pptx(out_png, args.pptx)
        print(f"[CASE_STUDY_3] Wrote {pptx_path}")
    else:
        default_pptx = str(Path(out_png).with_suffix(".pptx"))
        pptx_path = write_one_slide_pptx(out_png, default_pptx)
        print(f"[CASE_STUDY_3] Wrote {pptx_path}")


if __name__ == "__main__":
    main()
