#!/usr/bin/env python3
"""
Case Study Z: case_study (v1) layout plus a bottom CSP z-score vs nearest
interchain atom–atom distance scatter.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

try:
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
    from .config import paths
    from .plot_csp_z_vs_ca_distance import (
        ANY_ATOM_DISTANCE_COLUMN,
        ANY_ATOM_DISTANCE_XLABEL,
        collect_points_from_master_csv,
        _draw_csp_z_vs_distance_scatter,
    )
except Exception:
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
    from scripts.config import paths  # type: ignore
    from scripts.plot_csp_z_vs_ca_distance import (  # type: ignore
        ANY_ATOM_DISTANCE_COLUMN,
        ANY_ATOM_DISTANCE_XLABEL,
        collect_points_from_master_csv,
        _draw_csp_z_vs_distance_scatter,
    )


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


_REUSE_PANEL_STEMS: Tuple[Tuple[str, str, str], ...] = (
    (
        "case_study_color_csp_mask.png",
        "case_study_color_occlusion.png",
        "case_study_csp_classification_original.png",
    ),
    (
        "case_study3_color_csp_mask.png",
        "case_study3_color_occlusion.png",
        "case_study3_csp_classification_original.png",
    ),
)


def _default_panel_paths(target_dir: str, *, for_render: bool = False) -> Dict[str, str]:
    assets = os.path.join(target_dir, "case_study_assets")
    if for_render:
        return {
            "left": os.path.join(assets, "case_studyz_color_csp_mask.png"),
            "middle": os.path.join(assets, "case_studyz_color_occlusion.png"),
            "right": os.path.join(assets, "case_studyz_csp_classification_original.png"),
        }
    left, middle, right = _REUSE_PANEL_STEMS[0]
    return {
        "left": os.path.join(assets, left),
        "middle": os.path.join(assets, middle),
        "right": os.path.join(assets, right),
    }


def _panel_paths_if_complete(assets_dir: str, stems: Tuple[str, str, str]) -> Optional[Dict[str, str]]:
    paths = {
        "left": os.path.join(assets_dir, stems[0]),
        "middle": os.path.join(assets_dir, stems[1]),
        "right": os.path.join(assets_dir, stems[2]),
    }
    if all(os.path.exists(p) for p in paths.values()):
        return paths
    return None


def _find_reusable_v1_panels(
    target_dir: str,
    panel_search_dirs: Optional[Sequence[str]] = None,
) -> Optional[Dict[str, str]]:
    candidates: List[str] = [target_dir]
    if panel_search_dirs:
        candidates.extend(str(p) for p in panel_search_dirs)
    seen = set()
    for base in candidates:
        base = os.path.abspath(base)
        if base in seen:
            continue
        seen.add(base)
        assets = os.path.join(base, "case_study_assets")
        for stems in _REUSE_PANEL_STEMS:
            found = _panel_paths_if_complete(assets, stems)
            if found is not None:
                return found
    return None


_REQUIRED_INPUTS = (
    "hsqc_scatter.png",
    "csp_classification_bars_original.png",
    "master_alignment.csv",
)


def _ids_from_csp_table(target_dir: Path) -> Dict[str, str]:
    table = target_dir / "csp_table.csv"
    if not table.is_file():
        return {}
    with table.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle), None)
    if not row:
        return {}
    return {
        "apo_bmrb": (row.get("apo_bmrb") or "").strip(),
        "holo_bmrb": (row.get("holo_bmrb") or "").strip(),
        "holo_pdb": (row.get("holo_pdb") or "").strip(),
    }


def _load_apo_pdb_lookup(csv_path: Path) -> Dict[Tuple[str, str], str]:
    lookup: Dict[Tuple[str, str], str] = {}
    if not csv_path.is_file():
        return lookup
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            holo = (row.get("holo_pdb") or "").strip().upper()
            apo_bmrb = (row.get("apo_bmrb") or "").strip()
            apo_pdb = (row.get("apo_pdb") or "").strip()
            if holo and apo_bmrb and apo_pdb:
                lookup[(holo, apo_bmrb)] = apo_pdb
    return lookup


def _discover_output_targets(outputs_dir: Path) -> List[Path]:
    targets: List[Path] = []
    for path in sorted(outputs_dir.iterdir(), key=lambda p: p.name.lower()):
        if not path.is_dir() or path.name.startswith("."):
            continue
        if all((path / name).is_file() for name in _REQUIRED_INPUTS):
            targets.append(path)
    return targets


def _panel_search_dirs(label: str, panel_roots: Sequence[Path]) -> List[str]:
    found: List[str] = []
    for root in panel_roots:
        cand = root / label
        if cand.is_dir():
            found.append(str(cand))
        elif root.is_dir():
            found.append(str(root))
    return found


def _batch_generate_from_outputs(
    outputs_dir: Path,
    *,
    panel_roots: Sequence[Path],
    csv_path: Optional[Path],
    reuse_existing_panels: bool,
    allow_interactive_view: bool,
    force_view_reset: bool,
) -> int:
    apo_pdb_lookup = _load_apo_pdb_lookup(csv_path) if csv_path is not None else {}
    targets = _discover_output_targets(outputs_dir)
    if not targets:
        print(f"[CASE_STUDY_Z] No target directories with required inputs under {outputs_dir}")
        return 1

    n_ok = 0
    n_fail = 0
    failures: List[str] = []
    for i, target in enumerate(targets, start=1):
        label = target.name
        ids = _ids_from_csp_table(target)
        if "_" in label:
            dirname_pdb, dirname_apo = label.rsplit("_", 1)
        else:
            dirname_pdb, dirname_apo = label, ""
        holo_pdb = ids.get("holo_pdb") or dirname_pdb
        apo_bmrb = ids.get("apo_bmrb") or dirname_apo
        holo_bmrb = ids.get("holo_bmrb") or ""
        apo_pdb = apo_pdb_lookup.get((holo_pdb.upper(), apo_bmrb), "")
        print(f"\n[CASE_STUDY_Z] ({i}/{len(targets)}) {label}")
        try:
            out_png = generate_case_study_z_figure(
                target_dir=str(target),
                pdb_id=holo_pdb,
                apo_bmrb=apo_bmrb or None,
                holo_bmrb=holo_bmrb or None,
                apo_pdb=apo_pdb or None,
                force_view_reset=force_view_reset,
                view_key=label,
                reuse_existing_panels=reuse_existing_panels,
                panel_search_dirs=_panel_search_dirs(label, panel_roots) or None,
                allow_interactive_view=allow_interactive_view,
            )
            print(f"[CASE_STUDY_Z] Wrote {out_png}")
            n_ok += 1
        except Exception as exc:
            n_fail += 1
            failures.append(f"{label}: {exc}")
            print(f"[CASE_STUDY_Z] FAIL: {exc}")

    print(f"\n[CASE_STUDY_Z] Batch done: ok={n_ok} fail={n_fail} total={len(targets)}")
    if failures:
        print("[CASE_STUDY_Z] Failures:")
        for line in failures:
            print(f"  {line}")
    return 0 if n_fail == 0 else 1


def compose_case_study_z_figure(
    hsqc_bottom_right_panel,
    csp_classification_plot_path: str,
    rendered_panels: Dict[str, str],
    scatter_points,
    output_path: str,
    *,
    case_study_header: Optional[str] = None,
    hsqc_panel_title: Optional[str] = None,
    scatter_title: str = "CSP Z-Score vs Nearest Interchain Atom–Atom Distance",
) -> None:
    """
    Build case-study-Z figure:
      - top: HSQC | CSP classification bars
      - middle: three PyMOL renders
      - bottom: CSP z vs nearest interchain atom–atom distance scatter
    """
    csp_plot = _trim_whitespace(mpimg.imread(csp_classification_plot_path))
    panel_left = _trim_whitespace(mpimg.imread(rendered_panels["left"]))
    panel_middle = _trim_whitespace(mpimg.imread(rendered_panels["middle"]))
    panel_right = _trim_whitespace(mpimg.imread(rendered_panels["right"]))
    hsqc_panel = _trim_whitespace(hsqc_bottom_right_panel)

    fig = plt.figure(figsize=(22, 16))
    outer = GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[2.15, 1.0],
        hspace=0.12,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.05)
    if case_study_header:
        fig.suptitle(case_study_header, fontsize=18, fontweight="bold", x=0.5, y=0.985)

    main = outer[0].subgridspec(
        2,
        1,
        height_ratios=[1.0, 1.15],
        hspace=0.08,
    )
    top = main[0].subgridspec(1, 2, wspace=0.02)
    mid = main[1].subgridspec(1, 3, wspace=0.02)

    ax_hsqc = fig.add_subplot(top[0, 0])
    ax_csp = fig.add_subplot(top[0, 1])
    ax_bottom_left = fig.add_subplot(mid[0, 0])
    ax_bottom_mid = fig.add_subplot(mid[0, 1])
    ax_bottom_right = fig.add_subplot(mid[0, 2])
    ax_scatter = fig.add_subplot(outer[1])

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

    _draw_csp_z_vs_distance_scatter(
        ax_scatter,
        scatter_points,
        title=scatter_title,
        xlabel=ANY_ATOM_DISTANCE_XLABEL,
        legend_fontsize=10,
        label_fontsize=12,
        title_fontsize=14,
        marker_size=36,
        fp_text_fontsize=12,
    )
    _add_panel_letter(ax_scatter, "F.")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_case_study_z_figure(
    target_dir: str,
    pdb_id: str,
    apo_bmrb: Optional[str] = None,
    holo_bmrb: Optional[str] = None,
    apo_pdb: Optional[str] = None,
    force_view_reset: bool = False,
    view_key: Optional[str] = None,
    reuse_existing_panels: bool = True,
    panel_search_dirs: Optional[Sequence[str]] = None,
    allow_interactive_view: bool = True,
) -> str:
    """Build `{pdb_id}_case_study_z.png` in target_dir."""
    hsqc_scatter_path = os.path.join(target_dir, "hsqc_scatter.png")
    csp_bars_path = os.path.join(target_dir, "csp_classification_bars_original.png")
    master_csv = os.path.join(target_dir, "master_alignment.csv")
    for path in (hsqc_scatter_path, csp_bars_path, master_csv):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing case-study-Z input: {path}")

    reuse_paths = (
        _find_reusable_v1_panels(target_dir, panel_search_dirs)
        if reuse_existing_panels
        else None
    )
    if reuse_paths is not None:
        rendered = reuse_paths
        print(
            f"[CASE_STUDY_Z] Reusing existing case_study (v1) panels from "
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
                "Missing case-study-Z PyMOL inputs (and no reusable v1 panels):\n"
                + "\n".join(missing)
            )

        pymol_views_dir = paths.pymol_views_dir
        os.makedirs(pymol_views_dir, exist_ok=True)
        view_id = (view_key or os.path.basename(os.path.abspath(target_dir)) or pdb_id).strip()
        view_path_save = _case_study_view_save_path(pymol_views_dir, pdb_id, view_key)
        view: Optional[List[float]] = None
        if force_view_reset:
            print(f"[CASE_STUDY_Z] Forcing view recapture for {view_id}; ignoring saved view.")
        else:
            for cand in _case_study_view_load_candidates(pymol_views_dir, pdb_id, view_key):
                if not os.path.exists(cand):
                    continue
                try:
                    view = _load_view(cand)
                    print(f"[CASE_STUDY_Z] Reusing saved view: {cand}")
                    break
                except Exception as exc:
                    print(f"[CASE_STUDY_Z] WARNING: Saved view is invalid ({cand}): {exc}")

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

    points = collect_points_from_master_csv(
        master_csv,
        distance_column=ANY_ATOM_DISTANCE_COLUMN,
    )
    if not points:
        raise RuntimeError(f"No CSP z / distance points found in {master_csv}")

    hsqc_panel = extract_hsqc_bottom_right_panel(hsqc_scatter_path)
    out_path = os.path.join(target_dir, f"{pdb_id.upper()}_case_study_z.png")
    header = format_case_study_metadata_header(pdb_id, apo_bmrb, holo_bmrb, apo_pdb)
    compose_case_study_z_figure(
        hsqc_panel,
        csp_bars_path,
        rendered,
        points,
        out_path,
        case_study_header=header,
        hsqc_panel_title=format_hsqc_offset_panel_title(),
    )
    return out_path


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate Case Study Z (v1 panels + CSP z vs atom–atom distance)."
    )
    parser.add_argument(
        "--outputs-dir",
        default="",
        help="Batch mode: generate for every target directory under this outputs root.",
    )
    parser.add_argument(
        "--csv",
        default="",
        help="Optional CSP CSV used in batch mode to fill apo_pdb in figure headers.",
    )
    parser.add_argument(
        "--target-dir",
        default="",
        help="Single-target output directory with HSQC / CSP bars / master_alignment.csv",
    )
    parser.add_argument("--holo-pdb", default="")
    parser.add_argument("--apo-bmrb", default="")
    parser.add_argument("--holo-bmrb", default="")
    parser.add_argument("--apo-pdb", default="")
    parser.add_argument(
        "--panel-roots",
        default="outputs",
        help="Comma-separated roots to search for case_study_assets",
    )
    parser.add_argument("--force-view-reset", action="store_true")
    parser.add_argument("--no-reuse-panels", action="store_true")
    parser.add_argument(
        "--no-interactive-view",
        action="store_true",
        help="Fail instead of opening PyMOL if no saved view / panels",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    panel_roots = [Path(p.strip()) for p in args.panel_roots.split(",") if p.strip()]

    if args.outputs_dir or (not args.target_dir):
        outputs_dir = Path(args.outputs_dir) if args.outputs_dir else Path(paths.outputs_dir)
        csv_path = Path(args.csv) if args.csv else _REPO / "data" / "CSP_UBQ.csv"
        raise SystemExit(
            _batch_generate_from_outputs(
                outputs_dir,
                panel_roots=panel_roots,
                csv_path=csv_path if csv_path.is_file() else None,
                reuse_existing_panels=not args.no_reuse_panels,
                allow_interactive_view=not args.no_interactive_view,
                force_view_reset=args.force_view_reset,
            )
        )

    if not args.holo_pdb:
        raise SystemExit("Single-target mode requires --target-dir and --holo-pdb.")

    target_dir = Path(args.target_dir)
    label = target_dir.name
    out_png = generate_case_study_z_figure(
        target_dir=str(target_dir),
        pdb_id=args.holo_pdb,
        apo_bmrb=args.apo_bmrb or None,
        holo_bmrb=args.holo_bmrb or None,
        apo_pdb=args.apo_pdb or None,
        force_view_reset=args.force_view_reset,
        view_key=label,
        reuse_existing_panels=not args.no_reuse_panels,
        panel_search_dirs=_panel_search_dirs(label, panel_roots) or None,
        allow_interactive_view=not args.no_interactive_view,
    )
    print(f"[CASE_STUDY_Z] Wrote {out_png}")


if __name__ == "__main__":
    main()
