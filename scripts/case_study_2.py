"""
Alternative case-study layout (single-row, four panels).
"""

from __future__ import annotations

import json
import os
import sys
from typing import Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

try:
    from .case_study import capture_user_view_interactive, format_case_study_metadata_header, render_pymol_panel_with_view
except Exception:
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    from scripts.case_study import (  # type: ignore
        capture_user_view_interactive,
        format_case_study_metadata_header,
        render_pymol_panel_with_view,
    )


def _load_view(view_output_path: str) -> List[float]:
    with open(view_output_path, "r", encoding="utf-8") as handle:
        view = json.load(handle)
    if not isinstance(view, list) or len(view) != 18:
        raise ValueError(f"Invalid PyMOL view data in {view_output_path}; expected list of 18 floats.")
    return [float(x) for x in view]


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


def compose_case_study_2_figure(
    csp_classification_plot_path: str,
    rendered_panels: Dict[str, str],
    output_path: str,
    case_study_header: Optional[str] = None,
    panel_labels: Optional[Sequence[str]] = None,
) -> None:
    csp_plot = _trim_whitespace(mpimg.imread(csp_classification_plot_path))
    panel_mask = _trim_whitespace(mpimg.imread(rendered_panels["mask"]))
    panel_binding = _trim_whitespace(mpimg.imread(rendered_panels["binding"]))
    panel_confusion = _trim_whitespace(mpimg.imread(rendered_panels["confusion"]))

    # Give panel A more room while avoiding an overly wide final canvas.
    fig = plt.figure(figsize=(33, 7))
    grid = GridSpec(1, 4, figure=fig, wspace=0.005, width_ratios=[4.6, 1.8, 1.8, 1.8])
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.02)
    if case_study_header:
        fig.suptitle(case_study_header, fontsize=26, fontweight="bold", x=0.5, y=0.985)

    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[0, 3])

    ax_a.set_box_aspect(csp_plot.shape[0] / csp_plot.shape[1])
    ax_b.set_box_aspect(panel_mask.shape[0] / panel_mask.shape[1])
    ax_c.set_box_aspect(panel_binding.shape[0] / panel_binding.shape[1])
    ax_d.set_box_aspect(panel_confusion.shape[0] / panel_confusion.shape[1])

    ax_a.imshow(csp_plot)
    ax_b.imshow(panel_mask)
    ax_c.imshow(panel_binding)
    ax_d.imshow(panel_confusion)

    # Panel A: no title. Panels B, C, D: larger titles.
    ax_b.set_title("CSP Mask", fontsize=28, fontweight="bold", pad=4)
    ax_c.set_title("Binding Site Mask", fontsize=28, fontweight="bold", pad=4)
    ax_d.set_title("Confusion Matrix Classification", fontsize=28, fontweight="bold", pad=4)

    for axis in (ax_a, ax_b, ax_c, ax_d):
        axis.axis("off")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_case_study_2_figure(
    target_dir: str,
    pdb_id: str,
    apo_bmrb: Optional[str] = None,
    holo_bmrb: Optional[str] = None,
    apo_pdb: Optional[str] = None,
    force_view_reset: bool = False,
    panel_labels: Optional[Sequence[str]] = None,
    view_key: Optional[str] = None,
) -> str:
    csp_bars_path = os.path.join(target_dir, "csp_classification_bars_original.png")
    color_csp_mask_pml = os.path.join(target_dir, "color_csp_mask.pml")
    color_occlusion_pml = os.path.join(target_dir, "color_occlusion.pml")
    classification_pml = os.path.join(target_dir, "csp_classification_original.pml")

    required_files = [
        csp_bars_path,
        color_csp_mask_pml,
        color_occlusion_pml,
        classification_pml,
    ]
    missing = [path for path in required_files if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError("Missing case-study-2 input files:\n" + "\n".join(missing))

    try:
        from .config import paths
    except Exception:
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from scripts.config import paths  # type: ignore

    pymol_views_dir = paths.pymol_views_dir
    os.makedirs(pymol_views_dir, exist_ok=True)
    view_id = (view_key or os.path.basename(os.path.abspath(target_dir)) or pdb_id).strip()
    view_path = os.path.join(pymol_views_dir, f"{view_id}_case_study_view.json")

    assets_dir = os.path.join(target_dir, "case_study_assets")
    os.makedirs(assets_dir, exist_ok=True)
    view: Optional[List[float]] = None
    if force_view_reset:
        print(f"[CASE_STUDY_2] Forcing view recapture for {view_id}; ignoring saved view.")
    elif os.path.exists(view_path):
        try:
            view = _load_view(view_path)
            print(f"[CASE_STUDY_2] Reusing saved view: {view_path}")
        except Exception as exc:
            print(f"[CASE_STUDY_2] WARNING: Saved view is invalid ({exc}); recapturing.")

    if view is None:
        capture_user_view_interactive(
            color_csp_mask_pml_path=color_csp_mask_pml,
            view_output_path=view_path,
            pdb_id=pdb_id,
        )
        view = _load_view(view_path)

    rendered = {
        "mask": os.path.join(assets_dir, "case_study2_color_csp_mask.png"),
        "binding": os.path.join(assets_dir, "case_study2_color_occlusion.png"),
        "confusion": os.path.join(assets_dir, "case_study2_csp_classification_original.png"),
    }
    render_pymol_panel_with_view(color_csp_mask_pml, view, rendered["mask"])
    render_pymol_panel_with_view(color_occlusion_pml, view, rendered["binding"])
    render_pymol_panel_with_view(classification_pml, view, rendered["confusion"])

    out_path = os.path.join(target_dir, f"{pdb_id}_case_study_2.png")
    header = format_case_study_metadata_header(pdb_id, apo_bmrb, holo_bmrb, apo_pdb)
    compose_case_study_2_figure(
        csp_bars_path,
        rendered,
        out_path,
        case_study_header=header,
        panel_labels=panel_labels,
    )
    return out_path

