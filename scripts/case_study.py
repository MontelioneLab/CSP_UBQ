"""
Utilities for assembling per-target case-study figures.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np


def extract_hsqc_bottom_right_panel(hsqc_png_path: str):
    """
    Extract the bottom-right panel from the 2x2 HSQC comparison image.
    """
    image = mpimg.imread(hsqc_png_path)
    height, width = image.shape[0], image.shape[1]
    # Start slightly below the midpoint so the crop excludes labels/titles
    # from the upper row while preserving the lower-right subplot.
    row_start = int(height * 0.53)
    col_start = width // 2
    return image[row_start:height, col_start:width]


def capture_user_view_interactive(
    color_csp_mask_pml_path: str,
    view_output_path: str,
    pdb_id: str,
) -> None:
    """
    Open PyMOL GUI, let user orient the model, and save the camera with F5.
    """
    helper_script = f"""
from pymol import cmd
import json

VIEW_PATH = r\"\"\"{view_output_path}\"\"\"

def _save_case_study_view():
    view = list(cmd.get_view())
    with open(VIEW_PATH, "w", encoding="utf-8") as handle:
        json.dump(view, handle)
    print(f"[CASE_STUDY] Saved view to {{VIEW_PATH}}")
    cmd.quit()

cmd.set_key("F5", _save_case_study_view)
print("[CASE_STUDY] ------------------------------------------------------------")
print("[CASE_STUDY] Case-study view capture for {pdb_id}")
print("[CASE_STUDY] 1) Set your desired model orientation.")
print("[CASE_STUDY] 2) Press F5 to save the view and close PyMOL.")
print("[CASE_STUDY] ------------------------------------------------------------")
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix="_case_study_capture.py", delete=False) as tmp:
        tmp.write(helper_script)
        helper_path = tmp.name

    try:
        subprocess.run(
            ["pymol", color_csp_mask_pml_path, "-r", helper_path],
            check=True,
        )
    finally:
        if os.path.exists(helper_path):
            os.remove(helper_path)

    if not os.path.exists(view_output_path):
        raise RuntimeError("No view was captured. Re-run and press F5 in the PyMOL window.")


def _load_view(view_output_path: str) -> List[float]:
    with open(view_output_path, "r", encoding="utf-8") as handle:
        view = json.load(handle)
    if not isinstance(view, list) or len(view) != 18:
        raise ValueError(f"Invalid PyMOL view data in {view_output_path}; expected list of 18 floats.")
    return [float(x) for x in view]


def _view_command(view: List[float]) -> str:
    return "set_view (" + ", ".join(f"{value:.9f}" for value in view) + ")"


def render_pymol_panel_with_view(
    pml_path: str,
    view: List[float],
    output_png_path: str,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 300,
) -> None:
    """
    Render a PyMOL script using a fixed camera view.
    """
    subprocess.run(
        [
            "pymol",
            "-c",
            "-q",
            pml_path,
            "-d",
            _view_command(view),
            "-d",
            f"ray {width}, {height}",
            "-d",
            f"png {output_png_path}, dpi={dpi}",
            "-d",
            "quit",
        ],
        check=True,
    )

    if not os.path.exists(output_png_path):
        raise RuntimeError(f"Failed to render PyMOL panel: {output_png_path}")


def compose_case_study_figure(
    hsqc_bottom_right_panel,
    csp_classification_plot_path: str,
    rendered_panels: Dict[str, str],
    output_path: str,
    case_study_header: Optional[str] = None,
) -> None:
    """
    Build case-study figure:
      - top-left: HSQC bottom-right panel
      - top-right: CSP classification bars
      - bottom row: three PyMOL renders
    """
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

    csp_plot = _trim_whitespace(mpimg.imread(csp_classification_plot_path))
    panel_left = _trim_whitespace(mpimg.imread(rendered_panels["left"]))
    panel_middle = _trim_whitespace(mpimg.imread(rendered_panels["middle"]))
    panel_right = _trim_whitespace(mpimg.imread(rendered_panels["right"]))
    hsqc_panel = _trim_whitespace(hsqc_bottom_right_panel)

    fig = plt.figure(figsize=(22, 12))
    grid = GridSpec(2, 6, figure=fig, height_ratios=[1.0, 1.15], wspace=0.02, hspace=0.06)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)
    if case_study_header:
        fig.suptitle(case_study_header, fontsize=18, fontweight="bold", x=0.5, y=0.99)

    ax_top_left = fig.add_subplot(grid[0, 0:3])
    ax_top_right = fig.add_subplot(grid[0, 3:6])
    ax_bottom_left = fig.add_subplot(grid[1, 0:2])
    ax_bottom_mid = fig.add_subplot(grid[1, 2:4])
    ax_bottom_right = fig.add_subplot(grid[1, 4:6])

    # Preserve native image ratios while maximizing occupied subplot area.
    ax_top_left.set_box_aspect(hsqc_panel.shape[0] / hsqc_panel.shape[1])
    ax_top_right.set_box_aspect(csp_plot.shape[0] / csp_plot.shape[1])
    ax_bottom_left.set_box_aspect(panel_left.shape[0] / panel_left.shape[1])
    ax_bottom_mid.set_box_aspect(panel_middle.shape[0] / panel_middle.shape[1])
    ax_bottom_right.set_box_aspect(panel_right.shape[0] / panel_right.shape[1])

    ax_top_left.imshow(hsqc_panel)
    ax_top_left.set_title("Apo/Holo HSQC Offset", fontsize=14, fontweight="bold", pad=4)
    ax_top_right.imshow(csp_plot)
    ax_top_right.set_title("CSP Classification", fontsize=14, fontweight="bold", pad=4)

    ax_bottom_left.imshow(panel_left)
    ax_bottom_left.set_title("CSP Mask", fontsize=14, fontweight="bold", pad=4)
    ax_bottom_mid.imshow(panel_middle)
    ax_bottom_mid.set_title("Binding Site Mask", fontsize=14, fontweight="bold", pad=4)
    ax_bottom_right.imshow(panel_right)
    ax_bottom_right.set_title("Confusion Matrix Classification", fontsize=14, fontweight="bold", pad=4)

    labels = (
        (ax_top_left, "A."),
        (ax_top_right, "B."),
        (ax_bottom_left, "C."),
        (ax_bottom_mid, "D."),
        (ax_bottom_right, "E."),
    )
    for axis, label in labels:
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

    for axis in (ax_top_left, ax_top_right, ax_bottom_left, ax_bottom_mid, ax_bottom_right):
        axis.axis("off")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_case_study_figure(
    target_dir: str,
    pdb_id: str,
    apo_bmrb: Optional[str] = None,
    holo_bmrb: Optional[str] = None,
    force_view_reset: bool = False,
) -> str:
    """
    Generate <pdb_id>_case_study.png in a target output directory.
    """
    hsqc_scatter_path = os.path.join(target_dir, "hsqc_scatter.png")
    csp_bars_path = os.path.join(target_dir, "csp_classification_bars_original.png")
    color_csp_mask_pml = os.path.join(target_dir, "color_csp_mask.pml")
    color_occlusion_pml = os.path.join(target_dir, "color_occlusion.pml")
    classification_pml = os.path.join(target_dir, "csp_classification_original.pml")

    required_files = [
        hsqc_scatter_path,
        csp_bars_path,
        color_csp_mask_pml,
        color_occlusion_pml,
        classification_pml,
    ]
    missing = [path for path in required_files if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError("Missing case-study input files:\n" + "\n".join(missing))

    try:
        from .config import paths
    except Exception:
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from scripts.config import paths  # type: ignore

    pymol_views_dir = paths.pymol_views_dir
    os.makedirs(pymol_views_dir, exist_ok=True)
    view_path = os.path.join(pymol_views_dir, f"{pdb_id}_case_study_view.json")

    assets_dir = os.path.join(target_dir, "case_study_assets")
    os.makedirs(assets_dir, exist_ok=True)
    view: Optional[List[float]] = None
    if force_view_reset:
        print(f"[CASE_STUDY] Forcing view recapture for {pdb_id}; ignoring saved view.")
    elif os.path.exists(view_path):
        try:
            view = _load_view(view_path)
            print(f"[CASE_STUDY] Reusing saved view: {view_path}")
        except Exception as exc:
            print(f"[CASE_STUDY] WARNING: Saved view is invalid ({exc}); recapturing.")

    if view is None:
        capture_user_view_interactive(
            color_csp_mask_pml_path=color_csp_mask_pml,
            view_output_path=view_path,
            pdb_id=pdb_id,
        )
        view = _load_view(view_path)

    rendered = {
        "left": os.path.join(assets_dir, "case_study_color_csp_mask.png"),
        "middle": os.path.join(assets_dir, "case_study_color_occlusion.png"),
        "right": os.path.join(assets_dir, "case_study_csp_classification_original.png"),
    }

    render_pymol_panel_with_view(color_csp_mask_pml, view, rendered["left"])
    render_pymol_panel_with_view(color_occlusion_pml, view, rendered["middle"])
    render_pymol_panel_with_view(classification_pml, view, rendered["right"])

    hsqc_panel = extract_hsqc_bottom_right_panel(hsqc_scatter_path)
    out_path = os.path.join(target_dir, f"{pdb_id}_case_study.png")
    header = f"holo_pdb: {pdb_id} | apo_bmrb: {apo_bmrb or 'N/A'} | holo_bmrb: {holo_bmrb or 'N/A'}"
    compose_case_study_figure(hsqc_panel, csp_bars_path, rendered, out_path, case_study_header=header)
    return out_path
