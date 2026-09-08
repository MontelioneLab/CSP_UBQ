#!/usr/bin/env python3
"""
Create publication Figure 2: four stacked case-study-2 rows.

Each row is one target (A.–D.). Columns are CSP Classification, CSP Mask,
Binding Site, and TP/FP/TN/FN. PyMOL panels are sublabeled i, ii, iii.

Individual bar plots and PyMOL renders are written under
``figures/figure_2_assets/<HOLO>/`` and then composed into
``figures/figure_2.png``.

Example::

    python scripts/create_fig_2.py --outputs-dir outputs
    python scripts/create_fig_2.py --outputs-dir outputs --skip-pymol
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.ndimage import rotate as nd_rotate

_REPO_ROOT = Path(__file__).resolve().parent.parent

try:
    from .case_study import (
        _case_study_view_load_candidates,
        render_pymol_panel_with_view,
    )
    from .case_study_2 import _trim_whitespace
    from .config import paths
    from .plot_csp_classification_bars_simple_legend import (
        DEFAULT_YLABEL,
        plot_csp_classification_bars_simple_legend,
    )
    from .target_resolution import TargetRow, build_resolution_caches, resolve_row
except Exception:
    _sys_path = str(_REPO_ROOT)
    if _sys_path not in sys.path:
        sys.path.insert(0, _sys_path)
    from scripts.case_study import (  # type: ignore
        _case_study_view_load_candidates,
        render_pymol_panel_with_view,
    )
    from scripts.case_study_2 import _trim_whitespace  # type: ignore
    from scripts.config import paths  # type: ignore
    from scripts.plot_csp_classification_bars_simple_legend import (  # type: ignore
        DEFAULT_YLABEL,
        plot_csp_classification_bars_simple_legend,
    )
    from scripts.target_resolution import (  # type: ignore
        TargetRow,
        build_resolution_caches,
        resolve_row,
    )


COLUMN_TITLES: Sequence[str] = (
    "CSP Classification",
    "CSP Mask",
    "Binding Site",
    "TP/FP/TN/FN",
)
PYMOL_SUBLABELS: Sequence[str] = ("i", "ii", "iii")
XLABEL_COLUMN = "Receptor Sequence"
WIDTH_RATIOS: Sequence[float] = (4.6, 1.8, 1.8, 1.8)
ROTATION_STEP_DEG = 2.0
CONTENT_WHITE = 0.985

PML_MASK = "color_csp_mask.pml"
PML_BINDING = "color_occlusion.pml"
PML_CONFUSION = "csp_classification_original.pml"


@dataclass(frozen=True)
class Fig2Target:
    panel: str
    holo_pdb: str
    apo_bmrb: str
    holo_bmrb: str

    @property
    def row(self) -> TargetRow:
        return TargetRow(
            holo_pdb=self.holo_pdb,
            apo_bmrb=self.apo_bmrb,
            holo_bmrb=self.holo_bmrb,
        )


FIG2_TARGETS: Sequence[Fig2Target] = (
    Fig2Target("A.", "7JQ8", "30782", "30786"),
    Fig2Target("B.", "2M14", "6225", "18842"),
    Fig2Target("C.", "6FDT", "19757", "34224"),
    Fig2Target("D.", "2KWV", "17769", "16885"),
)


@dataclass
class TargetAssets:
    target: Fig2Target
    target_dir: Path
    bars: Path
    mask: Path
    binding: Path
    confusion: Path


def _as_float_image(image) -> np.ndarray:
    arr = np.asarray(image)
    if arr.dtype.kind in "ui":
        return arr.astype(np.float32) / 255.0
    arr = arr.astype(np.float32, copy=False)
    if arr.size and float(np.nanmax(arr)) > 1.5:
        return arr / 255.0
    return arr


def _content_mask(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 2:
        return arr < CONTENT_WHITE
    rgb = arr[..., :3]
    if arr.shape[-1] == 4:
        alpha = arr[..., 3]
        return (alpha > 0.01) & np.any(rgb < CONTENT_WHITE, axis=-1)
    return np.any(rgb < CONTENT_WHITE, axis=-1)


def _bbox_size(mask: np.ndarray) -> Tuple[int, int]:
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return int(mask.shape[0]), int(mask.shape[1])
    return int(rows[-1] - rows[0] + 1), int(cols[-1] - cols[0] + 1)


def _rotate_image(image: np.ndarray, angle: float) -> np.ndarray:
    if abs(angle) < 1e-6:
        return image
    return nd_rotate(image, angle, reshape=True, order=1, cval=1.0, prefilter=False)


def _pymol_cell_size_inches(
    *,
    fig_width: float,
    fig_height: float,
    n_rows: int,
    left: float,
    right: float,
    top: float,
    bottom: float,
) -> Tuple[float, float]:
    usable_w = fig_width * (right - left)
    usable_h = fig_height * (top - bottom)
    ratio_sum = float(sum(WIDTH_RATIOS))
    cell_w = usable_w * (float(WIDTH_RATIOS[1]) / ratio_sum)
    cell_h = usable_h / float(n_rows)
    return cell_w, cell_h


def optimal_pymol_rotation_angle(
    image: np.ndarray,
    cell_w: float,
    cell_h: float,
    *,
    step_deg: float = ROTATION_STEP_DEG,
) -> float:
    """Return the in-plane angle that maximizes scale of the content bbox in the cell."""
    mask = _content_mask(image)
    if not mask.any():
        return 0.0
    # Downsample the mask so a 2° sweep stays cheap.
    stride = max(1, int(min(mask.shape[:2]) / 200))
    small = mask[::stride, ::stride]
    best_angle = 0.0
    best_score = -1.0
    for angle in np.arange(0.0, 180.0, step_deg):
        rot = nd_rotate(small.astype(np.float32), float(angle), reshape=True, order=0, cval=0.0)
        h, w = _bbox_size(rot > 0.5)
        if w <= 0 or h <= 0:
            continue
        score = min(cell_w / float(w), cell_h / float(h))
        if score > best_score:
            best_score = score
            best_angle = float(angle)
    return best_angle


def rotate_and_trim(image: np.ndarray, angle: float) -> np.ndarray:
    return _trim_whitespace(_rotate_image(image, angle))


def _load_view(view_output_path: Path) -> List[float]:
    with view_output_path.open("r", encoding="utf-8") as handle:
        view = json.load(handle)
    if not isinstance(view, list) or len(view) != 18:
        raise ValueError(
            f"Invalid PyMOL view data in {view_output_path}; expected list of 18 floats."
        )
    return [float(x) for x in view]


def _find_saved_view(pdb_id: str, view_key: str) -> Path:
    candidates = _case_study_view_load_candidates(paths.pymol_views_dir, pdb_id, view_key)
    for cand in candidates:
        path = Path(cand)
        if path.is_file():
            return path
    tried = "\n".join(candidates)
    raise FileNotFoundError(
        f"No saved PyMOL view for {pdb_id} (view_key={view_key!r}). Tried:\n{tried}"
    )


def _resolve_target_dir(
    target: Fig2Target,
    outputs_index: Dict[str, Path],
    bmrb_cache,
) -> Path:
    path = resolve_row(target.row, outputs_index, bmrb_cache)
    if path is None:
        raise FileNotFoundError(
            "No outputs subdirectory for "
            f"holo_pdb={target.holo_pdb} apo_bmrb={target.apo_bmrb} "
            f"holo_bmrb={target.holo_bmrb}. "
            "Pass --outputs-dir pointing at the pipeline tree "
            "(e.g. outputs)."
        )
    return path


def _asset_paths(assets_dir: Path, target: Fig2Target) -> Tuple[Path, Path, Path, Path]:
    row_dir = assets_dir / target.holo_pdb.upper()
    return (
        row_dir / "csp_classification.png",
        row_dir / "csp_mask.png",
        row_dir / "binding_site.png",
        row_dir / "confusion.png",
    )


def _required_pml(target_dir: Path) -> Dict[str, Path]:
    files = {
        "mask": target_dir / PML_MASK,
        "binding": target_dir / PML_BINDING,
        "confusion": target_dir / PML_CONFUSION,
    }
    missing = [str(p) for p in files.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing PyMOL scripts for Figure 2:\n" + "\n".join(missing)
        )
    return files


# Top y-tick labels that collide with row letters on the composite.
_HIDE_YTICKS_BY_PANEL = {
    "B.": (0.6,),
    "D.": (0.8,),
}


def generate_bar_plot(
    target_dir: Path,
    output_path: Path,
    *,
    hide_yticks: Optional[Sequence[float]] = None,
) -> None:
    ok = plot_csp_classification_bars_simple_legend(
        target_dir,
        output_path,
        ylabel=DEFAULT_YLABEL,
        xlabel=None,
        hide_yticks=hide_yticks,
    )
    if not ok:
        raise RuntimeError(f"Failed to write CSP classification bars: {output_path}")


def render_pymol_assets(
    target: Fig2Target,
    target_dir: Path,
    mask_png: Path,
    binding_png: Path,
    confusion_png: Path,
    *,
    repo_root: Path,
) -> None:
    pml = _required_pml(target_dir)
    view_path = _find_saved_view(target.holo_pdb, target_dir.name)
    view = _load_view(view_path)
    print(f"[FIG2] Using PyMOL view {view_path} for {target.holo_pdb}")

    jobs = (
        (pml["mask"], mask_png),
        (pml["binding"], binding_png),
        (pml["confusion"], confusion_png),
    )
    old_cwd = os.getcwd()
    try:
        os.chdir(repo_root)
        for pml_path, out_png in jobs:
            out_png.parent.mkdir(parents=True, exist_ok=True)
            render_pymol_panel_with_view(str(pml_path), view, str(out_png))
    finally:
        os.chdir(old_cwd)


def collect_assets(
    outputs_dir: Path,
    assets_dir: Path,
    *,
    skip_pymol: bool,
    repo_root: Path,
) -> List[TargetAssets]:
    outputs_index, bmrb_cache = build_resolution_caches(outputs_dir)
    collected: List[TargetAssets] = []
    for target in FIG2_TARGETS:
        target_dir = _resolve_target_dir(target, outputs_index, bmrb_cache)
        bars, mask, binding, confusion = _asset_paths(assets_dir, target)
        print(f"[FIG2] Panel {target.panel}: {target_dir}")
        generate_bar_plot(
            target_dir,
            bars,
            hide_yticks=_HIDE_YTICKS_BY_PANEL.get(target.panel),
        )
        if skip_pymol:
            missing = [str(p) for p in (mask, binding, confusion) if not p.is_file()]
            if missing:
                raise FileNotFoundError(
                    "--skip-pymol set but rendered panels are missing:\n"
                    + "\n".join(missing)
                )
            print(f"[FIG2] Reusing PyMOL panels under {mask.parent}")
        else:
            render_pymol_assets(
                target,
                target_dir,
                mask,
                binding,
                confusion,
                repo_root=repo_root,
            )
        collected.append(
            TargetAssets(
                target=target,
                target_dir=target_dir,
                bars=bars,
                mask=mask,
                binding=binding,
                confusion=confusion,
            )
        )
    return collected


def compose_figure_2(
    assets: Sequence[TargetAssets],
    output_path: Path,
    *,
    dpi: int,
    fig_width: float,
    fig_height: float,
) -> None:
    n_rows = len(assets)
    left, right, top, bottom = 0.05, 0.99, 0.97, 0.06
    cell_w, cell_h = _pymol_cell_size_inches(
        fig_width=fig_width,
        fig_height=fig_height,
        n_rows=n_rows,
        left=left,
        right=right,
        top=top,
        bottom=bottom,
    )

    fig = plt.figure(figsize=(fig_width, fig_height))
    grid = GridSpec(
        n_rows,
        4,
        figure=fig,
        wspace=0.02,
        hspace=0.10,
        width_ratios=list(WIDTH_RATIOS),
    )
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

    row_axes: List[List[plt.Axes]] = []
    for r, item in enumerate(assets):
        bars = _trim_whitespace(_as_float_image(mpimg.imread(str(item.bars))))
        pymol_raw = [
            _as_float_image(mpimg.imread(str(item.mask))),
            _as_float_image(mpimg.imread(str(item.binding))),
            _as_float_image(mpimg.imread(str(item.confusion))),
        ]
        angle = optimal_pymol_rotation_angle(pymol_raw[0], cell_w, cell_h)
        print(f"[FIG2] {item.target.holo_pdb} PyMOL rotation: {angle:.0f} deg")
        pymol_images = [rotate_and_trim(img, angle) for img in pymol_raw]
        images = [bars, *pymol_images]

        axes_row: List[plt.Axes] = []
        for c, image in enumerate(images):
            ax = fig.add_subplot(grid[r, c])
            ax.imshow(image, aspect="equal")
            if c == 0:
                ax.set_box_aspect(image.shape[0] / image.shape[1])
                ax.set_anchor("N")
            ax.axis("off")
            axes_row.append(ax)
        row_axes.append(axes_row)

        axes_row[0].text(
            0.0,
            0.99,
            item.target.panel,
            transform=axes_row[0].transAxes,
            ha="left",
            va="top",
            fontsize=22,
            fontweight="bold",
            color="black",
            clip_on=False,
        )

    fig.canvas.draw()
    for c, title in enumerate(COLUMN_TITLES):
        ax = row_axes[0][c]
        ax.text(
            0.5,
            1.035,
            title,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=16,
            fontweight="bold",
            color="black",
            clip_on=False,
        )
        if c > 0:
            ax.text(
                0.0,
                1.012,
                PYMOL_SUBLABELS[c - 1],
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=14,
                fontweight="bold",
                color="black",
                clip_on=False,
            )

    pos = row_axes[-1][0].get_position()
    fig.text(
        (pos.x0 + pos.x1) / 2.0,
        pos.y0 - 0.008,
        XLABEL_COLUMN,
        ha="center",
        va="top",
        fontsize=16,
        color="black",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"[FIG2] Wrote {output_path.resolve()}")


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create Figure 2: four stacked case-study rows (CSP bars + PyMOL panels)."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=_REPO_ROOT / "outputs",
        help="Pipeline outputs root (default: outputs).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_REPO_ROOT / "figures" / "figure_2.png",
        help="Composite PNG path (default: figures/figure_2.png).",
    )
    parser.add_argument(
        "--assets-dir",
        type=Path,
        default=_REPO_ROOT / "figures" / "figure_2_assets",
        help="Directory for per-target bar plots and PyMOL renders.",
    )
    parser.add_argument(
        "--skip-pymol",
        action="store_true",
        help="Reuse existing PyMOL panel PNGs under --assets-dir.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--fig-width", type=float, default=18.0)
    parser.add_argument("--fig-height", type=float, default=20.0)
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    output_path = args.output if args.output.is_absolute() else _REPO_ROOT / args.output
    assets_dir = args.assets_dir if args.assets_dir.is_absolute() else _REPO_ROOT / args.assets_dir

    if not outputs_dir.is_dir():
        print(f"Error: outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1

    assets = collect_assets(
        outputs_dir,
        assets_dir,
        skip_pymol=args.skip_pymol,
        repo_root=_REPO_ROOT,
    )
    compose_figure_2(
        assets,
        output_path,
        dpi=args.dpi,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
