#!/usr/bin/env python3
"""
SI Fig. S16 — buffer-threshold sweep heatmaps.

Regenerates the buffer sweep metrics via functions in
``scripts/sweep_buffer_thresholds.py`` and writes:

  - ``outputs/buffer_threshold_sweep/sweep_metrics.csv``
  - ``outputs/buffer_threshold_sweep/heatmap_n.png``
  - ``outputs/buffer_threshold_sweep/heatmap_pct_allosteric.png``
  - ``figures/SF16_buffer_sweep.png``

The supplementary figure is a 2-panel layout:
  - Panel A: number of targets in subset
  - Panel B: mean FP / (TP + FP) (% allosteric CSPs)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Pillow is required. Install with: pip install Pillow"
    ) from exc

try:
    from .sweep_buffer_thresholds import (
        SWEEP_HEATMAP_TITLE_FONTSIZE,
        _save_heatmap,
        compute_sweep_metrics,
    )
except Exception:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scripts.sweep_buffer_thresholds import (  # type: ignore
        SWEEP_HEATMAP_TITLE_FONTSIZE,
        _save_heatmap,
        compute_sweep_metrics,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create SI Fig. S16 (buffer threshold sweep).")
    parser.add_argument("--csp", type=Path, default=Path("data/CSP_UBQ.csv"))
    parser.add_argument(
        "--exp",
        type=Path,
        default=Path("data/apo_holo_exp_conditions.csv"),
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--sweep-out-dir",
        type=Path,
        default=Path("outputs/buffer_threshold_sweep"),
        help="Directory for regenerated sweep CSV and panel heatmaps.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figures") / "SF16_buffer_sweep.png",
        help="Destination for the 2-panel supplementary figure.",
    )
    parser.add_argument("--ph-min", type=float, default=0.1)
    parser.add_argument("--ph-max", type=float, default=2.0)
    parser.add_argument("--ph-step", type=float, default=0.1)
    parser.add_argument("--temp-min", type=float, default=2.0)
    parser.add_argument("--temp-max", type=float, default=20.0)
    parser.add_argument("--temp-step", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--no-annot",
        action="store_true",
        help="Skip per-cell numeric labels in the regenerated heatmaps and figure panels.",
    )
    return parser.parse_args()


def _resolve_path(project_root: Path, path: Path) -> Path:
    return path if path.is_absolute() else project_root / path


def _load_font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", size)
    except OSError:
        try:
            return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size)
        except OSError:
            return ImageFont.load_default()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parent.parent

    csp_path = _resolve_path(project_root, args.csp)
    exp_path = _resolve_path(project_root, args.exp)
    outputs_dir = _resolve_path(project_root, args.outputs_dir)
    sweep_out_dir = _resolve_path(project_root, args.sweep_out_dir)
    output_path = _resolve_path(project_root, args.output)

    for label, path in (("--csp", csp_path), ("--exp", exp_path)):
        if not path.is_file():
            print(f"Error: {label} file not found: {path}", file=sys.stderr)
            return 1
    if not outputs_dir.is_dir():
        print(f"Error: --outputs-dir not found: {outputs_dir}", file=sys.stderr)
        return 1

    ph_values = np.round(
        np.arange(args.ph_min, args.ph_max + args.ph_step / 2.0, args.ph_step),
        6,
    )
    temp_values = np.round(
        np.arange(args.temp_min, args.temp_max + args.temp_step / 2.0, args.temp_step),
        6,
    )
    annotate = not args.no_annot

    sweep_out_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = compute_sweep_metrics(
        csp_path=csp_path,
        exp_path=exp_path,
        outputs_dir=outputs_dir,
        ph_values=ph_values,
        temp_values=temp_values,
    )
    csv_path = sweep_out_dir / "sweep_metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path} ({len(df)} rows).")

    _save_heatmap(
        df=df,
        value_col="n",
        out_path=sweep_out_dir / "heatmap_n.png",
        title="Number of targets in subset",
        cbar_label="n targets",
        fmt="d",
        cmap="viridis",
        annotate=annotate,
        dpi=args.dpi,
    )
    _save_heatmap(
        df=df,
        value_col="pct_allosteric",
        out_path=sweep_out_dir / "heatmap_pct_allosteric.png",
        title="Mean FP / (TP + FP)  (% allosteric CSPs)",
        cbar_label="FP / (TP + FP)",
        fmt=".3f",
        cmap="viridis",
        annotate=annotate,
        dpi=args.dpi,
    )

    panel_paths = [
        sweep_out_dir / "heatmap_n.png",
        sweep_out_dir / "heatmap_pct_allosteric.png",
    ]
    images = [Image.open(path).convert("RGB") for path in panel_paths]
    max_height = max(image.height for image in images)

    padded_images = []
    for image in images:
        if image.height == max_height:
            padded_images.append(image)
            continue
        canvas = Image.new("RGB", (image.width, max_height), color="white")
        canvas.paste(image, (0, (max_height - image.height) // 2))
        padded_images.append(canvas)

    gap = 30
    # Panel letters (A, B): slightly larger than matplotlib heatmap title at export DPI.
    dpi = float(args.dpi)
    panel_px = max(
        28,
        int(SWEEP_HEATMAP_TITLE_FONTSIZE * (dpi / 72.0) * 1.15),
    )
    font = _load_font(panel_px)
    label_margin = panel_px + 28

    total_width = sum(image.width for image in padded_images) + gap
    composite = Image.new("RGB", (total_width, max_height + label_margin), color="white")

    draw = ImageDraw.Draw(composite)
    x = 0
    for label, image in zip(("A", "B"), padded_images):
        composite.paste(image, (x, label_margin))
        draw.text(
            (x, label_margin // 2),
            label,
            fill=(0, 0, 0),
            font=font,
            anchor="lm",
        )
        x += image.width + gap

    d = int(args.dpi)
    composite.save(output_path, format="PNG", dpi=(d, d))
    print(f"SI Fig. S16 saved to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
