#!/usr/bin/env python3
"""
Arrange the four buffer-threshold sweep heatmaps into one 2×2 figure.

Expects PNGs written by ``scripts/sweep_buffer_thresholds.py``:

  heatmap_n.png, heatmap_mean_f1.png, heatmap_pct_allosteric.png, heatmap_fp_pct.png

Uses Pillow only (no NumPy/Matplotlib) so compositing works in minimal environments.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Pillow is required. Install with: pip install Pillow\n"
        "(It is usually installed alongside matplotlib.)"
    ) from exc

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

PANEL_FILENAMES = (
    "heatmap_n.png",
    "heatmap_mean_f1.png",
    "heatmap_pct_allosteric.png",
    "heatmap_fp_pct.png",
)


def _load_font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", size)
    except OSError:
        try:
            return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size)
        except OSError:
            return ImageFont.load_default()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs/buffer_threshold_sweep"),
        help="Directory containing the four heatmap PNGs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: <input-dir>/buffer_threshold_sweep_4panel.png).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="PNG DPI metadata for saved file (default: 200). Does not rescale pixels.",
    )
    parser.add_argument(
        "--no-panel-labels",
        action="store_true",
        help="Do not draw (a)–(d) labels on the composite.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    in_dir = args.input_dir if args.input_dir.is_absolute() else repo_root / args.input_dir
    if not in_dir.is_dir():
        print(f"Error: --input-dir not found: {in_dir}", file=sys.stderr)
        return 1

    paths = [in_dir / name for name in PANEL_FILENAMES]
    for p in paths:
        if not p.is_file():
            print(f"Error: missing expected file: {p}", file=sys.stderr)
            return 1

    out_path = args.output
    if out_path is None:
        out_path = in_dir / "buffer_threshold_sweep_4panel.png"
    elif not out_path.is_absolute():
        out_path = repo_root / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    images = [Image.open(p).convert("RGB") for p in paths]
    widths = [im.width for im in images]
    heights = [im.height for im in images]
    cell_w = max(widths[0], widths[1])
    cell_h = max(heights[0], heights[2])
    cell_w = max(cell_w, max(widths[2], widths[3]))
    cell_h = max(cell_h, max(heights[1], heights[3]))

    def fitted(im: Image.Image) -> Image.Image:
        if im.width == cell_w and im.height == cell_h:
            return im
        out = Image.new("RGB", (cell_w, cell_h), color="white")
        out.paste(im, ((cell_w - im.width) // 2, (cell_h - im.height) // 2))
        return out

    fitted_imgs = [fitted(im) for im in images]
    composite = Image.new("RGB", (cell_w * 2, cell_h * 2), color="white")
    composite.paste(fitted_imgs[0], (0, 0))
    composite.paste(fitted_imgs[1], (cell_w, 0))
    composite.paste(fitted_imgs[2], (0, cell_h))
    composite.paste(fitted_imgs[3], (cell_w, cell_h))

    if not args.no_panel_labels:
        draw = ImageDraw.Draw(composite)
        font = _load_font(max(14, cell_w // 55))
        labels = ("(a)", "(b)", "(c)", "(d)")
        offsets = ((8, 8), (cell_w + 8, 8), (8, cell_h + 8), (cell_w + 8, cell_h + 8))
        for label, (ox, oy) in zip(labels, offsets):
            draw.text((ox, oy), label, fill=(0, 0, 0), font=font)

    d = int(args.dpi)
    composite.save(out_path, format="PNG", dpi=(d, d))
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
