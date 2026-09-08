#!/usr/bin/env python3
"""
Compose SI Fig. S26: CSP vs interchain distance (two panels).

  a) P(significant CSP | CA–CA distance) stacked by CSP z-score
     (from plot_p_significant_vs_ca_distance.py)
  b) CSP z vs nearest atom–atom distance scatter (confusion-colored)
     — bottom-left panel of csp_z_vs_distance_scatter_ca_vs_any_atom_2x2.png
     (same content as csp_z_vs_any_atom_distance_scatter_max05.png)

Default sources live under outputs/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

_REPO = Path(__file__).resolve().parents[1]

DEFAULT_PANEL_A = (
    _REPO / "outputs" / "p_significant_vs_ca_distance_ph05.png"
)
# Equivalent to any-atom max(0.05, cleaned mean) scatter for the n=145
# buffer-similar set (not the all-targets scatter).
DEFAULT_PANEL_B = (
    _REPO
    / "outputs"
    / "csp_z_vs_any_atom_distance_scatter_max05_ph05.png"
)
DEFAULT_OUT = _REPO / "figures" / "SF26_csp_vs_distance_panels.png"


def compose(
    panel_a: Path,
    panel_b: Path,
    out_path: Path,
) -> None:
    img_a = mpimg.imread(str(panel_a))
    img_b = mpimg.imread(str(panel_b))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16.5, 7.2),
        gridspec_kw={"wspace": 0.08},
    )
    for ax, img, label in ((axes[0], img_a, "a"), (axes[1], img_b, "b")):
        ax.imshow(img)
        ax.axis("off")
        ax.text(
            0.01,
            1.01,
            label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=18,
            fontweight="bold",
            color="black",
            clip_on=False,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel-a", type=Path, default=DEFAULT_PANEL_A)
    ap.add_argument("--panel-b", type=Path, default=DEFAULT_PANEL_B)
    ap.add_argument(
        "--from-2x2",
        type=Path,
        default=None,
        help=(
            "Optional: crop bottom-left panel from "
            "csp_z_vs_distance_scatter_ca_vs_any_atom_2x2.png instead of --panel-b."
        ),
    )
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    panel_a = args.panel_a if args.panel_a.is_absolute() else _REPO / args.panel_a
    panel_b = args.panel_b if args.panel_b.is_absolute() else _REPO / args.panel_b
    out = args.output if args.output.is_absolute() else _REPO / args.output

    if args.from_2x2 is not None:
        src = args.from_2x2 if args.from_2x2.is_absolute() else _REPO / args.from_2x2
        if not src.is_file():
            print(f"Error: missing 2x2 figure {src}", file=sys.stderr)
            return 1
        # Approximate bottom-left quadrant (includes axis labels / title for that panel).
        import numpy as np

        quad = mpimg.imread(str(src))
        h, w = quad.shape[:2]
        # Slight inset to drop outer "c" row/col gutters from constrained_layout.
        y0, y1 = int(0.50 * h), h
        x0, x1 = 0, int(0.50 * w)
        cropped = quad[y0:y1, x0:x1]
        tmp = out.parent / "_sf25_panel_b_crop.png"
        plt.imsave(tmp, cropped)
        panel_b = tmp

    for path, name in ((panel_a, "panel a"), (panel_b, "panel b")):
        if not path.is_file():
            print(f"Error: missing {name}: {path}", file=sys.stderr)
            return 1

    compose(panel_a, panel_b, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
