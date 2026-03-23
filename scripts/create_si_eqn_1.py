#!/usr/bin/env python3
"""
SI Eqn. S1 — N/H CSP formula.

Renders the N/H chemical shift perturbation formula as a PNG image.
Formula: CSP = sqrt(1/2*(dH^2+(0.14*dN)^2))

Output: figures/SE1_nh_csp.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create SI Eqn. S1 (N/H CSP formula)."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figures") / "SE1_nh_csp.png",
        help="Output PNG path (default: figures/SE1_nh_csp.png).",
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    output_path = args.output if args.output.is_absolute() else project_root / args.output

    fig, ax = plt.subplots(figsize=(8, 2))
    ax.axis("off")

    eqn = (
        r"$\mathrm{CSP} = \sqrt{\frac{1}{2}\left(\Delta\delta_\mathrm{H}^2 + \left(0.14\,\Delta\delta_\mathrm{N}\right)^2\right)}$"
    )
    ax.text(0.5, 0.5, eqn, fontsize=24, ha="center", va="center")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        output_path,
        dpi=args.dpi,
        bbox_inches="tight",
        pad_inches=0.2,
    )
    plt.close()

    print(f"SI Eqn. S1 written to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
