#!/usr/bin/env python3
"""Compare closest-peak RMSD H/N offsets vs global CSP-count grid offsets.

For selected targets, the 10 residues with the smallest raw (zero-offset) HN CSP
are used as anchors. Analytic least-squares offsets minimize per-nucleus RMSD
on those anchors and are compared to the stored global grid-search offsets
(maximize count of CSP < 0.05 ppm).
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _HAS_PLT = True
except Exception:
    _HAS_PLT = False

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.compare_terminal_anchor_offsets import (  # noqa: E402
    _anchor_holo_resi,
    _draw_csp_classification_bars,
    _draw_overlay_scheme,
    _holo_hn,
    _parse_float,
    _read_csv_rows,
    _scheme_plot_data,
    count_aligned_hn,
    load_global_offsets,
    terminal_anchor_offsets,
)
from scripts.config import Referencing  # noqa: E402
from scripts.csp import _csp_n_term  # noqa: E402

DEFAULT_TARGETS = (
    "2MRE_18610",
    "7ZEY_34724",
    "2M0U_18824",
)

DEFAULT_N_ANCHORS = 10
_ANCHOR_LABEL = "Closest-peak anchors"


def closest_peak_anchor_rows(
    rows: Sequence[Dict[str, str]],
    *,
    n: int = DEFAULT_N_ANCHORS,
) -> List[Dict[str, str]]:
    """Return the n residues with smallest zero-offset HN CSP."""
    scored: List[Tuple[float, Dict[str, str]]] = []
    for row in rows:
        h_apo = _parse_float(row.get("H_apo"))
        n_apo = _parse_float(row.get("N_apo"))
        h_holo, n_holo = _holo_hn(row)
        if None in (h_apo, n_apo, h_holo, n_holo):
            continue
        dH = h_holo - h_apo
        dN = n_holo - n_apo
        csp = math.sqrt(0.5 * (dH * dH + _csp_n_term(dN)))
        scored.append((csp, row))
    scored.sort(key=lambda item: item[0])
    return [row for _, row in scored[: max(0, n)]]


def process_target(
    tgt_dir: Path,
    *,
    n_anchors: int = DEFAULT_N_ANCHORS,
) -> Dict[str, object]:
    master_path = tgt_dir / "master_alignment.csv"
    if not master_path.is_file():
        raise FileNotFoundError(f"missing {master_path}")
    master_rows = _read_csv_rows(master_path)

    csp_path = tgt_dir / "csp_table.csv"
    shift_rows = _read_csv_rows(csp_path) if csp_path.is_file() else master_rows

    anchors = closest_peak_anchor_rows(shift_rows, n=n_anchors)
    h_anchor, n_anchor, anchor_rmsd = terminal_anchor_offsets(anchors)
    h_global, n_global, best_count = load_global_offsets(tgt_dir, shift_rows)

    cutoff = float(Referencing().grid_cutoff)
    aligned_anchor, n_hn = count_aligned_hn(shift_rows, h_anchor, n_anchor, cutoff=cutoff)
    aligned_global, _ = count_aligned_hn(shift_rows, h_global, n_global, cutoff=cutoff)

    return {
        "target": tgt_dir.name,
        "tgt_dir": tgt_dir,
        "shift_rows": shift_rows,
        "anchors": anchors,
        "n_anchors": len(anchors),
        "n_hn_pairs": n_hn,
        "H_offset_anchor": h_anchor,
        "N_offset_anchor": n_anchor,
        "anchor_csp_rmsd": anchor_rmsd,
        "aligned_count_anchor": aligned_anchor,
        "H_offset_global": h_global,
        "N_offset_global": n_global,
        "global_best_count": best_count,
        "aligned_count_global": aligned_global,
        "delta_H": h_anchor - h_global,
        "delta_N": n_anchor - n_global,
        "cutoff": cutoff,
    }


def plot_target_overlay(
    result: Dict[str, object],
    out_path: Path,
) -> None:
    """Write a 2x2 figure: HSQC overlays (top) + CSP classification bars (bottom)."""
    if not _HAS_PLT:
        raise RuntimeError("matplotlib is required to generate HSQC plots")

    data = _scheme_plot_data(result)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    _draw_overlay_scheme(
        axes[0, 0],
        pairs=data["pairs"],
        h_offset=float(data["h_a"]),
        n_offset=float(data["n_a"]),
        aligned_count=int(result["aligned_count_anchor"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Closest-peak RMSD — HSQC",
        highlight_pairs=data["highlight"],
        x_limits=data["x_limits"],
        y_limits=data["y_limits"],
        anchor_label=_ANCHOR_LABEL,
    )
    _draw_overlay_scheme(
        axes[0, 1],
        pairs=data["pairs"],
        h_offset=float(data["h_g"]),
        n_offset=float(data["n_g"]),
        aligned_count=int(result["aligned_count_global"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Global CSP-count grid — HSQC",
        highlight_pairs=None,
        x_limits=data["x_limits"],
        y_limits=data["y_limits"],
    )
    _draw_csp_classification_bars(
        axes[1, 0],
        entries=data["entries_a"],
        threshold=float(data["thr_a"]),
        h_offset=float(data["h_a"]),
        n_offset=float(data["n_a"]),
        aligned_count=int(result["aligned_count_anchor"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Closest-peak RMSD — CSP bars",
        mark_anchor_resi=data["anchor_resi"],
        anchor_label=_ANCHOR_LABEL,
    )
    _draw_csp_classification_bars(
        axes[1, 1],
        entries=data["entries_g"],
        threshold=float(data["thr_g"]),
        h_offset=float(data["h_g"]),
        n_offset=float(data["n_g"]),
        aligned_count=int(result["aligned_count_global"]),
        n_hn=int(result["n_hn_pairs"]),
        scheme_title="Global CSP-count grid — CSP bars",
        mark_anchor_resi=None,
    )

    fig.suptitle(str(result["target"]), fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_aggregate_overlay(
    results: Sequence[Dict[str, object]],
    out_path: Path,
) -> None:
    """Write a rows×4 summary: HSQC_anchor | HSQC_global | CSP_anchor | CSP_global."""
    if not _HAS_PLT:
        raise RuntimeError("matplotlib is required to generate HSQC plots")
    if not results:
        return

    n_rows = len(results)
    fig, axes = plt.subplots(n_rows, 4, figsize=(22, 4.8 * n_rows), squeeze=False)

    for row_idx, result in enumerate(results):
        data = _scheme_plot_data(result)
        target = str(result["target"])

        _draw_overlay_scheme(
            axes[row_idx, 0],
            pairs=data["pairs"],
            h_offset=float(data["h_a"]),
            n_offset=float(data["n_a"]),
            aligned_count=int(result["aligned_count_anchor"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Anchor HSQC",
            highlight_pairs=data["highlight"],
            x_limits=data["x_limits"],
            y_limits=data["y_limits"],
            legend_fontsize=6,
            anchor_label=_ANCHOR_LABEL,
        )
        _draw_overlay_scheme(
            axes[row_idx, 1],
            pairs=data["pairs"],
            h_offset=float(data["h_g"]),
            n_offset=float(data["n_g"]),
            aligned_count=int(result["aligned_count_global"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Global HSQC",
            highlight_pairs=None,
            x_limits=data["x_limits"],
            y_limits=data["y_limits"],
            legend_fontsize=6,
        )
        _draw_csp_classification_bars(
            axes[row_idx, 2],
            entries=data["entries_a"],
            threshold=float(data["thr_a"]),
            h_offset=float(data["h_a"]),
            n_offset=float(data["n_a"]),
            aligned_count=int(result["aligned_count_anchor"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Anchor CSP bars",
            mark_anchor_resi=data["anchor_resi"],
            legend_fontsize=6,
            anchor_label=_ANCHOR_LABEL,
        )
        _draw_csp_classification_bars(
            axes[row_idx, 3],
            entries=data["entries_g"],
            threshold=float(data["thr_g"]),
            h_offset=float(data["h_g"]),
            n_offset=float(data["n_g"]),
            aligned_count=int(result["aligned_count_global"]),
            n_hn=int(result["n_hn_pairs"]),
            scheme_title=f"{target} — Global CSP bars",
            mark_anchor_resi=None,
            legend_fontsize=6,
        )

    fig.suptitle(
        "Closest-peak vs global offsets: HSQC overlays and CSP classification bars",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputs",
        type=Path,
        default=_ROOT / "outputs",
        help="Per-target output root",
    )
    parser.add_argument(
        "--ids",
        nargs="*",
        default=list(DEFAULT_TARGETS),
        help="Target directory names (default: 2MRE_18610 7ZEY_34724 2M0U_18824)",
    )
    parser.add_argument(
        "--n-anchors",
        type=int,
        default=DEFAULT_N_ANCHORS,
        help=f"Number of closest-peak anchors (default: {DEFAULT_N_ANCHORS})",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Output CSV path (default: <outputs>/closest_peak_anchor_offset_compare.csv)",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip HSQC overlay figure generation",
    )
    args = parser.parse_args(argv)
    out_csv = args.out_csv or (args.outputs / "closest_peak_anchor_offset_compare.csv")

    results: List[Dict[str, object]] = []
    failures: List[Tuple[str, str]] = []
    for name in args.ids:
        tgt = args.outputs / name
        try:
            results.append(process_target(tgt, n_anchors=args.n_anchors))
        except Exception as e:
            failures.append((name, str(e)))
            print(f"FAIL {name}: {e}", flush=True)

    fieldnames = [
        "target",
        "n_anchors",
        "n_hn_pairs",
        "H_offset_anchor",
        "N_offset_anchor",
        "anchor_csp_rmsd",
        "aligned_count_anchor",
        "H_offset_global",
        "N_offset_global",
        "global_best_count",
        "aligned_count_global",
        "delta_H",
        "delta_N",
        "anchor_holo_resi",
    ]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            anchor_resi = sorted(_anchor_holo_resi(row["anchors"]))  # type: ignore[arg-type]
            writer.writerow(
                {
                    "target": row["target"],
                    "n_anchors": row["n_anchors"],
                    "n_hn_pairs": row["n_hn_pairs"],
                    "H_offset_anchor": f"{float(row['H_offset_anchor']):.6f}",
                    "N_offset_anchor": f"{float(row['N_offset_anchor']):.6f}",
                    "anchor_csp_rmsd": f"{float(row['anchor_csp_rmsd']):.6f}",
                    "aligned_count_anchor": row["aligned_count_anchor"],
                    "H_offset_global": f"{float(row['H_offset_global']):.6f}",
                    "N_offset_global": f"{float(row['N_offset_global']):.6f}",
                    "global_best_count": row["global_best_count"],
                    "aligned_count_global": row["aligned_count_global"],
                    "delta_H": f"{float(row['delta_H']):.6f}",
                    "delta_N": f"{float(row['delta_N']):.6f}",
                    "anchor_holo_resi": ";".join(str(r) for r in anchor_resi),
                }
            )

    hdr = (
        f"{'target':<14} {'nA':>3} {'nHN':>4} "
        f"{'H_anc':>8} {'N_anc':>8} {'rmsd':>8} {'alnA':>5} "
        f"{'H_glb':>8} {'N_glb':>8} {'alnG':>5} {'dH':>8} {'dN':>8}"
    )
    print(hdr, flush=True)
    print("-" * len(hdr), flush=True)
    for row in results:
        print(
            f"{row['target']:<14} {row['n_anchors']:>3} {row['n_hn_pairs']:>4} "
            f"{float(row['H_offset_anchor']):>8.4f} {float(row['N_offset_anchor']):>8.4f} "
            f"{float(row['anchor_csp_rmsd']):>8.4f} {row['aligned_count_anchor']:>5} "
            f"{float(row['H_offset_global']):>8.4f} {float(row['N_offset_global']):>8.4f} "
            f"{row['aligned_count_global']:>5} "
            f"{float(row['delta_H']):>8.4f} {float(row['delta_N']):>8.4f}",
            flush=True,
        )
        resi = sorted(_anchor_holo_resi(row["anchors"]))  # type: ignore[arg-type]
        print(f"  anchors (holo_resi): {', '.join(str(r) for r in resi)}", flush=True)
    tot_a = sum(int(r["aligned_count_anchor"]) for r in results)
    tot_g = sum(int(r["aligned_count_global"]) for r in results)
    tot_hn = sum(int(r["n_hn_pairs"]) for r in results)
    print(
        f"\nTotals (CSP < {float(Referencing().grid_cutoff):g}): "
        f"anchor={tot_a}/{tot_hn}  global={tot_g}/{tot_hn}",
        flush=True,
    )
    print(f"Wrote {out_csv}", flush=True)

    if not args.no_plots:
        if not _HAS_PLT:
            print("WARNING: matplotlib unavailable; skipping plots", flush=True)
        else:
            for row in results:
                out_png = Path(row["tgt_dir"]) / "hsqc_overlay_closest_peak_vs_global.png"  # type: ignore[arg-type]
                plot_target_overlay(row, out_png)
                print(f"Wrote {out_png}", flush=True)
            agg_png = args.outputs / "hsqc_overlay_closest_peak_vs_global_all.png"
            plot_aggregate_overlay(results, agg_png)
            print(f"Wrote {agg_png}", flush=True)

    if failures:
        print(f"Failures: {len(failures)}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
