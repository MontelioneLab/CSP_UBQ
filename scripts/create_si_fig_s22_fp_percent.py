#!/usr/bin/env python3
"""
Create SI Fig. S22: among-significant FP% = 100 × FP/(TP+FP) for four cohorts.

Cohorts (sizes taken from the CSVs at runtime):
  - ideal sequence-match ∩ ph0.5
  - non-ideal sequence-match ∩ ph0.5
  - original pH/temp-matched set (n=136; excludes rematch pairs later
    appended to CSP_UBQ_ph0.5_temp5C.csv)
  - full CSP_UBQ.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.create_fp_percent_distribution import (  # noqa: E402
    _plot_metric_row,
    collect_fp_pct_table,
)

# (apo_bmrb, holo_pdb.upper()) rematch pairs appended to the ph0.5 CSV.
_PH05_REMATCH_KEYS = frozenset(
    {
        ("16925", "2RSE"),
        ("19826", "2KID"),
        ("17760", "2LY4"),
        ("17769", "2LZ6"),
        ("17769", "2MBH"),
        ("36005", "2MG5"),
        ("16411", "5LVF"),
        ("16411", "5M9D"),
        ("34725", "7ZEY"),
    }
)

_PH05_N137_CSV = Path("data/CSP_UBQ_ph0.5_temp5C.csv")

DEFAULT_COHORTS: list[tuple[str, Path]] = [
    ("Ideal sequence match", Path("data/CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv")),
    ("Non-ideal sequence match", Path("data/CSP_UBQ_ph0.5_temp5C_non_ideal_sequence_match.csv")),
    ("pH/temp matched", _PH05_N137_CSV),
    ("Full CSPdb", Path("data/CSP_UBQ.csv")),
]


def _write_ph05_original_n137(src: Path) -> Path:
    """Write ph0.5 rows excluding rematch pairs (original n=136 set)."""
    with src.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        kept = []
        for row in reader:
            key = (
                (row.get("apo_bmrb") or "").strip(),
                (row.get("holo_pdb") or "").strip().upper(),
            )
            if key not in _PH05_REMATCH_KEYS:
                kept.append(row)
    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        encoding="utf-8",
        suffix="_ph05_n137.csv",
        delete=False,
    )
    with tmp:
        writer = csv.DictWriter(tmp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(kept)
    return Path(tmp.name)

# Among-significant FP% only: 100 × FP / (TP+FP)
_METRIC_COL = "fp_of_sig_pct"
_METRIC_COLOR = "#7d3c98"


def plot_multi_cohort(
    cohort_dfs: list[tuple[str, pd.DataFrame]],
    out_path: Path,
    *,
    n_extreme: int = 2,
) -> None:
    """One row per cohort: hist + strip for FP/(TP+FP) among significant CSPs."""
    n = len(cohort_dfs)
    fig, axes = plt.subplots(
        n,
        2,
        figsize=(12.5, 3.2 * n),
        squeeze=False,
        gridspec_kw={"width_ratios": [1.25, 1.15], "wspace": 0.32, "hspace": 0.42},
    )

    for i, (title, df) in enumerate(cohort_dfs):
        _plot_metric_row(
            axes[i],
            df,
            _METRIC_COL,
            n_extreme=n_extreme,
            color=_METRIC_COLOR,
            legend_loc="upper left",
        )
        axes[i][0].set_title(
            f"{title} — 100 × FP / (TP+FP)",
            fontsize=10,
            pad=6,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    ap.add_argument(
        "--output-image",
        type=Path,
        default=Path("figures") / "SF22_fp_percent_distribution_cleaned_mean.png",
    )
    ap.add_argument("--n-extreme", type=int, default=2)
    args = ap.parse_args(list(argv) if argv is not None else None)

    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO / args.outputs_dir
    out_img = args.output_image if args.output_image.is_absolute() else _REPO / args.output_image

    cohort_dfs: list[tuple[str, pd.DataFrame]] = []
    tmp_n137: Path | None = None
    try:
        for title, rel in DEFAULT_COHORTS:
            path = rel if rel.is_absolute() else _REPO / rel
            if not path.is_file():
                print(f"Error: missing cohort CSV {path}", file=sys.stderr)
                return 1
            if rel == _PH05_N137_CSV:
                tmp_n137 = _write_ph05_original_n137(path)
                path = tmp_n137
            df = collect_fp_pct_table(outputs_dir, path)
            if df.empty:
                print(f"Error: no FP data for {path}", file=sys.stderr)
                return 1
            labeled = f"{title} (n={len(df)})"
            print(f"[SF21] {labeled}: collected from {rel.name}")
            cohort_dfs.append((labeled, df))
    finally:
        if tmp_n137 is not None:
            tmp_n137.unlink(missing_ok=True)

    plot_multi_cohort(cohort_dfs, out_img, n_extreme=args.n_extreme)
    # Also keep a copy under the legacy name for convenience.
    legacy = _REPO / "figures" / "fp_percent_distribution_cleaned_mean.png"
    if out_img.resolve() != legacy.resolve():
        import shutil

        shutil.copy2(out_img, legacy)
        shutil.copy2(out_img.with_suffix(".pdf"), legacy.with_suffix(".pdf"))
    print(f"Wrote {out_img}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
