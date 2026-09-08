#!/usr/bin/env python3
"""Offline FP% sensitivity: exclude residue classes from confusion-matrix counts.

Recomputes the published per-target FP rates from existing master_alignment.csv
files (no pipeline rerun), under:

1. baseline (all recorded-CSP residues)
2. exclude histidine
3. exclude His/Lys/Arg/Glu/Asp and literal N/C termini (min/max holo_resi)

Metrics match ``create_fp_percent_distribution.py``:
  fp_pct         = 100 × FP / (TP+FP+TN+FN)
  fp_of_sig_pct  = 100 × FP / (FP+TP)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

import pandas as pd

try:
    from .analyze_targets import (
        AlignmentParsingError,
        PREDICTOR_COLUMNS,
        load_alignment,
    )
    from .create_fp_percent_distribution import (
        _load_targets_meta,
        _pair_label,
        summarize,
    )
    from .target_resolution import load_target_rows, resolve_target_rows
except ImportError:
    from analyze_targets import (  # type: ignore
        AlignmentParsingError,
        PREDICTOR_COLUMNS,
        load_alignment,
    )
    from create_fp_percent_distribution import (  # type: ignore
        _load_targets_meta,
        _pair_label,
        summarize,
    )
    from target_resolution import load_target_rows, resolve_target_rows  # type: ignore

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SIGNIFICANT_COL = "significant"

_AA_1LETTER = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

SCENARIOS = ("baseline", "no_his", "no_charged_termini")


def normalize_aa(raw: object) -> str:
    """Map holo_aa to a single uppercase 1-letter code (or '')."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return ""
    text = str(raw).strip().upper()
    if not text:
        return ""
    if len(text) == 1:
        return text
    return _AA_1LETTER.get(text, text)


def _confusion_counts(df: pd.DataFrame) -> Optional[dict]:
    """TP/FP/TN/FN and both FP rates for a residue table."""
    if _SIGNIFICANT_COL not in df.columns:
        return None
    available = [c for c in PREDICTOR_COLUMNS if c in df.columns]
    if not available:
        return None
    total = int(len(df))
    if total == 0:
        return None

    predicted = df[available].any(axis=1)
    actual = df[_SIGNIFICANT_COL]
    n_fp = int((actual & ~predicted).sum())
    n_tp = int((actual & predicted).sum())
    n_fn = int((~actual & predicted).sum())
    n_tn = total - n_tp - n_fp - n_fn
    sig = n_tp + n_fp
    return {
        "n_tp": n_tp,
        "n_fp": n_fp,
        "n_tn": n_tn,
        "n_fn": n_fn,
        "n_total": total,
        "n_significant": sig,
        "fp_pct": 100.0 * n_fp / total,
        "fp_of_sig_pct": (100.0 * n_fp / sig) if sig > 0 else float("nan"),
    }


def _holo_termini_mask(df: pd.DataFrame) -> pd.Series:
    """True for residues at min/max numeric holo_resi in this table."""
    if "holo_resi" not in df.columns:
        return pd.Series(False, index=df.index)
    resi = pd.to_numeric(df["holo_resi"], errors="coerce")
    valid = resi.notna()
    if not valid.any():
        return pd.Series(False, index=df.index)
    lo = float(resi[valid].min())
    hi = float(resi[valid].max())
    return valid & ((resi == lo) | (resi == hi))


def apply_exclusion(df: pd.DataFrame, scenario: str) -> pd.DataFrame:
    """Return a copy of df with excluded residues removed."""
    if scenario == "baseline":
        return df
    aa = df["holo_aa"].map(normalize_aa) if "holo_aa" in df.columns else pd.Series("", index=df.index)
    if scenario == "no_his":
        keep = aa != "H"
        return df.loc[keep].copy()
    if scenario == "no_charged_termini":
        charged: Set[str] = {"H", "K", "R", "E", "D"}
        keep = ~aa.isin(charged) & ~_holo_termini_mask(df)
        return df.loc[keep].copy()
    raise ValueError(f"Unknown scenario: {scenario}")


def compute_scenarios_for_target(target_dir: Path) -> Optional[dict]:
    """Baseline + exclusion FP metrics for one target directory."""
    alignment_path = target_dir / "master_alignment.csv"
    if not alignment_path.is_file():
        return None
    try:
        df = load_alignment(alignment_path)
    except AlignmentParsingError as exc:
        print(f"[WARN] Skipping {alignment_path}: {exc}", file=sys.stderr)
        return None
    if "holo_aa" not in df.columns:
        print(f"[WARN] Skipping {alignment_path}: missing holo_aa", file=sys.stderr)
        return None

    record: dict = {"system_id": target_dir.name}
    baseline_total = len(df)
    for scenario in SCENARIOS:
        filtered = apply_exclusion(df, scenario)
        counts = _confusion_counts(filtered)
        if counts is None:
            return None
        prefix = scenario
        for key, value in counts.items():
            record[f"{prefix}_{key}"] = value
        record[f"{prefix}_n_excluded"] = baseline_total - int(counts["n_total"])
        if scenario != "baseline":
            # Residues dropped that were FP under the unfiltered matrix
            base_pred = df[[c for c in PREDICTOR_COLUMNS if c in df.columns]].any(axis=1)
            base_actual = df[_SIGNIFICANT_COL]
            was_fp = base_actual & ~base_pred
            dropped_idx = df.index.difference(filtered.index)
            record[f"{prefix}_n_fp_excluded"] = int(was_fp.loc[dropped_idx].sum())
    return record


def collect_exclusion_table(outputs_dir: Path, targets_csv: Path) -> pd.DataFrame:
    """Resolve targets and compute FP rates under each exclusion scenario."""
    rows = load_target_rows(targets_csv)
    target_paths = resolve_target_rows(rows, outputs_dir)
    meta = _load_targets_meta(targets_csv)

    records: List[dict] = []
    for target_dir in target_paths:
        metrics = compute_scenarios_for_target(target_dir)
        if metrics is None:
            continue
        info = meta.get(target_dir.name, {})
        apo_pdb = info.get("apo_pdb", "")
        apo_bmrb = info.get("apo_bmrb", "")
        holo_pdb = info.get("holo_pdb", target_dir.name.split("_")[0])
        if not apo_bmrb and "_" in target_dir.name:
            apo_bmrb = target_dir.name.split("_", 1)[1]
        records.append(
            {
                **metrics,
                "apo_pdb": apo_pdb,
                "apo_bmrb": apo_bmrb,
                "holo_pdb": holo_pdb,
                "pair_label": _pair_label(apo_pdb, apo_bmrb, holo_pdb),
            }
        )

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values("baseline_fp_pct", ascending=False).reset_index(drop=True)


def _print_summary(df: pd.DataFrame) -> None:
    """Print mean/median FP rates and deltas vs baseline."""
    labels = {
        "baseline": "Baseline",
        "no_his": "Exclude His",
        "no_charged_termini": "Exclude H/K/R/E/D + holo termini",
    }
    print("FP% residue-exclusion sensitivity (offline; primary significant cutoff)")
    print("=" * 72)
    print(f"n systems: {len(df)}")
    print()

    base_fp = summarize(df["baseline_fp_pct"])
    base_sig = summarize(df["baseline_fp_of_sig_pct"])

    for scenario in SCENARIOS:
        fp = summarize(df[f"{scenario}_fp_pct"])
        sig = summarize(df[f"{scenario}_fp_of_sig_pct"])
        print(labels[scenario])
        print("-" * 72)
        print(
            f"  all-residue FP%     mean={fp['mean']:6.2f}  median={fp['median']:6.2f}  "
            f"(n={int(fp['n'])})"
        )
        print(
            f"  among-sig FP%       mean={sig['mean']:6.2f}  median={sig['median']:6.2f}  "
            f"(n={int(sig['n'])})"
        )
        if scenario != "baseline":
            print(
                f"  Δ all-residue FP%   mean={fp['mean'] - base_fp['mean']:+6.2f}  "
                f"median={fp['median'] - base_fp['median']:+6.2f}"
            )
            print(
                f"  Δ among-sig FP%     mean={sig['mean'] - base_sig['mean']:+6.2f}  "
                f"median={sig['median'] - base_sig['median']:+6.2f}"
            )
            n_excl = int(df[f"{scenario}_n_excluded"].sum())
            n_fp_excl = int(df[f"{scenario}_n_fp_excluded"].sum())
            print(f"  residues excluded (sum over systems): {n_excl}")
            print(f"  of which were baseline FP:            {n_fp_excl}")
        print()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Recompute per-target FP% after excluding His or charged+holo-termini "
            "residues from confusion-matrix counts (offline)."
        )
    )
    p.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help=(
            "Root with per-target master_alignment.csv "
            "(default: %(default)s)."
        ),
    )
    p.add_argument(
        "--targets-csv",
        type=Path,
        default=Path("data/CSP_UBQ_ph0.5_temp5C.csv"),
        help="Targets CSV with apo_bmrb / holo_pdb (default: %(default)s).",
    )
    p.add_argument(
        "--output-csv",
        type=Path,
        default=Path("figures") / "fp_percent_exclusion_sensitivity.csv",
        help="Per-target rates under each scenario (default: %(default)s).",
    )
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir if args.outputs_dir.is_absolute() else _REPO_ROOT / args.outputs_dir
    targets_csv = args.targets_csv if args.targets_csv.is_absolute() else _REPO_ROOT / args.targets_csv
    out_csv = args.output_csv if args.output_csv.is_absolute() else _REPO_ROOT / args.output_csv

    if not targets_csv.is_file():
        print(f"Error: targets CSV not found: {targets_csv}", file=sys.stderr)
        return 1

    df = collect_exclusion_table(outputs_dir, targets_csv)
    if df.empty:
        print("No receptors with primary-cutoff classification data found.", file=sys.stderr)
        return 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    _print_summary(df)
    print(f"Wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
