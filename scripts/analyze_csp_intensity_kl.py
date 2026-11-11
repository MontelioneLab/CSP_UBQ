#!/usr/bin/env python3
"""Reproduce the within-target permutation test reported in SI Fig. S17D."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analyze_csp_thresholds import sha256


def kl_directions(tp_counts, fp_counts):
    tp = tp_counts + 0.5
    fp = fp_counts + 0.5
    tp /= tp.sum(axis=-1, keepdims=True)
    fp /= fp.sum(axis=-1, keepdims=True)
    log_ratio = np.log(tp) - np.log(fp)
    return np.sum(tp * log_ratio, axis=-1), np.sum(-fp * log_ratio, axis=-1)


def permutation_test(values, targets, is_tp, *, n_bins, permutations, seed):
    _, target_index = np.unique(targets, return_inverse=True)
    edges = np.geomspace(values.min(), values.max(), n_bins + 1)
    bins = np.searchsorted(edges[1:-1], values, side="right")
    tp = np.zeros((target_index.max() + 1, n_bins), dtype=np.int64)
    fp = np.zeros_like(tp)
    np.add.at(tp, (target_index[is_tp], bins[is_tp]), 1)
    np.add.at(fp, (target_index[~is_tp], bins[~is_tp]), 1)
    tp_to_fp, fp_to_tp = kl_directions(tp.sum(axis=0), fp.sum(axis=0))
    observed = (tp_to_fp + fp_to_tp) / 2

    rng = np.random.default_rng(seed)
    totals = tp + fp
    permuted_tp = np.zeros((permutations, n_bins), dtype=np.int64)
    for counts, n_tp in zip(totals, tp.sum(axis=1)):
        if n_tp == counts.sum():
            permuted_tp += counts
        elif n_tp > 0:
            # This draws exactly the bin counts from a within-target label shuffle.
            nonempty = counts > 0
            permuted_tp[:, nonempty] += rng.multivariate_hypergeometric(
                counts[nonempty], n_tp, size=permutations,
            )
    permuted_fp = totals.sum(axis=0) - permuted_tp
    null_tp_to_fp, null_fp_to_tp = kl_directions(permuted_tp, permuted_fp)
    null = (null_tp_to_fp + null_fp_to_tp) / 2
    exceedances = int(np.count_nonzero(null >= observed - 1e-12))
    return {
        "n_bins": n_bins, "permutations": permutations, "seed": seed,
        "tp_to_fp_nats": float(tp_to_fp), "fp_to_tp_nats": float(fp_to_tp),
        "symmetric_kl_nats": float(observed), "exceedances": exceedances,
        "permutation_p": (exceedances + 1) / (permutations + 1),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest_path = args.analysis_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    input_path = args.analysis_dir / "violin_residues.csv"
    input_sha256 = sha256(input_path)
    if input_sha256 != manifest["data_sha256"][input_path.name]:
        raise ValueError(f"Frozen analysis input changed: {input_path}")

    frame = pd.read_csv(input_path)
    values = frame.csp_ppm.to_numpy()
    if (
        set(frame.baseline_class) != {"TP", "FP"}
        or not frame.eligible.eq(1).all()
        or not frame.baseline_significant.eq(1).all()
        or not np.all(np.isfinite(values) & (values > 0))
    ):
        raise ValueError("Expected eligible significant TP/FP residues with positive finite CSPs.")
    targets = frame.target_id.to_numpy()
    is_tp = frame.baseline_class.eq("TP").to_numpy()
    result = {
        "created_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "threshold_analysis": os.path.relpath(args.analysis_dir, args.output.parent),
        "threshold_analysis_manifest_sha256": sha256(manifest_path),
        "violin_residues_sha256": input_sha256,
        "source_sha256": sha256(__file__),
        "numpy_version": np.__version__,
        "n_targets": int(frame.target_id.nunique()),
        "n_targets_with_both_classes": int(frame.groupby("target_id").baseline_class.nunique().eq(2).sum()),
        "log_spaced_range_ppm": [float(values.min()), float(values.max())],
        "pseudocount_per_bin": 0.5,
        "null": "Shuffle TP/FP labels within targets, holding per-target class counts fixed.",
        "groups": {
            label: {"n": len(group), "mean_ppm": float(group.mean()),
                    "sample_sd_ppm": float(group.std(ddof=1)), "median_ppm": float(group.median())}
            for label, group in frame.groupby("baseline_class").csp_ppm
        },
        "primary": permutation_test(values, targets, is_tp, n_bins=20, permutations=99_999, seed=20261111),
        "bin_sensitivity": [
            permutation_test(values, targets, is_tp, n_bins=bins, permutations=9_999, seed=seed)
            for bins, seed in ((12, 20261112), (30, 20261113))
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
