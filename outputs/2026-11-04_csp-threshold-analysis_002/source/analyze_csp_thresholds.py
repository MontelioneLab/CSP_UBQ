#!/usr/bin/env python3
"""Immutable, read-only post-processing of continuous HN CSP threshold sweeps.

Run ``python -B scripts/analyze_csp_thresholds.py --help`` for the CLI.
No upstream tables, cached offsets, or existing analysis runs are overwritten.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.config import thresholds
from scripts.csp import compute_threshold_with_outlier_removal
from scripts.merge_csv import (
    exceeds_max_classification_csp_z, is_binding_site_row, parse_optional_bool,
)
from scripts.target_resolution import build_resolution_caches, load_target_rows, resolve_row

REPO = Path(__file__).resolve().parents[1]
STRUCTURAL_FLAGS = (
    "is_occluded_occlusion", "passes_filter_distance", "has_hbond_interaction",
    "has_charge_complement_interaction", "has_pi_contact_interaction",
)
CURVE_METRICS = ("tp", "fp", "tn", "fn", "tpr", "fpr", "specificity", "precision", "fdr", "tp_fp_ratio")


def finite_number(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def structural_label(row):
    """Positive evidence establishes membership; a negative needs complete evidence.

    Reuse the canonical union, preserving its strict <2 A fallback. Missing
    flags cannot silently become negative labels in the new evaluation.
    """
    if is_binding_site_row(row):
        return True
    flags_known = all(parse_optional_bool(row.get(k)) is not None for k in STRUCTURAL_FLAGS)
    direct_known = (
        parse_optional_bool(row.get("passes_sub_2A_filter_any_atom")) is not None
        or finite_number(row.get("min_any_atom_distance")) is not None
    )
    return False if flags_known and direct_known else None


def eligibility_reason(row):
    csp = finite_number(row.get("csp_A"))
    if csp is None:
        return "missing_or_nonfinite_csp"
    if csp < 0:
        return "negative_csp"
    if (row.get("holo_aa") or "").strip().upper() == "P":
        return "proline"
    if exceeds_max_classification_csp_z(row.get("csp_z")):
        return "extreme_csp_z"
    if parse_optional_bool(row.get("significant")) is None:
        return "missing_baseline_mask"
    if structural_label(row) is None:
        return "unknown_structural_label"
    return "eligible"


def counts_and_rates(tp, fp, tn, fn):
    """Undefined ratios stay NaN, including precision when nothing is selected."""
    tp, fp, tn, fn = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in (tp, fp, tn, fn)])
    def divide(a, b):
        return np.divide(a, b, out=np.full(a.shape, np.nan), where=b != 0)
    return dict(zip(CURVE_METRICS, (
        tp, fp, tn, fn, divide(tp, tp + fn), divide(fp, fp + tn),
        divide(tn, tn + fp), divide(tp, tp + fp), divide(fp, tp + fp), divide(tp, fp),
    )))


def roc_table(labels, scores):
    """Exact empirical ROC with inclusive thresholds and tied scores grouped.

    The +inf row selects nothing; the minimum score selects everything.
    AUC is undefined for a single structural class (rates remain explicit).
    """
    labels, scores = np.asarray(labels, dtype=bool), np.asarray(scores, dtype=float)
    if labels.ndim != 1 or scores.shape != labels.shape or not len(scores):
        raise ValueError("ROC requires nonempty, matching one-dimensional labels/scores")
    if not np.isfinite(scores).all():
        raise ValueError("ROC scores must be finite")
    order = np.argsort(-scores, kind="stable")
    s, y = scores[order], labels[order]
    ends = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s) - 1]
    tp = np.r_[0, np.cumsum(y)[ends]]
    fp = np.r_[0, np.cumsum(~y)[ends]]
    table = pd.DataFrame(counts_and_rates(tp, fp, (~y).sum() - fp, y.sum() - tp))
    table.insert(0, "threshold", np.r_[np.inf, s[ends]])
    auc = float(trapezoid(table.tpr, table.fpr)) if y.any() and (~y).any() else math.nan
    return table, auc


def at_thresholds(curve, cutoffs):
    """Look up inclusive >= operating points without a residues x thresholds matrix."""
    cutoffs = np.asarray(cutoffs, dtype=float)
    indices = np.searchsorted(-curve.threshold.to_numpy()[1:], -cutoffs, side="right")
    result = curve.iloc[indices].reset_index(drop=True).copy()
    result["threshold"] = cutoffs
    return result


def interpolate_roc(curve, grid):
    """Retain vertical segments; interpolate the upper TPR at duplicate FPRs."""
    upper = curve.groupby("fpr", sort=True).tpr.max()
    return np.interp(grid, upper.index.to_numpy(), upper.to_numpy())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def default_run_dir(outputs_dir):
    date = datetime.now(ZoneInfo("America/New_York")).date().isoformat()
    for number in range(1, 10000):
        path = outputs_dir / f"{date}_csp-threshold-analysis_{number:03d}"
        if not path.exists():
            return path
    raise RuntimeError("No free analysis run name")


def read_residues(targets_csv, outputs_dir):
    targets = load_target_rows(targets_csv)
    if not targets:
        raise ValueError("Empty target selection")
    caches = build_resolution_caches(outputs_dir)
    records, summaries, inputs, seen = [], [], [], set()
    for target in targets:
        directory = resolve_row(target, *caches)
        if directory is None:
            raise ValueError(f"Unresolved target: {target}")
        if directory.name in seen:
            raise ValueError(f"Duplicate resolved target: {directory.name}")
        seen.add(directory.name)
        path = directory / "master_alignment.csv"
        if caches[1][directory] != (target.apo_bmrb, target.holo_bmrb):
            raise ValueError(f"BMRB pair mismatch: {directory.name}")
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"csp_A", "csp_z", "significant", "holo_aa", "sequential_position", *STRUCTURAL_FLAGS}
            if required - set(reader.fieldnames or []):
                raise ValueError(f"Missing columns in {path}: {required - set(reader.fieldnames or [])}")
            rows = list(reader)
        keys = [(r.get("chain"), r.get("sequential_position")) for r in rows]
        if len(set(keys)) != len(keys):
            raise ValueError(f"Duplicate alignment positions: {path}")
        values = [v for r in rows if (v := finite_number(r.get("csp_A"))) is not None and v >= 0]
        if not values:
            raise ValueError(f"No recorded CSPs: {path}")
        stats = compute_threshold_with_outlier_removal(
            values, thresholds.outlier_z_score, thresholds.significance_z_score,
            thresholds.max_outlier_iterations, thresholds.max_outlier_fraction,
        )
        cutoff = max(stats.threshold, thresholds.cutoff_05_ppm)
        if cutoff <= 0:
            raise ValueError("Relative scores require a positive primary cutoff")
        local = []
        for r in rows:
            csp, label = finite_number(r.get("csp_A")), structural_label(r)
            baseline = parse_optional_bool(r.get("significant"))
            reason = eligibility_reason(r)
            reconstructed = csp >= cutoff if csp is not None else None
            rec = {
                "target_id": directory.name, "holo_pdb": target.holo_pdb.upper(),
                "apo_bmrb": target.apo_bmrb, "holo_bmrb": target.holo_bmrb,
                "chain": r.get("chain", ""), "sequential_position": r.get("sequential_position", ""),
                "pdb_residue_number": r.get("pdb_residue_number", ""), "holo_aa": r.get("holo_aa", ""),
                "csp_ppm": csp, "csp_z": finite_number(r.get("csp_z")),
                "primary_cutoff_ppm": cutoff,
                "relative_score": csp / cutoff if csp is not None else None,
                "binding_site": int(label) if label is not None else None,
                "baseline_significant": int(baseline) if baseline is not None else None,
                "reconstructed_significant": int(reconstructed) if reconstructed is not None else None,
                "cutoff_boundary_disagreement": int(baseline is not None and reconstructed is not None and baseline != reconstructed),
                "eligible": int(reason == "eligible"), "exclusion_reason": reason,
                "baseline_class": ("TP" if label else "FP") if baseline and label is not None else (
                    ("FN" if label else "TN") if baseline is False and label is not None else ""
                ),
                "saved_classification": r.get("classification", ""),
            }
            records.append(rec)
            local.append(rec)
        eligible = [r for r in local if r["eligible"]]
        if not eligible:
            raise ValueError(f"No eligible residues: {directory.name}")
        summary = {"target_id": directory.name, "holo_pdb": target.holo_pdb.upper(),
                   "apo_bmrb": target.apo_bmrb, "holo_bmrb": target.holo_bmrb,
                   "n_rows": len(rows), "n_recorded_csp": len(values), "n_eligible": len(eligible),
                   "cleaned_mean_ppm": stats.mean, "cleaned_sd_ppm": stats.sd,
                   "primary_cutoff_ppm": cutoff, "n_cutoff_estimation_outliers": stats.outliers_removed,
                   "cutoff_boundary_disagreements": sum(r["cutoff_boundary_disagreement"] for r in local)}
        for reason in ("missing_or_nonfinite_csp", "negative_csp", "proline", "extreme_csp_z", "missing_baseline_mask", "unknown_structural_label"):
            summary["n_" + reason] = sum(r["exclusion_reason"] == reason for r in local)
        summaries.append(summary)
        inputs.append({"target_id": directory.name, "path": str(path.resolve()), "sha256": sha256(path)})
    return pd.DataFrame(records), pd.DataFrame(summaries), inputs


def _write_csv(frame, path):
    frame.to_csv(path, index=False, float_format="%.12g")


def analyze(targets_csv, outputs_dir, run_dir, *, bootstrap=2000, seed=271828, predecessor=None):
    if run_dir.exists():
        raise FileExistsError(f"Immutable run already exists: {run_dir}")
    if bootstrap < 1:
        raise ValueError("--bootstrap must be positive")
    residues, summary, inputs = read_residues(targets_csv, outputs_dir)
    eligible = residues.loc[residues.eligible == 1].copy()
    all_curves, target_curves, points, intensities = [], {}, [], []
    multipliers = np.unique(np.r_[0, np.geomspace(0.05, 20, 160), 0.5, 1, 1.5, 2, 3, 5, 10, np.inf])
    for target_id, group in eligible.groupby("target_id", sort=False):
        curve, auc = roc_table(group.binding_site, group.csp_ppm)
        relative_curve, _ = roc_table(group.binding_site, group.relative_score)
        target_curves[target_id] = (curve, relative_curve)
        idx = summary.target_id == target_id
        summary.loc[idx, "auc"] = auc
        summary.loc[idx, "n_positive"] = int(group.binding_site.sum())
        summary.loc[idx, "n_negative"] = int((group.binding_site == 0).sum())
        for mode, table in (("absolute_ppm", curve), ("relative", relative_curve)):
            all_curves.append(table.assign(target_id=target_id, score_mode=mode))
        y, pred = group.binding_site.to_numpy(dtype=bool), group.baseline_significant.to_numpy(dtype=bool)
        counts = counts_and_rates((y & pred).sum(), (~y & pred).sum(), (~y & ~pred).sum(), (y & ~pred).sum())
        summary.loc[idx, ["baseline_" + key for key in CURVE_METRICS]] = [float(counts[k]) for k in CURVE_METRICS]
        for cls in ("TP", "FP"):
            v = group.loc[group.baseline_class == cls, "csp_ppm"].to_numpy()
            rv = group.loc[group.baseline_class == cls, "relative_score"].to_numpy()
            intensities.append({"target_id": target_id, "class": cls, "n": len(v),
                                "median_ppm": float(np.median(v)) if len(v) else math.nan,
                                "q25_ppm": float(np.quantile(v, 0.25)) if len(v) else math.nan,
                                "q75_ppm": float(np.quantile(v, 0.75)) if len(v) else math.nan,
                                "median_relative": float(np.median(rv)) if len(v) else math.nan})
        op = at_thresholds(relative_curve, multipliers).rename(columns={"threshold": "multiplier"})
        op.insert(0, "target_id", target_id)
        op["cutoff_ppm"] = multipliers * float(group.primary_cutoff_ppm.iloc[0])
        points.append(op)
    pooled = {}
    for mode, column in (("absolute_ppm", "csp_ppm"), ("relative", "relative_score")):
        curve, auc = roc_table(eligible.binding_site, eligible[column])
        pooled[mode] = {"auc": auc, "curve": curve}
        all_curves.append(curve.assign(target_id="POOLED", score_mode=mode))
    grid = np.linspace(0, 1, 1001)
    valid_ids = [t for t, (c, _) in target_curves.items() if c.tpr.notna().all() and c.fpr.notna().all()]
    if not valid_ids:
        raise ValueError("No targets have both structural classes")
    interpolated = np.array([interpolate_roc(target_curves[t][0], grid) for t in valid_ids])
    macro = pd.DataFrame({"fpr": grid, "mean_tpr": interpolated.mean(axis=0)})
    # All uncertainty resamples whole target pairs; residues stay clustered.
    rng = np.random.default_rng(seed)
    ids = summary.target_id.tolist()
    weights = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=bootstrap)
    valid_indices = [ids.index(t) for t in valid_ids]
    valid_weights = weights[:, valid_indices]
    denom = valid_weights.sum(axis=1)
    boot_roc = valid_weights[denom > 0] @ interpolated / denom[denom > 0, None]
    macro["tpr_ci_low"], macro["tpr_ci_high"] = np.quantile(boot_roc, [0.025, 0.975], axis=0)
    aucs = summary.set_index("target_id").loc[valid_ids, "auc"].to_numpy()
    macro_boot_auc = valid_weights[denom > 0] @ aucs / denom[denom > 0]
    # Rank-sum contributions allow exact pooled AUC bootstrap without sorting
    # all residue scores again for every resampled set of target pairs.
    bootstrap_rows = {"replicate": np.arange(bootstrap)}
    aggregate_metrics = []
    for mode, column in (("absolute_ppm", "csp_ppm"), ("relative", "relative_score")):
        neg_scores = [np.sort(eligible.loc[(eligible.target_id == t) & (eligible.binding_site == 0), column]) for t in ids]
        pos_scores = [eligible.loc[(eligible.target_id == t) & (eligible.binding_site == 1), column].to_numpy() for t in ids]
        wins = np.array([[np.sum(np.searchsorted(n, p, side="left") + np.searchsorted(n, p, side="right")) / 2
                          for n in neg_scores] for p in pos_scores])
        numerator = np.einsum("bi,ij,bj->b", weights, wins, weights, optimize=True)
        denominator = (weights @ np.array([len(p) for p in pos_scores])) * (weights @ np.array([len(n) for n in neg_scores]))
        boot_auc = np.divide(numerator, denominator, out=np.full(bootstrap, np.nan), where=denominator > 0)
        bootstrap_rows["pooled_auc_" + mode] = boot_auc
        lo, hi = np.nanquantile(boot_auc, [0.025, 0.975])
        aggregate_metrics.append({"aggregation": "pooled", "score_mode": mode,
                                  "auc": pooled[mode]["auc"], "ci_low": lo, "ci_high": hi,
                                  "n_targets": len(ids), "n_residues": len(eligible)})
    bootstrap_rows["mean_target_auc"] = np.divide(valid_weights @ aucs, denom, out=np.full(bootstrap, np.nan), where=denom > 0)
    lo, hi = np.quantile(macro_boot_auc, [0.025, 0.975])
    aggregate_metrics.append({"aggregation": "equal_target_mean", "score_mode": "either",
                              "auc": float(aucs.mean()), "ci_low": lo, "ci_high": hi,
                              "n_targets": len(valid_ids), "n_residues": len(eligible)})
    operating = pd.concat(points, ignore_index=True)
    pooled_op = operating.groupby("multiplier", sort=True)[["tp", "fp", "tn", "fn"]].sum().reset_index()
    rates = counts_and_rates(*(pooled_op[k] for k in ("tp", "fp", "tn", "fn")))
    for key, value in rates.items():
        pooled_op[key] = value
    pooled_op.insert(0, "target_id", "POOLED")
    operating = pd.concat([operating, pooled_op], ignore_index=True)
    boot_counts = {k: weights @ np.array([p[k].to_numpy() for p in points]) for k in ("tp", "fp", "tn", "fn")}
    bands = pd.DataFrame({"multiplier": multipliers})
    for key, values in counts_and_rates(*(boot_counts[k] for k in ("tp", "fp", "tn", "fn"))).items():
        if key in ("tpr", "fdr", "precision"):
            bands[key + "_n_bootstrap_defined"] = np.isfinite(values).sum(axis=0)
            # Explicitly skip thresholds selecting nothing, keeping undefined
            # uncertainty blank without all-NaN percentile warnings.
            for suffix, q in (("ci_low", 0.025), ("ci_high", 0.975)):
                bands[key + "_" + suffix] = [np.quantile(col[np.isfinite(col)], q) if np.isfinite(col).any() else math.nan for col in values.T]
    run_dir.mkdir(parents=True, exist_ok=False)
    _write_csv(residues, run_dir / "residue_eligibility.csv")
    _write_csv(eligible.loc[eligible.baseline_class.isin(["TP", "FP"])], run_dir / "violin_residues.csv")
    _write_csv(summary, run_dir / "target_summary.csv")
    _write_csv(pd.DataFrame(intensities), run_dir / "intensity_summary.csv")
    _write_csv(pd.concat(all_curves, ignore_index=True), run_dir / "roc_curves.csv")
    _write_csv(macro, run_dir / "macro_roc.csv")
    _write_csv(operating, run_dir / "threshold_operating_points.csv")
    _write_csv(bands, run_dir / "threshold_bootstrap_bands.csv")
    _write_csv(pd.DataFrame(aggregate_metrics), run_dir / "aggregate_metrics.csv")
    _write_csv(pd.DataFrame(bootstrap_rows), run_dir / "bootstrap_auc.csv")
    summary.loc[summary.auc.isna(), ["target_id", "n_positive", "n_negative"]].to_csv(run_dir / "undefined_auc_targets.csv", index=False)
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    source_paths = [Path(__file__), REPO / "scripts/config.py",
                    REPO / "scripts/csp.py", REPO / "scripts/merge_csv.py", REPO / "scripts/target_resolution.py"]
    snapshot = run_dir / "source"
    snapshot.mkdir()
    for path in source_paths:
        (snapshot / path.name).write_bytes(path.read_bytes())
    # Retain exact input tables, so a future upstream regeneration cannot make
    # an old analysis irreproducible. Original scientific target IDs are kept.
    for entry in inputs:
        contents = Path(entry["path"]).read_bytes()
        if hashlib.sha256(contents).hexdigest() != entry["sha256"]:
            raise ValueError(f"Input changed during analysis: {entry['path']}")
        relative = Path("inputs") / entry["target_id"] / "master_alignment.csv"
        destination = run_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=False)
        destination.write_bytes(contents)
        entry["snapshot_path"] = str(relative)
    manifest = {
        "schema_version": 1, "created_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "run_id": run_dir.name, "predecessor": str(predecessor.resolve()) if predecessor else None,
        "source_outputs_dir": str(outputs_dir.resolve()), "git_head": git_head,
        "git_dirty_at_analysis": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True)),
        "target_selection": {"path": str(targets_csv.resolve()), "sha256": sha256(targets_csv)},
        "n_targets": len(ids), "n_eligible_residues": len(eligible),
        "exclusion_counts": residues.exclusion_reason.value_counts().to_dict(),
        "threshold_config": asdict(thresholds), "bootstrap_replicates": bootstrap, "seed": seed,
        "bootstrap_unit": "apo/holo target pair; repeated receptor families are not independent clusters",
        "intervals": "pointwise percentile 95%; AUC scalar percentile 95%; whole-target bootstrap",
        "score_modes": {"absolute_ppm": "csp_A", "relative": "csp_A / frozen reconstructed primary cutoff"},
        "baseline": "saved significant flags, with the same fixed eligibility and structural labels",
        "rounding_note": "Saved CSPs are rounded to four decimals. Boundary disagreements are exported, never silently corrected.",
        "macro_note": "Mean target AUC is exact; mean ROC is interpolated at 1001 FPRs with upper vertical endpoints.",
        "sources": inputs, "source_sha256": {str(p.relative_to(REPO)): sha256(p) for p in source_paths},
        "python": platform.python_version(), "versions": {p: version(p) for p in ("numpy", "scipy", "pandas", "matplotlib", "seaborn")},
    }
    (run_dir / targets_csv.name).write_bytes(targets_csv.read_bytes())
    manifest["data_sha256"] = {p.name: sha256(p) for p in sorted(run_dir.glob("*.csv"))}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(f"Analyzed {len(ids)} targets / {len(eligible)} eligible residues -> {run_dir}")
    print(pd.DataFrame(aggregate_metrics).to_string(index=False))
    return run_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets-csv", type=Path, default=REPO / "data/CSP_UBQ_ph0.5_temp5C.csv")
    parser.add_argument("--outputs-dir", type=Path, default=REPO / "outputs")
    parser.add_argument("--run-dir", type=Path, help="New immutable directory; default date-first numbered run")
    parser.add_argument("--predecessor", type=Path, help="Previous run manifest or directory, when settings/inputs change")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=271828)
    args = parser.parse_args(argv)
    analyze(args.targets_csv.resolve(), args.outputs_dir.resolve(),
            (args.run_dir or default_run_dir(args.outputs_dir)).resolve(),
            bootstrap=args.bootstrap, seed=args.seed, predecessor=args.predecessor)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
