#!/usr/bin/env python3
"""Add an immutable ROC export to each target represented by a frozen analysis.

Each outputs/<target>/<date-first-export-id>/ contains roc_curve.pdf/.png,
roc_curve.csv (absolute and relative thresholds), and provenance.json. A matching
run-level directory contains an index and source snapshots. No pipeline rerun
or overwriting is permitted; identical exports are verified and reused.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analyze_csp_thresholds import sha256
from scripts.create_si_fig_s27 import BLUE, GRAY, INK, roc_axis, set_style, validate_analysis

REPO = Path(__file__).resolve().parents[1]
STYLE_SOURCE = REPO / "scripts/create_si_fig_s27.py"


def draw_target_roc(row, curve, output_dir, *, dpi=300):
    """Draw one target's exact empirical ROC and its saved baseline point."""
    set_style()
    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.21, top=0.78)
    fig.text(0.15, 0.95, f"{row.holo_pdb} | apo {row.apo_bmrb} / holo {row.holo_bmrb}",
             fontsize=11, fontweight="bold")
    fig.text(0.15, 0.91, f"HN CSP ROC | {int(row.n_eligible)} eligible residues", fontsize=9, color=GRAY)
    fig.text(0.15, 0.87, f"Binding site: {int(row.n_positive)} | Distal: {int(row.n_negative)}", fontsize=8, color=GRAY)
    roc_axis(ax)
    auc_label = f"ROC AUC = {row.auc:.3f}" if math.isfinite(row.auc) else "AUC undefined: one structural class"
    ax.plot(curve.fpr, curve.tpr, color=BLUE, lw=1.6, label=auc_label)
    if math.isfinite(row.baseline_fpr) and math.isfinite(row.baseline_tpr):
        ax.scatter(row.baseline_fpr, row.baseline_tpr, marker="*", s=65,
                   facecolor="#bd751b", edgecolor=INK, lw=0.5, zorder=5,
                   clip_on=False, label="Saved primary rule")
    ax.legend(loc="lower right", fontsize=8)
    fig.text(0.15, 0.11, f"Primary cutoff: {row.primary_cutoff_ppm:.4f} ppm", fontsize=8)
    fig.text(0.15, 0.077,
             f"TP {int(row.baseline_tp)} | FP {int(row.baseline_fp)} | TN {int(row.baseline_tn)} | FN {int(row.baseline_fn)}",
             fontsize=8)
    fig.text(0.15, 0.043, "Fixed structural labels and eligibility; all CSP cutoffs swept.", fontsize=7, color=GRAY)
    fig.text(0.15, 0.02, "AUC measures binding-site localization of CSP magnitude.", fontsize=7, color=GRAY)
    for extension in ("pdf", "png"):
        fig.savefig(output_dir / f"roc_curve.{extension}", dpi=dpi)
    plt.close(fig)


def verify_receipt(directory, receipt_name, signature):
    """An existing export may only be read and reused when fully identical."""
    path = directory / receipt_name
    if not path.is_file():
        raise FileExistsError(f"Incomplete existing export preserved: {directory}; choose a new --export-id")
    receipt = json.loads(path.read_text())
    if receipt["signature"] != signature:
        raise FileExistsError(f"Existing export uses different inputs/code/settings: {directory}; choose a new --export-id")
    for name, expected in receipt["artifacts"].items():
        if sha256(directory / name) != expected:
            raise ValueError(f"Existing ROC artifact changed: {directory / name}")
    return receipt


def export_target_rocs(analysis_dir, outputs_dir=None, *, export_id=None, dpi=300, predecessor=None):
    analysis_dir = analysis_dir.resolve()
    analysis = validate_analysis(analysis_dir)
    outputs_dir = (outputs_dir or Path(analysis["source_outputs_dir"])).resolve()
    date = datetime.now(ZoneInfo("America/New_York")).date().isoformat()
    export_id = export_id or f"{date}_csp-roc_001"
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}_[A-Za-z0-9_-]+", export_id):
        raise ValueError("--export-id must be a date-first directory basename (YYYY-MM-DD_name)")
    if dpi < 72:
        raise ValueError("--dpi must be at least 72")
    signature = {
        "analysis_manifest_sha256": sha256(analysis_dir / "manifest.json"),
        "exporter_sha256": sha256(__file__), "style_source_sha256": sha256(STYLE_SOURCE),
        "export_id": export_id, "dpi": dpi,
        "predecessor_export": str(predecessor.resolve()) if predecessor else None,
    }
    summary = pd.read_csv(analysis_dir / "target_summary.csv")
    curves = pd.read_csv(analysis_dir / "roc_curves.csv")
    if summary.target_id.duplicated().any():
        raise ValueError("Duplicate targets in frozen summary")
    sources = {entry["target_id"]: entry for entry in analysis["sources"]}
    if set(summary.target_id) != set(sources):
        raise ValueError("Frozen source table and target summary disagree")
    plans = []
    # Validate the entire selection before writing anything. A target's current
    # table must match its original input, preventing stale curves being added
    # alongside a newly regenerated or rematched target.
    for row in summary.itertuples(index=False):
        if Path(row.target_id).name != row.target_id or row.target_id in (".", ".."):
            raise ValueError(f"Invalid target directory name: {row.target_id}")
        target_dir = outputs_dir / row.target_id
        if not (target_dir / "master_alignment.csv").is_file():
            raise FileNotFoundError(f"Target output table missing: {target_dir}")
        if sha256(target_dir / "master_alignment.csv") != sources[row.target_id]["sha256"]:
            raise ValueError(f"Target input changed since the frozen analysis: {row.target_id}")
        subset = curves.loc[curves.target_id == row.target_id].copy()
        if set(subset.score_mode) != {"absolute_ppm", "relative"}:
            raise ValueError(f"Missing score mode(s): {row.target_id}")
        curve = subset.loc[subset.score_mode == "absolute_ppm"]
        if len(curve) < 2 or not curve.tp.is_monotonic_increasing or not curve.fp.is_monotonic_increasing:
            raise ValueError(f"Invalid frozen ROC: {row.target_id}")
        if not ((curve.tp + curve.fn == row.n_positive) & (curve.fp + curve.tn == row.n_negative)).all():
            raise ValueError(f"ROC/summary counts disagree: {row.target_id}")
        destination = target_dir / export_id
        target_signature = {**signature, "target_id": row.target_id,
                            "master_alignment_sha256": sources[row.target_id]["sha256"]}
        if destination.exists():
            verify_receipt(destination, "provenance.json", target_signature)
        plans.append((row, subset, curve, destination, target_signature))
    index_dir = outputs_dir / export_id
    if index_dir.exists():
        previous = verify_receipt(index_dir, "manifest.json", signature)
        for entry in previous["target_exports"]:
            if sha256(Path(entry["provenance_path"])) != entry["provenance_sha256"]:
                raise ValueError(f"Target ROC provenance changed: {entry['target_id']}")
        if not all(destination.is_dir() for _, _, _, destination, _ in plans):
            raise ValueError("Completed export index references missing target exports")
        print(f"Verified and reused {len(plans)} target ROC exports: {index_dir / 'index.csv'}")
        return index_dir
    index_rows = []
    created = 0
    for row, subset, curve, destination, target_signature in plans:
        if not destination.exists():
            destination.mkdir(exist_ok=False)
            draw_target_roc(row, curve, destination, dpi=dpi)
            subset.to_csv(destination / "roc_curve.csv", index=False, float_format="%.12g")
            summary_row = {key: (None if pd.isna(value) else value) for key, value in row._asdict().items()}
            receipt = {
                "created_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
                "signature": target_signature, "analysis_dir": str(analysis_dir),
                "summary": summary_row,
                "curve": "inclusive CSP >= threshold; grouped ties; frozen eligibility and structural union",
                "score_modes": {"absolute_ppm": "HN CSP magnitude in ppm", "relative": "CSP / frozen target-specific primary cutoff"},
                "baseline": "saved primary significant flags; rounding disagreements retained in source analysis",
                "artifacts": {p.name: sha256(p) for p in sorted(destination.iterdir()) if p.is_file()},
            }
            (destination / "provenance.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
            created += 1
        index_rows.append({"target_id": row.target_id, "auc": row.auc,
                           "n_eligible": row.n_eligible, "n_positive": row.n_positive, "n_negative": row.n_negative,
                           **{name: str(destination / name) for name in ("roc_curve.pdf", "roc_curve.png", "roc_curve.csv", "provenance.json")}})
        if len(index_rows) % 25 == 0:
            print(f"[ROC] {len(index_rows)}/{len(plans)} targets ready", flush=True)
    index_dir.mkdir(exist_ok=False)
    pd.DataFrame(index_rows).to_csv(index_dir / "index.csv", index=False, float_format="%.12g")
    source_dir = index_dir / "source"
    source_dir.mkdir()
    for src in (Path(__file__), STYLE_SOURCE):
        shutil.copyfile(src, source_dir / src.name)
    (index_dir / "analysis_manifest.json").write_bytes((analysis_dir / "manifest.json").read_bytes())
    record = {
        "created_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "signature": signature, "analysis_dir": str(analysis_dir),
        "n_targets": len(plans), "n_new_target_exports": created,
        "target_exports": [{"target_id": row.target_id,
                            "provenance_path": str(destination / "provenance.json"),
                            "provenance_sha256": sha256(destination / "provenance.json")}
                           for row, _, _, destination, _ in plans],
        "artifacts": {str(p.relative_to(index_dir)): sha256(p) for p in sorted(index_dir.rglob("*")) if p.is_file()},
    }
    (index_dir / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Exported {len(plans)} target ROC curves ({created} new): {index_dir / 'index.csv'}")
    return index_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--outputs-dir", type=Path, help="Default: source outputs recorded in the analysis manifest")
    parser.add_argument("--export-id", help="Date-first version name; default: today's YYYY-MM-DD_csp-roc_001")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--predecessor-export", type=Path, help="Previous export manifest or directory for a new version")
    args = parser.parse_args(argv)
    export_target_rocs(args.analysis_dir, args.outputs_dir, export_id=args.export_id,
                       dpi=args.dpi, predecessor=args.predecessor_export)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
