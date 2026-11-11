#!/usr/bin/env python3
"""Render selected SI figures S27/S28 or the original overview and target atlas.

Use --layout si for the ROC and pooled-intensity figures only. The original
four-panel overview and atlas remain available with --layout full (default).
All figure versions are immutable; publication adds new canonical asset names.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.integrate import trapezoid

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analyze_csp_thresholds import analyze, default_run_dir, sha256
from scripts.config import classification_colors

REPO = Path(__file__).resolve().parents[1]
COLORS = {"TP": classification_colors.TP, "FP": classification_colors.FP}
BLUE, INK, GRAY = "#245a81", "#24272b", "#71777d"
MAIN_STEM = "SF27_csp_roc_intensity"
ATLAS_STEM = "per_target_roc_violin_atlas"
SI_STEMS = ("SF27_csp_roc", "SF28_csp_intensity")


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5,
        "axes.titlesize": 9, "axes.titleweight": "bold", "axes.linewidth": 0.65,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
        "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white", "axes.spines.top": False,
        "axes.spines.right": False, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
    })


def log_density(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    logv = np.log10(values)
    if len(values) < 5 or np.ptp(logv) <= 1e-10:
        return None
    support = np.linspace(logv.min(), logv.max(), 180)
    density = gaussian_kde(logv, bw_method="scott")(support)
    return support, density / trapezoid(density, support)


def density_peak(groups):
    densities = [log_density(values) for values in groups]
    return max([float(d[1].max()) for d in densities if d is not None] or [1.0])


def violin(ax, values, position, color, *, width=0.32, rng=None, show_points=False, density_scale=1):
    """Bounded log-domain KDE; small or constant groups show their actual points.

    Each nondegenerate density has equal area in log-CSP space. No extension
    beyond observed values, no KDE for n<5, and no invented empty distribution.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    rng = np.random.default_rng(0) if rng is None else rng
    if not len(values):
        ax.text(position, 0.5, "No observations", rotation=90, ha="center", va="center",
                transform=ax.get_xaxis_transform(), fontsize=7, color=GRAY)
        return
    logv = np.log10(values)
    kde = log_density(values)
    if kde is not None:
        support, density = kde
        # Normalize truncated densities to unit area and share one scale per pair.
        extent = width * density / density_scale
        ax.fill_betweenx(10 ** support, position - extent, position + extent,
                         facecolor=color, alpha=0.65, edgecolor=INK, linewidth=0.6)
    if show_points or len(values) < 5 or np.ptp(logv) <= 1e-10:
        jitter = rng.uniform(-0.085, 0.085, len(values))
        ax.scatter(position + jitter, values, s=9, c=color, edgecolors=INK,
                   linewidths=0.25, alpha=0.7, zorder=3)
    low, median, high = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot([position, position], [low, high], color=INK, lw=2.5, solid_capstyle="round", zorder=4)
    ax.scatter([position], [median], s=16, facecolor="white", edgecolor=INK, lw=0.65, zorder=5)


def intensity_axis(ax, limits):
    ax.set_yscale("log")
    ax.set_ylim(limits)
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1], ["TP", "FP"])
    ax.set_ylabel("HN CSP magnitude (ppm)")
    ax.grid(axis="y", which="major", color="#e6e8ea", linewidth=0.5)
    ax.set_axisbelow(True)


def roc_axis(ax):
    ax.plot([0, 1], [0, 1], color=GRAY, linestyle=(0, (3, 3)), lw=0.8, zorder=0)
    ax.set(xlim=(0, 1), ylim=(0, 1.02), xlabel="False-positive rate", ylabel="True-positive rate")
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])


def validate_analysis(analysis_dir, *, check_upstream=False):
    manifest = json.loads((analysis_dir / "manifest.json").read_text())
    for name, expected in manifest["data_sha256"].items():
        if sha256(analysis_dir / name) != expected:
            raise ValueError(f"Frozen analysis file changed: {name}")
    for entry in manifest["sources"]:
        if entry.get("snapshot_path") and sha256(analysis_dir / entry["snapshot_path"]) != entry["sha256"]:
            raise ValueError(f"Frozen input snapshot changed: {entry['target_id']}")
    if check_upstream:
        entries = [manifest["target_selection"], *manifest["sources"]]
        for entry in entries:
            if sha256(Path(entry["path"])) != entry["sha256"]:
                raise ValueError("Source inputs changed; create a new run linked with --predecessor")
        for name, expected in manifest["source_sha256"].items():
            if sha256(REPO / name) != expected:
                raise ValueError("Analysis code/settings changed; create a new run linked with --predecessor")
    return manifest


def _load(analysis_dir, name):
    return pd.read_csv(analysis_dir / name, keep_default_na=True)


def draw_roc(ax, curves, macro, metrics, summary, *, legend_size=6.5):
    """Shared plotting of the original panel A and standalone SI Fig. S27."""
    roc_axis(ax)
    for _, c in curves.loc[(curves.score_mode == "relative") & (curves.target_id != "POOLED")].groupby("target_id"):
        ax.plot(c.fpr, c.tpr, color="#91979d", lw=0.4, alpha=0.24, zorder=1)
    pooled = curves.loc[(curves.score_mode == "relative") & (curves.target_id == "POOLED")]
    ax.fill_between(macro.fpr, macro.tpr_ci_low, macro.tpr_ci_high, color=INK, alpha=0.12, linewidth=0)
    m = metrics.loc[metrics.score_mode == "relative"].iloc[0]
    ax.plot(pooled.fpr, pooled.tpr, color=BLUE, lw=1.7,
            label=f"Pooled AUC {m.auc:.3f} [{m.ci_low:.3f}, {m.ci_high:.3f}]")
    m = metrics.loc[metrics.aggregation == "equal_target_mean"].iloc[0]
    ax.plot(macro.fpr, macro.mean_tpr, color=INK, lw=1.4, linestyle="--",
            label=f"Mean AUC {m.auc:.3f} [{m.ci_low:.3f}, {m.ci_high:.3f}]")
    tp, fp, tn, fn = (summary["baseline_" + k].sum() for k in ("tp", "fp", "tn", "fn"))
    ax.scatter(fp / (fp + tn), tp / (tp + fn), marker="*", s=65, c="#bd751b", edgecolors=INK,
               lw=0.5, zorder=6, label="Saved primary rule")
    ax.legend(loc="lower right", fontsize=legend_size, handlelength=1.8, labelspacing=0.45)


def draw_pooled_intensity(ax, residues, *, rng=None):
    """Shared plotting of the original panel C and standalone SI Fig. S28."""
    limits = (float(residues.csp_ppm.min()) / 1.3, float(residues.csp_ppm.max()) * 1.35)
    scale = density_peak([residues.loc[residues.baseline_class == c, "csp_ppm"] for c in ("TP", "FP")])
    for position, cls in enumerate(("TP", "FP")):
        values = residues.loc[residues.baseline_class == cls, "csp_ppm"].to_numpy()
        violin(ax, values, position, COLORS[cls], rng=rng, density_scale=scale)
        ax.text(position, 0.98, f"n = {len(values):,} residues", ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=7)
    intensity_axis(ax, limits)


def create_main(analysis_dir, output_dir, manifest):
    summary = _load(analysis_dir, "target_summary.csv")
    curves = _load(analysis_dir, "roc_curves.csv")
    macro = _load(analysis_dir, "macro_roc.csv")
    metrics = _load(analysis_dir, "aggregate_metrics.csv")
    operating = _load(analysis_dir, "threshold_operating_points.csv")
    bands = _load(analysis_dir, "threshold_bootstrap_bands.csv")
    residues = _load(analysis_dir, "violin_residues.csv")
    intensity = _load(analysis_dir, "intensity_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.10, top=0.88, hspace=0.45, wspace=0.35)
    fig.suptitle("CSP threshold sensitivity and intensity distributions", x=0.09, ha="left", y=0.985, fontsize=12, fontweight="bold")
    fig.text(0.09, 0.95, f"SI Figure S27  |  {len(summary)} matched apo/holo targets  |  {manifest['n_eligible_residues']:,} eligible residues", fontsize=8.5, color=GRAY)
    ax = axes[0, 0]
    draw_roc(ax, curves, macro, metrics, summary)
    ax.set_title("A  Binding-site localization", loc="left", pad=9)
    ax = axes[0, 1]
    tp, fp, tn, fn = (summary["baseline_" + k].sum() for k in ("tp", "fp", "tn", "fn"))
    op = operating.loc[(operating.target_id == "POOLED") & np.isfinite(operating.multiplier)].copy()
    bd = bands.loc[np.isfinite(bands.multiplier)]
    for key, color, label, style in (("fdr", COLORS["FP"], "Distal fraction: FP/(TP+FP)", "-"),
                                      ("tpr", BLUE, "Sensitivity: TP/(TP+FN)", "--")):
        ax.plot(op.multiplier, op[key], color=color, lw=1.6, linestyle=style, label=label)
        ax.fill_between(bd.multiplier, bd[key + "_ci_low"], bd[key + "_ci_high"], color=color, alpha=0.12, lw=0)
    ax.axvline(1, color=GRAY, lw=0.9, linestyle=":")
    ax.scatter([1, 1], [fp / (tp + fp), tp / (tp + fn)], s=30, marker="*", c="#bd751b", edgecolors=INK, lw=0.4, zorder=6)
    ax.set(xscale="log", xlim=(0.2, 10), ylim=(0, 1.02), xlabel="Cutoff multiplier, m", ylabel="Fraction of residues")
    ax.set_xticks([0.25, 0.5, 1, 2, 5, 10], ["0.25", "0.5", "1", "2", "5", "10"])
    ax.legend(loc="upper right", fontsize=6.7)
    ax.set_title("B  Threshold tradeoff", loc="left", pad=9)
    limits = (float(residues.csp_ppm.min()) / 1.3, float(residues.csp_ppm.max()) * 1.35)
    rng = np.random.default_rng(manifest["seed"])
    draw_pooled_intensity(axes[1, 0], residues, rng=rng)
    median_scale = density_peak([intensity.loc[(intensity["class"] == c) & (intensity.n > 0), "median_ppm"] for c in ("TP", "FP")])
    for position, cls in enumerate(("TP", "FP")):
        medians = intensity.loc[(intensity["class"] == cls) & (intensity.n > 0), "median_ppm"].to_numpy()
        violin(axes[1, 1], medians, position, COLORS[cls], rng=rng, density_scale=median_scale)
        axes[1, 1].text(position, 0.98, f"n = {len(medians)} targets", ha="center", va="top", transform=axes[1, 1].get_xaxis_transform(), fontsize=7)
    paired = intensity.pivot(index="target_id", columns="class", values="median_ppm").dropna(subset=["TP", "FP"])
    offsets = {target: rng.uniform(-0.085, 0.085) for target in summary.target_id}
    for row in paired.itertuples():
        offset = offsets[row.Index]
        axes[1, 1].plot([offset, 1 + offset], [row.TP, row.FP], color=GRAY, alpha=0.20, lw=0.5, zorder=2)
    for position, cls in enumerate(("TP", "FP")):
        points = intensity.loc[(intensity["class"] == cls) & (intensity.n > 0)]
        axes[1, 1].scatter([position + offsets[t] for t in points.target_id], points.median_ppm,
                           s=9, c=COLORS[cls], edgecolors=INK, linewidths=0.25, alpha=0.7, zorder=3)
    for ax in axes[1]:
        intensity_axis(ax, limits)
    axes[1, 0].set_title("C  Residue-level CSP intensity", loc="left", pad=9)
    axes[1, 1].set_title("D  Per-target median CSP intensity", loc="left", pad=9)
    fig.text(0.09, 0.043, "TP: binding-site CSP   |   FP: distal CSP   |   Primary groups use saved significance masks.", fontsize=7, color=GRAY)
    fig.text(0.09, 0.022, "95% intervals resample whole targets; gray ROC lines show individual targets. Intensity axes are logarithmic.", fontsize=7, color=GRAY)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(output_dir / f"{MAIN_STEM}.{ext}", dpi=600)
    plt.close(fig)
    caption = (
        f"SI Figure S27. CSP threshold sensitivity and intensity distributions for {len(summary)} "
        f"pH/temperature-matched apo/holo targets ({manifest['n_eligible_residues']:,} eligible residues). "
        "Binding-site labels use the fixed canonical structural union. Residues without finite CSPs, "
        "structural labels or saved baseline significance masks, prolines, and CSP z > 100 are excluded. "
        "A: Individual target ROC curves (gray), pooled ROC from CSP/primary cutoff (blue), and equal-target "
        "mean ROC (dashed). Mean AUC is the arithmetic mean of target AUCs; the curve is interpolated. "
        "B: Pooled distal fraction and sensitivity as each frozen target cutoff is multiplied by m. "
        "Stars in A/B use the original saved primary significance flags. C: Residue-level TP/FP CSP "
        "distributions at the saved primary rule. D: Per-target median TP/FP CSP distributions with paired "
        f"lines for {len(paired)} targets; targets lacking either class are unpaired. "
        "Violins use bounded Scott-bandwidth KDE in log10 CSP space with equal area for the displayed "
        "truncated densities; dots and thick bars show medians and interquartile ranges. Groups with <5 "
        "observations or zero variance show observations only. Axes retain ppm labels on logarithmic scales. "
        f"All intervals are pointwise/scalar percentile 95% intervals from {manifest['bootstrap_replicates']} "
        "whole-target bootstrap replicates; related receptor families may remain correlated. "
        f"Reconstructed m=1 and saved significance disagree for {int(summary.cutoff_boundary_disagreements.sum())} "
        "rounded CSP boundary observation(s), exported in residue_eligibility.csv. AUC measures binding-site "
        "localization; distal CSPs can be biologically meaningful. The companion atlas includes every target.\n"
    )
    (output_dir / "S27_caption.txt").write_text(caption)


def create_si_panels(analysis_dir, output_dir, manifest):
    """Two standalone SI figures, omitting original panels B and D."""
    summary = _load(analysis_dir, "target_summary.csv")
    curves = _load(analysis_dir, "roc_curves.csv")
    macro = _load(analysis_dir, "macro_roc.csv")
    metrics = _load(analysis_dir, "aggregate_metrics.csv")
    residues = _load(analysis_dir, "violin_residues.csv")
    captions = {}
    fig, ax = plt.subplots(figsize=(5.4, 5.1))
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.81)
    fig.text(0.14, 0.955, "Binding-site localization by CSP magnitude", fontsize=11, fontweight="bold")
    fig.text(0.14, 0.913, f"SI Figure S27  |  {len(summary)} targets  |  {manifest['n_eligible_residues']:,} eligible residues", fontsize=8.5, color=GRAY)
    draw_roc(ax, curves, macro, metrics, summary, legend_size=7.5)
    ax.set_aspect("equal", adjustable="box")
    fig.text(0.14, 0.042, "ROC score = CSP / target-specific primary cutoff", fontsize=7.5, color=GRAY)
    fig.text(0.14, 0.019, "95% intervals from whole-target bootstrap resampling", fontsize=7.5, color=GRAY)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(output_dir / f"{SI_STEMS[0]}.{ext}", dpi=600)
    plt.close(fig)
    pooled = metrics.loc[metrics.score_mode == "relative"].iloc[0]
    mean = metrics.loc[metrics.aggregation == "equal_target_mean"].iloc[0]
    captions["S27"] = (
        f"SI Fig. S27. Receiver operating characteristic analysis of CSP magnitude and binding-site membership. "
        f"ROC curves for {len(summary)} pH/temperature-matched apo/holo targets "
        f"({manifest['n_eligible_residues']:,} eligible residues) were obtained by varying the HN CSP significance "
        "threshold while keeping structural binding-site labels and residue eligibility fixed. "
        "Gray lines show individual targets. The blue curve pools all eligible residues using CSP divided "
        "by the target-specific primary cutoff, max(cleaned mean CSP, 0.05 ppm). The black dashed curve "
        "is the equal-target mean ROC on a common false-positive-rate grid; mean AUC is the arithmetic "
        f"mean of the {len(summary)} individual target AUCs. Pooled AUC is {pooled.auc:.3f} "
        f"(95% CI {pooled.ci_low:.3f}-{pooled.ci_high:.3f}); mean target AUC is {mean.auc:.3f} "
        f"(95% CI {mean.ci_low:.3f}-{mean.ci_high:.3f}). The shaded band is the pointwise 95% interval "
        f"for the mean ROC; intervals use {manifest['bootstrap_replicates']:,} whole-target bootstrap "
        "replicates. Related receptor families may remain correlated. The gray diagonal indicates chance "
        "ranking; the star marks the pooled operating point from saved primary significance masks. "
        "Residues lacking valid CSPs or structural labels and the extreme CSP-z outlier (z > 100) are excluded."
    )
    fig, ax = plt.subplots(figsize=(4.4, 5.1))
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.15, top=0.81)
    fig.text(0.20, 0.955, "Significant CSP magnitudes", fontsize=11, fontweight="bold")
    fig.text(0.20, 0.913, f"SI Figure S28  |  {len(summary)} matched apo/holo targets", fontsize=8.5, color=GRAY)
    draw_pooled_intensity(ax, residues, rng=np.random.default_rng(manifest["seed"]))
    ax.set_xticks([0, 1], ["TP\nBinding site", "FP\nDistal"])
    fig.text(0.20, 0.031, "Saved primary significance masks; logarithmic ppm axis", fontsize=7, color=GRAY)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(output_dir / f"{SI_STEMS[1]}.{ext}", dpi=600)
    plt.close(fig)
    n_tp = int((residues.baseline_class == "TP").sum())
    n_fp = int((residues.baseline_class == "FP").sum())
    captions["S28"] = (
        "SI Fig. S28. Magnitudes of significant HN chemical-shift perturbations within and outside the binding site. "
        f"Distributions pool significant CSPs across the same {len(summary)} pH/temperature-matched targets "
        f"and eligibility rules used in SI Fig. S27: {n_tp:,} binding-site residues (TP, green) and "
        f"{n_fp:,} distal residues (FP, purple). Group membership uses the saved primary significance "
        "masks, corresponding to CSP >= max(cleaned mean CSP, 0.05 ppm) for each target. "
        "These distributions contain only significant CSPs. Magnitudes are in ppm on a logarithmic axis. "
        "Violins show Scott-bandwidth kernel density estimates in log10 CSP space, bounded by the "
        "observed ranges and normalized to equal area. White dots show medians; thick black lines "
        "show interquartile ranges. Each residue has equal weight in the pooled distributions; "
        "violin width represents density. TP and FP refer to the structural binding-site definition "
        "(Methods); distal CSPs can reflect real biological responses."
    )
    for name, caption in captions.items():
        (output_dir / f"{name}_caption.txt").write_text(caption + "\n")
    (output_dir / "S27-S28_captions.txt").write_text("\n\n".join(captions.values()) + "\n")


def create_atlas(analysis_dir, output_dir, manifest):
    summary = _load(analysis_dir, "target_summary.csv").sort_values("target_id").reset_index(drop=True)
    curves = _load(analysis_dir, "roc_curves.csv")
    residues = _load(analysis_dir, "violin_residues.csv")
    curve_groups = {key: group for key, group in curves.loc[curves.score_mode == "absolute_ppm"].groupby("target_id")}
    residue_groups = {key: group for key, group in residues.groupby("target_id")}
    limits = (float(residues.csp_ppm.min()) / 1.3, float(residues.csp_ppm.max()) * 1.35)
    pages = math.ceil(len(summary) / 4) + 1
    rng = np.random.default_rng(manifest["seed"])
    with PdfPages(output_dir / f"{ATLAS_STEM}.pdf", metadata={"Title": "CSP ROC and TP/FP intensity atlas", "Author": "CSP_UBQ", "Subject": manifest["run_id"]}) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.07, 0.95, "Per-target ROC / CSP intensity atlas", fontsize=17, fontweight="bold")
        fig.text(0.07, 0.918, f"{len(summary)} targets | {manifest['n_eligible_residues']:,} eligible residues | {manifest['run_id']}", fontsize=8.5, color=GRAY)
        lines = [
            "ROC curves sweep all distinct recorded CSP magnitudes; gray diagonal denotes chance.",
            "Stars mark saved primary significance. Structural labels and eligibility stay fixed.",
            "TP/FP groups use the saved primary rule. Green = TP; purple = FP; counts are shown.",
            "Intensity axes share the same logarithmic ppm range across every target and page.",
            "White dots / thick lines show median / IQR. Sparse or constant groups show actual points.",
            "AUC measures binding-site localization. Distal CSPs can be biologically meaningful.",
            "Missing structural labels and CSP z > 100 are excluded. See S27_caption.txt for methods.",
        ]
        for i, line in enumerate(lines):
            fig.text(0.07, 0.88 - i * 0.022, line, fontsize=8)
        fig.text(0.07, 0.70, "Target index (holo PDB / apo BMRB) and PDF page", fontsize=10, fontweight="bold")
        rows_per_col = math.ceil(len(summary) / 3)
        for i, row in summary.iterrows():
            column, offset = divmod(i, rows_per_col)
            fig.text(0.07 + column * 0.30, 0.675 - offset * 0.0135,
                     f"{row.target_id.replace('_', ' / '):<19}  {2 + i // 4:>2}", fontsize=7.5, family="DejaVu Sans Mono")
        fig.text(0.5, 0.025, f"1 / {pages}", ha="center", fontsize=8, color=GRAY)
        pdf.savefig(fig)
        plt.close(fig)
        for start in range(0, len(summary), 4):
            page_rows = summary.iloc[start:start + 4]
            fig, axes = plt.subplots(4, 2, figsize=(8.27, 11.69))
            fig.subplots_adjust(left=0.10, right=0.96, top=0.89, bottom=0.06, hspace=0.70, wspace=0.38)
            fig.text(0.10, 0.965, "Per-target ROC and significant-CSP intensity", fontsize=13, fontweight="bold")
            fig.text(0.10, 0.945, "ROC: all eligible CSPs    |    Violins: saved primary TP/FP groups", fontsize=8.5, color=GRAY)
            for slot in range(4):
                if slot >= len(page_rows):
                    for ax in axes[slot]:
                        ax.axis("off")
                    continue
                row = page_rows.iloc[slot]
                roc, ax = curve_groups[row.target_id], axes[slot, 0]
                roc_axis(ax)
                ax.plot(roc.fpr, roc.tpr, color=BLUE, lw=1.3)
                if row.n_positive and row.n_negative:
                    ax.scatter(row.baseline_fpr, row.baseline_tpr, marker="*", s=40, c="#bd751b", edgecolors=INK, lw=0.4, zorder=5)
                ax.set_title(f"{row.holo_pdb} | apo {row.apo_bmrb} / holo {row.holo_bmrb}\nAUC {row.auc:.3f} | P={int(row.n_positive)}, N={int(row.n_negative)}", loc="left", fontsize=8.5, pad=6)
                ax = axes[slot, 1]
                group = residue_groups.get(row.target_id, pd.DataFrame(columns=["baseline_class", "csp_ppm"]))
                scale = density_peak([group.loc[group.baseline_class == c, "csp_ppm"] for c in ("TP", "FP")])
                for position, cls in enumerate(("TP", "FP")):
                    values = group.loc[group.baseline_class == cls, "csp_ppm"].to_numpy()
                    violin(ax, values, position, COLORS[cls], rng=rng, show_points=len(values) <= 40, density_scale=scale)
                intensity_axis(ax, limits)
                ax.set_title(f"Primary cutoff {row.primary_cutoff_ppm:.4f} ppm\nTP={int(row.baseline_tp)}, FP={int(row.baseline_fp)} | eligible n={int(row.n_eligible)}", loc="left", fontsize=8.5, pad=6)
            fig.text(0.5, 0.022, f"{2 + start // 4} / {pages}", ha="center", fontsize=8, color=GRAY)
            pdf.savefig(fig)
            plt.close(fig)


def publish_assets(output_dir, publish_dir, *, layout="full"):
    """Only add new canonical assets, or verify an identical existing publication."""
    publish_dir.mkdir(parents=True, exist_ok=True)
    names = [f"{stem}.{ext}" for stem in (SI_STEMS if layout == "si" else (MAIN_STEM,)) for ext in ("pdf", "png", "svg")]
    for name in names:
        src, dst = output_dir / name, publish_dir / name
        if dst.exists():
            if sha256(src) != sha256(dst):
                raise FileExistsError(f"Existing SI asset differs; preserved: {dst}")
        else:
            shutil.copyfile(src, dst)


def render(analysis_dir, output_dir=None, *, publish_dir=None, layout="full"):
    manifest = validate_analysis(analysis_dir)
    output_dir = output_dir or analysis_dir / ("si_figures" if layout == "si" else "figures")
    figure_manifest = output_dir / "figure_manifest.json"
    if output_dir.exists():
        if not figure_manifest.is_file():
            raise FileExistsError(f"Figure directory exists without a complete manifest: {output_dir}")
        previous = json.loads(figure_manifest.read_text())
        if previous.get("layout", "full") != layout or previous["analysis_manifest_sha256"] != sha256(analysis_dir / "manifest.json") or previous["renderer_sha256"] != sha256(__file__):
            raise FileExistsError("Existing immutable figures use different inputs/code; choose a new figure directory")
        for name, expected in previous["artifacts"].items():
            if sha256(output_dir / name) != expected:
                raise ValueError(f"Rendered artifact changed: {name}")
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
        set_style()
        if layout == "si":
            create_si_panels(analysis_dir, output_dir, manifest)
        else:
            create_main(analysis_dir, output_dir, manifest)
            create_atlas(analysis_dir, output_dir, manifest)
        shutil.copyfile(__file__, output_dir / "create_si_fig_s27.py")
        metadata = {"created_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
                    "analysis_manifest_sha256": sha256(analysis_dir / "manifest.json"),
                    "renderer_sha256": sha256(__file__), "layout": layout,
                    "predecessor_figure_manifest": str((analysis_dir / "figures/figure_manifest.json").resolve()) if layout == "si" else None,
                    "artifacts": {p.name: sha256(p) for p in sorted(output_dir.iterdir()) if p.is_file()}}
        figure_manifest.write_text(json.dumps(metadata, indent=2) + "\n")
    if publish_dir:
        publish_assets(output_dir, publish_dir, layout=layout)
    print(f"{'S27/S28' if layout == 'si' else 'Overview and atlas'}: {output_dir}")
    return output_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, help="Frozen run; if absent create a new dated run")
    parser.add_argument("--layout", choices=("full", "si"), default="full", help="si: standalone S27 ROC and S28 intensity; full: original overview/atlas")
    parser.add_argument("--outputs-dir", type=Path, help="Default: repository outputs")
    parser.add_argument("--targets-csv", type=Path, help="Default: CSP_UBQ_ph0.5_temp5C.csv")
    parser.add_argument("--figures-dir", type=Path, help="New figure directory; default <analysis-dir>/figures")
    parser.add_argument("--publish-dir", type=Path, help="Add new canonical S27 assets here; never overwrite")
    parser.add_argument("--bootstrap", type=int, help="New runs only (default 2000)")
    parser.add_argument("--seed", type=int, help="New runs only (default 271828)")
    parser.add_argument("--predecessor", type=Path)
    args = parser.parse_args(argv)
    outputs_dir = args.outputs_dir or REPO / "outputs"
    targets_csv = args.targets_csv or REPO / "data/CSP_UBQ_ph0.5_temp5C.csv"
    analysis_dir = args.analysis_dir or default_run_dir(outputs_dir)
    if not analysis_dir.exists():
        analyze(targets_csv.resolve(), outputs_dir.resolve(), analysis_dir.resolve(),
                bootstrap=args.bootstrap if args.bootstrap is not None else 2000,
                seed=args.seed if args.seed is not None else 271828, predecessor=args.predecessor)
    else:
        previous = validate_analysis(analysis_dir, check_upstream=True)
        changed = (
            (args.bootstrap is not None and args.bootstrap != previous["bootstrap_replicates"])
            or (args.seed is not None and args.seed != previous["seed"])
            or (args.targets_csv is not None and sha256(args.targets_csv) != previous["target_selection"]["sha256"])
            or (args.outputs_dir is not None and str(args.outputs_dir.resolve()) != previous["source_outputs_dir"])
            or args.predecessor is not None
        )
        if changed:
            raise ValueError("Requested inputs/settings differ; choose a new analysis directory and --predecessor")
    render(analysis_dir.resolve(), args.figures_dir.resolve() if args.figures_dir else None,
           publish_dir=args.publish_dir.resolve() if args.publish_dir else None, layout=args.layout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
