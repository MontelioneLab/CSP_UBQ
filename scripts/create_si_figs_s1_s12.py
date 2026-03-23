#!/usr/bin/env python3
"""
Generate SI Figs S1–S12 (Confusion Matrix Histograms) and SI Tables S2–S13 (Data).

Naming scheme:
  SI Table S2 / SI Fig S1: Hydrolase Receptors
  SI Table S3 / SI Fig S2: Isomerase Receptors
  SI Table S4 / SI Fig S3: Oxidoreductase Receptors
  SI Table S5 / SI Fig S4: Transferase Receptors
  SI Table S6 / SI Fig S5: Translocase Receptors
  SI Table S7 / SI Fig S6: All Alpha Receptors
  SI Table S8 / SI Fig S7: All Beta Receptors
  SI Table S9 / SI Fig S8: Alpha and Beta (a+b) Receptors
  SI Table S10 / SI Fig S9: CBP Domain Receptors
  SI Table S11 / SI Fig S10: BET-ET Domain Receptors
  SI Table S12 / SI Fig S11: TFIIH Domain Receptors
  SI Table S13 / SI Fig S12: Ubiquitin Domain Receptors

For each configured selection, this script:
1) Loads holo_pdb IDs either from a targets CSV or from CSP_UBQ.csv class columns.
2) Writes/uses a targets CSV.
3) Runs analyze_targets.py for that subset.
4) Saves stacked histogram outputs to figures/:
   - <SF_ID>_<file_token>_stacked_hist.png
   - <SF_ID>_<file_token>_confusion_stacked_hist.png
5) Composes a multi-panel figure: A (left) = stacked_hist, B (right) = confusion_stacked_hist
   - <SF_ID>_<file_token>.png
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Set, TypedDict

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


class ClassSpec(TypedDict):
    column: Optional[str]
    class_label: Optional[str]
    file_token: str
    st_id: str
    sf_id: str
    display_title: str
    targets_csv: Optional[Path]


# Order and IDs match SI Tables S2–S13 and SI Figs S1–S12:
# S2/S1 Hydrolases, S3/S2 Isomerases, S4/S3 Oxidoreductases, S5/S4 Transferases, S6/S5 Translocases
# S7/S6 All Alpha, S8/S7 All Beta, S9/S8 Alpha and Beta (a+b)
# S10/S9 CBP, S11/S10 BET-ET, S12/S11 TFIIH, S13/S12 Ubiquitin
CLASS_SPECS: List[ClassSpec] = [
    {
        "column": "ec_classes",
        "class_label": "hydrolases",
        "file_token": "hydrolases",
        "st_id": "ST2",
        "sf_id": "SF1",
        "display_title": "Hydrolase Receptors",
        "targets_csv": None,
    },
    {
        "column": "ec_classes",
        "class_label": "isomerases",
        "file_token": "isomerases",
        "st_id": "ST3",
        "sf_id": "SF2",
        "display_title": "Isomerase Receptors",
        "targets_csv": None,
    },
    {
        "column": "ec_classes",
        "class_label": "oxidoreductases",
        "file_token": "oxidoreductases",
        "st_id": "ST4",
        "sf_id": "SF3",
        "display_title": "Oxidoreductase Receptors",
        "targets_csv": None,
    },
    {
        "column": "ec_classes",
        "class_label": "transferases",
        "file_token": "transferases",
        "st_id": "ST5",
        "sf_id": "SF4",
        "display_title": "Transferase Receptors",
        "targets_csv": None,
    },
    {
        "column": "ec_classes",
        "class_label": "translocases",
        "file_token": "translocases",
        "st_id": "ST6",
        "sf_id": "SF5",
        "display_title": "Translocase Receptors",
        "targets_csv": None,
    },
    {
        "column": "scope_fold_type",
        "class_label": "all alpha proteins",
        "file_token": "all_alpha_proteins",
        "st_id": "ST7",
        "sf_id": "SF6",
        "display_title": "All Alpha Receptors",
        "targets_csv": None,
    },
    {
        "column": "scope_fold_type",
        "class_label": "all beta proteins",
        "file_token": "all_beta_proteins",
        "st_id": "ST8",
        "sf_id": "SF7",
        "display_title": "All Beta Receptors",
        "targets_csv": None,
    },
    {
        "column": "scope_fold_type",
        "class_label": "alpha and beta proteins (a+b)",
        "file_token": "alpha_and_beta_proteins_a_plus_b",
        "st_id": "ST9",
        "sf_id": "SF8",
        "display_title": "Alpha and Beta (a+b) Receptors",
        "targets_csv": None,
    },
    {
        "column": None,
        "class_label": None,
        "file_token": "CBP",
        "st_id": "ST10",
        "sf_id": "SF9",
        "display_title": "CBP Domain Receptors",
        "targets_csv": Path("targets_CBP.csv"),
    },
    {
        "column": None,
        "class_label": None,
        "file_token": "BET_ET",
        "st_id": "ST11",
        "sf_id": "SF10",
        "display_title": "BET-ET Domain Receptors",
        "targets_csv": Path("targets_BET_ET.csv"),
    },
    {
        "column": None,
        "class_label": None,
        "file_token": "TFIIH",
        "st_id": "ST12",
        "sf_id": "SF11",
        "display_title": "TFIIH Domain Receptors",
        "targets_csv": Path("targets_TFIIH.csv"),
    },
    {
        "column": None,
        "class_label": None,
        "file_token": "ubiquitin",
        "st_id": "ST13",
        "sf_id": "SF12",
        "display_title": "Ubiquitin Domain Receptors",
        "targets_csv": Path("targets_ubiquitin.csv"),
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate supplementary stacked histograms by class.")
    parser.add_argument("--csv", type=Path, default=Path("CSP_UBQ.csv"))
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--aux-dir", type=Path, default=Path("outputs") / "si_figs_s1_s12_aux")
    parser.add_argument("--python", default=sys.executable)
    return parser.parse_args()


def _split_classes(raw: str) -> List[str]:
    return [part.strip().lower() for part in (raw or "").split(",") if part.strip()]


def collect_targets(csprank_csv: Path, column: str, class_label: str) -> List[str]:
    class_lower = class_label.lower()
    targets: Set[str] = set()
    with open(csprank_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            holo_pdb = (row.get("holo_pdb") or "").strip().lower()
            if not holo_pdb:
                continue
            class_tokens = _split_classes(row.get(column) or "")
            if class_lower in class_tokens:
                targets.add(holo_pdb)
    return sorted(targets)


def load_targets_csv(targets_csv: Path) -> List[str]:
    targets: Set[str] = set()
    with open(targets_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            holo_pdb = (row.get("holo_pdb") or "").strip().lower()
            if holo_pdb:
                targets.add(holo_pdb)
    return sorted(targets)


def write_targets_csv(targets_csv: Path, targets: Iterable[str]) -> None:
    targets_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(targets_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["holo_pdb"])
        for target in targets:
            writer.writerow([target])


def compose_two_panel_figure(
    panel_a_path: Path,
    panel_b_path: Path,
    output_image: Path,
) -> None:
    """Create a two-panel figure: A (left) = stacked_hist, B (right) = confusion_stacked_hist."""
    image_a = mpimg.imread(panel_a_path)
    image_b = mpimg.imread(panel_b_path)

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for ax in axes:
        ax.set_axis_off()

    axes[0].imshow(image_a)
    axes[1].imshow(image_b)

    axes[0].text(
        0.01,
        0.99,
        "A.",
        transform=axes[0].transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    axes[1].text(
        0.01,
        0.99,
        "B.",
        transform=axes[1].transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
    )

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02)
    output_image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_image, dpi=300)
    plt.close(fig)


def run_subset(
    *,
    repo_root: Path,
    python_executable: str,
    outputs_dir: Path,
    figures_dir: Path,
    aux_dir: Path,
    file_token: str,
    sf_id: str,
    targets_csv: Path,
) -> int:
    analyze_script = repo_root / "scripts" / "analyze_targets.py"
    output_image = aux_dir / f"si_figs_s1_s12_{file_token}.png"
    summary_csv = aux_dir / f"si_figs_s1_s12_{file_token}_summary.csv"
    histogram_image = aux_dir / f"si_figs_s1_s12_{file_token}_hist.png"
    stacked_histogram_image = figures_dir / f"{sf_id}_{file_token}_stacked_hist.png"
    confusion_histogram_image = figures_dir / f"{sf_id}_{file_token}_confusion_stacked_hist.png"
    summary_dir = aux_dir / f"si_figs_s1_s12_{file_token}_summary_statistics"

    cmd = [
        python_executable,
        str(analyze_script),
        "--outputs-dir",
        str(outputs_dir),
        "--targets-csv",
        str(targets_csv),
        "--output-image",
        str(output_image),
        "--summary-csv",
        str(summary_csv),
        "--histogram-image",
        str(histogram_image),
        "--stacked-histogram-image",
        str(stacked_histogram_image),
        "--confusion-matrix-stacked-histogram-image",
        str(confusion_histogram_image),
        "--summary-dir",
        str(summary_dir),
        "--no-plot-titles",
        "--axis-fontsize",
        "16",
        "--tick-fontsize",
        "14",
        "--legend-fontsize",
        "14",
        "--legend-title-fontsize",
        "14",
        "--compact-legend-labels",
    ]

    result = subprocess.run(cmd, cwd=str(repo_root))
    if result.returncode != 0:
        print(
            f"[SI_FIGS_S1_S12] analyze_targets.py failed for {sf_id} / {file_token} "
            f"(exit {result.returncode}). "
            "Ensure each target has outputs/<id>/master_alignment.csv (run the pipeline).",
            file=sys.stderr,
        )
    return result.returncode


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    csprank_csv = (repo_root / args.csv).resolve()
    outputs_dir = (repo_root / args.outputs_dir).resolve()
    figures_dir = (repo_root / args.figures_dir).resolve()
    aux_dir = (repo_root / args.aux_dir).resolve()

    figures_dir.mkdir(parents=True, exist_ok=True)
    aux_dir.mkdir(parents=True, exist_ok=True)

    for spec in CLASS_SPECS:
        st_id = spec["st_id"]
        sf_id = spec["sf_id"]
        file_token = spec["file_token"]
        display_title = spec["display_title"]

        targets_csv_override = spec["targets_csv"]
        if targets_csv_override is not None:
            targets_csv_override = (repo_root / targets_csv_override).resolve()
            if not targets_csv_override.exists():
                print(
                    f"[SI_FIGS_S1_S12] Skipping '{display_title}' "
                    f"(missing targets CSV: {targets_csv_override})."
                )
                continue
            targets = load_targets_csv(targets_csv_override)
            targets_csv = targets_csv_override
        else:
            column = spec["column"]
            class_label = spec["class_label"]
            if column is None or class_label is None:
                print(f"[SI_FIGS_S1_S12] Skipping '{display_title}' (invalid class spec).")
                continue
            targets = collect_targets(csprank_csv, column, class_label)
            targets_csv = aux_dir / f"targets_{file_token}.csv"
            write_targets_csv(targets_csv, targets)

        if not targets:
            print(f"[SI_FIGS_S1_S12] Skipping '{display_title}' (no targets found).")
            continue

        print(
            f"[SI_FIGS_S1_S12] {st_id} / {sf_id} {display_title}: "
            f"{len(targets)} targets"
        )
        rc = run_subset(
            repo_root=repo_root,
            python_executable=args.python,
            outputs_dir=outputs_dir,
            figures_dir=figures_dir,
            aux_dir=aux_dir,
            file_token=file_token,
            sf_id=sf_id,
            targets_csv=targets_csv,
        )
        if rc != 0:
            return rc

        # Compose multi-panel figure: A (left) = stacked_hist, B (right) = confusion_stacked_hist
        panel_a = figures_dir / f"{sf_id}_{file_token}_stacked_hist.png"
        panel_b = figures_dir / f"{sf_id}_{file_token}_confusion_stacked_hist.png"
        collated_path = figures_dir / f"{sf_id}_{file_token}.png"
        if panel_a.exists() and panel_b.exists():
            compose_two_panel_figure(panel_a, panel_b, collated_path)
            print(f"  -> Collated: {collated_path.name}")

    print("[SI_FIGS_S1_S12] Completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
