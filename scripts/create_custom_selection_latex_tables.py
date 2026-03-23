#!/usr/bin/env python3
"""
Create ST2-ST13 LaTeX tables for predefined custom target selections.

This script runs `scripts/create_csp_latex_table.py` for each selection CSV and
writes one `.tex` file per selection to `./figures/`. The generated `.tex` files
include the corresponding supplementary figure block with SF1-SF12 filenames.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Iterable


# Order matches SI Tables S2–S13
SELECTIONS = [
    {
        "st_id": "ST2",
        "title": "Hydrolase Receptors",
        "selection_name": "hydrolases",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_hydrolases.csv",
        "output_tex": "ST2_hydrolase_receptors.tex",
    },
    {
        "st_id": "ST3",
        "title": "Isomerase Receptors",
        "selection_name": "isomerases",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_isomerases.csv",
        "output_tex": "ST3_isomerase_receptors.tex",
    },
    {
        "st_id": "ST4",
        "title": "Oxidoreductase Receptors",
        "selection_name": "oxidoreductases",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_oxidoreductases.csv",
        "output_tex": "ST4_oxidoreductase_receptors.tex",
    },
    {
        "st_id": "ST5",
        "title": "Transferase Receptors",
        "selection_name": "transferases",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_transferases.csv",
        "output_tex": "ST5_transferase_receptors.tex",
    },
    {
        "st_id": "ST6",
        "title": "Translocase Receptors",
        "selection_name": "translocases",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_translocases.csv",
        "output_tex": "ST6_translocase_receptors.tex",
    },
    {
        "st_id": "ST7",
        "title": "All Alpha Receptors",
        "selection_name": "all_alpha_proteins",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_all_alpha_proteins.csv",
        "output_tex": "ST7_all_alpha_receptors.tex",
    },
    {
        "st_id": "ST8",
        "title": "All Beta Receptors",
        "selection_name": "all_beta_proteins",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_all_beta_proteins.csv",
        "output_tex": "ST8_all_beta_receptors.tex",
    },
    {
        "st_id": "ST9",
        "title": "Alpha and Beta (a+b) Receptors",
        "selection_name": "alpha_and_beta_proteins_a_plus_b",
        "targets_csv": Path("outputs") / "si_figs_s1_s12_aux" / "targets_alpha_and_beta_proteins_a_plus_b.csv",
        "output_tex": "ST9_alpha_and_beta_a_plus_b_receptors.tex",
    },
    {
        "st_id": "ST10",
        "title": "CBP Domain Receptors",
        "selection_name": "CBP",
        "targets_csv": Path("targets_CBP.csv"),
        "output_tex": "ST10_cbp_domain_receptors.tex",
    },
    {
        "st_id": "ST11",
        "title": "BET-ET Domain Receptors",
        "selection_name": "BET_ET",
        "targets_csv": Path("targets_BET_ET.csv"),
        "output_tex": "ST11_bet_et_domain_receptors.tex",
    },
    {
        "st_id": "ST12",
        "title": "TFIIH Domain Receptors",
        "selection_name": "TFIIH",
        "targets_csv": Path("targets_TFIIH.csv"),
        "output_tex": "ST12_tfiih_domain_receptors.tex",
    },
    {
        "st_id": "ST13",
        "title": "Ubiquitin Domain Receptors",
        "selection_name": "ubiquitin",
        "targets_csv": Path("targets_ubiquitin.csv"),
        "output_tex": "ST13_ubiquitin_domain_receptors.tex",
    },
]


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create ST2–ST13 LaTeX tables (SI Table S14 is create_si_st14_fig_s13_dissimilar_conditions.py)."
    )
    parser.add_argument(
        "--csp-csv",
        type=Path,
        default=Path("CSP_UBQ.csv"),
        help="Path to CSP_UBQ.csv (default: CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--confusion-csv",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv.",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("figures"),
        help="Directory for output .tex files (default: figures).",
    )
    parser.add_argument(
        "--python",
        type=str,
        default=sys.executable,
        help="Python executable for invoking create_csp_latex_table.py.",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parents[1]
    latex_script = repo_root / "scripts" / "create_csp_latex_table.py"
    csp_csv = args.csp_csv if args.csp_csv.is_absolute() else repo_root / args.csp_csv
    confusion_csv = args.confusion_csv if args.confusion_csv.is_absolute() else repo_root / args.confusion_csv
    figures_dir = args.figures_dir if args.figures_dir.is_absolute() else repo_root / args.figures_dir

    if not latex_script.exists():
        print(f"Error: Missing script: {latex_script}", file=sys.stderr)
        return 1
    if not csp_csv.exists():
        print(f"Error: Missing CSP CSV: {csp_csv}", file=sys.stderr)
        return 1
    if not confusion_csv.exists():
        print(f"Error: Missing confusion CSV: {confusion_csv}", file=sys.stderr)
        return 1

    figures_dir.mkdir(parents=True, exist_ok=True)

    for spec in SELECTIONS:
        targets_csv = spec["targets_csv"]
        targets_csv = targets_csv if targets_csv.is_absolute() else repo_root / targets_csv
        if not targets_csv.exists():
            print(
                f"Error: Missing targets CSV for {spec['st_id']}: {targets_csv}",
                file=sys.stderr,
            )
            return 1

        out_tex = figures_dir / spec["output_tex"]
        cmd = [
            args.python,
            str(latex_script),
            "--csp-csv",
            str(csp_csv),
            "--confusion-csv",
            str(confusion_csv),
            "--targets-csv",
            str(targets_csv),
            "--selection-name",
            str(spec["selection_name"]),
            "--output",
            str(out_tex),
        ]
        print(f"[{spec['st_id']}] Generating {out_tex.name} ({spec['title']})")
        subprocess.run(cmd, check=True, cwd=str(repo_root))

    print(f"Wrote {len(SELECTIONS)} LaTeX files to {figures_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
