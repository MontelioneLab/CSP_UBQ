#!/usr/bin/env python3
"""
Wrapper script to run all scripts that create supplemental figures and tables.

SI figure generators (one entry point per SI Fig. S1–S27):
  - S1–S8:  create_si_figs_s1_s9.py
  - S9:     create_si_st10_fig_s9_dissimilar_conditions.py
  - S10:    create_si_fig_s10.py  (HA/CA confusion histograms)
  - S11:    create_si_fig_s11.py  (CA-inclusive confusion histograms)
  - S12:    create_si_fig_s12.py
  - S13:    create_si_fig_s13.py
  - S14:    create_si_fig_s14.py
  - S15:    create_si_fig_s15.py
  - S16:    create_si_fig_s16.py  (static PDB search screenshot)
  - S17:    create_si_fig_s17.py
  - S18:    create_si_fig_s18.py
  - S19:    create_si_fig_s19.py
  - S20:    create_si_fig_s20.py
  - S21:    create_si_fig_s21.py  (wraps create_fig_3_thresholds.py)
  - S22:    create_si_fig_s22.py  (wraps create_si_fig_s22_fp_percent.py)
  - S23:    create_si_fig_s23.py  (wraps create_fig_3.py)
  - S24:    create_si_fig_s24.py  (apo–apo controls S24A–E)
  - S25:    create_si_fig_s25.py  (terminal-anchor vs global offsets)
  - S26:    create_si_fig_s26.py  (wraps create_si_fig_s26_csp_vs_distance.py)
  - S27:    create_si_fig_s27.py  (2FIN case-study z panel)

Also: SI Table S1 (create_si_table_s1.py), ST2–ST9 LaTeX
(create_custom_selection_latex_tables.py), equations S1–S3
(create_si_eqn_*.py), and SI merged PDF (build_si_merged_pdf.py).

Prerequisites:
  - Pair edits in data/CSP_UBQ_ph0.5_temp5C.csv synced via scripts/sync_pairs_from_ph05.py
  - Pipeline outputs under --outputs-dir (master_alignment.csv per target, etc.)
  - confusion_matrix_per_system.csv (e.g. scripts/confusion_matrix_analysis.py)
  - Domain targets CSVs for ST7–ST9 (targets_BET_ET.csv, targets_TFIIH.csv, targets_ubiquitin.csv)
  - latexmk/pdflatex and pypdf for the final SI PDF merge
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_script(script_name: str, args: list[str], repo_root: Path) -> int:
    """Run a script and return its exit code."""
    script_path = repo_root / "scripts" / script_name
    if not script_path.exists():
        print(f"[WRAPPER] Skipping {script_name} (not found)")
        return 0
    cmd = [sys.executable, str(script_path)] + args
    result = subprocess.run(cmd, cwd=str(repo_root))
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run all supplemental figure and table generation scripts."
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Outputs directory (default: outputs).",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("figures"),
        help="Figures directory (default: figures).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("data/CSP_UBQ.csv"),
        help="CSP_UBQ.csv path (default: data/CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop on first script failure (default: continue).",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    outputs = str(args.outputs_dir)
    figures = str(args.figures_dir)
    csv_path = str(args.csv)

    confusion_csv = str(Path(outputs) / "confusion_matrix_per_system.csv")
    ph05 = "data/CSP_UBQ_ph0.5_temp5C.csv"
    # Figure 3 / S12 / S13: original pH/temp-matched n=136 (excludes rematch pairs)
    ph05_n137 = str(Path(outputs) / "ph05_n137_targets.csv")
    if not Path(ph05_n137).is_file():
        ph05_n137 = ph05

    scripts = [
        (
            "create_si_table_s1.py",
            [
                "--csp-csv",
                ph05,
                "--output",
                str(Path(figures) / "ST1_all_receptors.tex"),
            ],
        ),
        # SI Figs S1–S8 (+ ST2–ST9 data)
        ("create_si_figs_s1_s9.py", ["--outputs-dir", outputs, "--figures-dir", figures, "--csv", csv_path]),
        # SI Fig. S9 (+ ST10)
        (
            "create_si_st10_fig_s9_dissimilar_conditions.py",
            [
                "--full-csp-csv",
                csv_path,
                "--outputs-dir",
                outputs,
                "--figures-dir",
                figures,
                "--confusion-csv",
                confusion_csv,
            ],
        ),
        (
            "create_si_fig_s10.py",
            [
                "--outputs-dir",
                outputs,
                "--output-image",
                str(Path(figures) / "SF10_ha_ca.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s11.py",
            [
                "--outputs-dir",
                outputs,
                "--output-image",
                str(Path(figures) / "SF11_ca_inclusive.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s12.py",
            ["--outputs-dir", outputs, "--targets-csv", ph05],
        ),
        (
            "create_si_fig_s13.py",
            [
                "--outputs-dir",
                outputs,
                "--figures-dir",
                figures,
                "--targets-csv",
                ph05_n137,
            ],
        ),
        (
            "create_si_fig_s14.py",
            [
                "--outputs-dir",
                outputs,
                "--figures-dir",
                figures,
                "--targets-csv",
                ph05_n137,
            ],
        ),
        (
            "create_si_fig_s15.py",
            [
                "--outputs-dir",
                outputs,
                "--output-image",
                str(Path(figures) / "SF15_1d_CSP_boxplot.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s16.py",
            ["--figures-dir", figures],
        ),
        (
            "create_si_fig_s17.py",
            [
                "--outputs-dir",
                outputs,
                "--sweep-out-dir",
                str(Path(outputs) / "buffer_threshold_sweep"),
                "--output",
                str(Path(figures) / "SF17_buffer_sweep.png"),
            ],
        ),
        (
            "create_si_fig_s18.py",
            [
                "--outputs-dir",
                outputs,
                "--figures-dir",
                figures,
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s19.py",
            [
                "--outputs-dir",
                outputs,
                "--output",
                str(Path(figures) / "SF19_significance_threshold.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        ("create_si_eqn_1.py", ["--output", str(Path(figures) / "SE1_nh_csp.png")]),
        ("create_si_eqn_2.py", ["--output", str(Path(figures) / "SE2_nh_ca_csp.png")]),
        ("create_si_eqn_3.py", ["--output", str(Path(figures) / "SE3_ha_ca_csp.png")]),
        (
            "create_si_fig_s20.py",
            [
                "--outputs-dir",
                outputs,
                "--input",
                confusion_csv,
                "--output-image",
                str(Path(figures) / "SF20_f1_vs_mcc.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s21.py",
            [
                "--outputs-dir",
                outputs,
                "--figures-dir",
                figures,
                "--output",
                str(Path(figures) / "SF21_figure_3_thresholds.png"),
                "--targets-csv",
                ph05,
            ],
        ),
        (
            "create_si_fig_s22.py",
            [
                "--outputs-dir",
                outputs,
                "--output-image",
                str(Path(figures) / "SF22_fp_percent_distribution_cleaned_mean.png"),
            ],
        ),
        (
            "create_si_fig_s23.py",
            [
                "--outputs-dir",
                outputs,
                "--targets-csv",
                "data/CSP_UBQ_ph0.5_temp5C_same_author_list.csv",
                "--output-a",
                str(Path(figures) / "figure_3_a_same_author_list.png"),
                "--output-b",
                str(Path(figures) / "figure_3_b_same_author_list.png"),
                "--output-combined",
                str(Path(figures) / "SF23_figure_3_combined_same_author_list.png"),
            ],
        ),
        (
            "create_si_fig_s24.py",
            ["--output-dir", str(Path(figures))],
        ),
        (
            "create_si_fig_s25.py",
            ["--output", str(Path(figures) / "SF25_terminal_anchor_vs_global.png")],
        ),
        (
            "create_si_fig_s26.py",
            ["--output", str(Path(figures) / "SF26_csp_vs_distance_panels.png")],
        ),
        (
            "create_si_fig_s27.py",
            ["--output", str(Path(figures) / "SF27_2FIN_case_study_z.png")],
        ),
        ("create_custom_selection_latex_tables.py", ["--figures-dir", figures]),
        (
            "build_si_merged_pdf.py",
            [
                "--figures-dir",
                figures,
                "--si-dir",
                str(repo_root / "SI_documents"),
            ],
        ),
    ]

    failed = []
    for script_name, script_args in scripts:
        print(f"\n[WRAPPER] Running {script_name} ...")
        code = run_script(script_name, script_args, repo_root)
        if code != 0:
            print(f"[WRAPPER] {script_name} exited with code {code}")
            failed.append(script_name)
            if args.stop_on_error:
                return code

    if failed:
        print(f"\n[WRAPPER] Failed: {', '.join(failed)}")
        return 1

    print("\n[WRAPPER] All supplemental figures and tables generated successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
