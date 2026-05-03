#!/usr/bin/env python3
"""
Wrapper script to run all scripts that create supplemental figures and tables.

SI numbering (aligned with each script’s own docstring and default output paths):
  - SI Table S1: create_si_table_s1.py → figures/ST1_all_receptors.tex
  - SI Tables S2–S9 / SI Figs S1–S8: create_si_figs_s1_s9.py
        (selections: hydrolases, transferases, all-α, all-β, a+b, BET-ET, TFIIH, ubiquitin;
         figures SF1–SF8 stacked histograms plus per-selection data/CSVs — see script)
  - SI Table S10 + SI Fig. S9: create_si_st10_fig_s9_dissimilar_conditions.py
        → figures/ST10_dissimilar_apo_holo_conditions.tex, SF9_dissimilar_apo_holo_conditions.png
  - SI Fig. S10: create_si_fig_s10.py → figures/SF10_ca_inclusive.png (same 1d CA row-coverage gate as S11/S14)
  - SI Fig. S11: create_si_fig_s11.py → figures/SF11_f1_ca_vs_exclusive.png
  - SI Fig. S12: create_si_fig_s12.py → figures/SF12_nn_distance.png
  - SI Fig. S13: create_si_fig_s13.py → figures/SF13_any_atom_distance.png
  - SI Fig. S14: create_si_fig_s14.py → figures/SF14_1d_CSP_boxplot.png
        (per-target F1 scores for 1D H/N/Cα CSPs; Holm-adjusted pairwise Wilcoxon stats CSV alongside)
  - SI Fig. S16: create_si_fig_s16.py → figures/SF16_ideal_offsets.png
  - SI Fig. S17: create_si_fig_s17.py → figures/SF17_significance_threshold.png
  - SI Fig. S18: create_si_fig_s18.py → figures/SF18_f1_vs_mcc.png
  - SI Eqn. S1 / S2: create_si_eqn_1.py, create_si_eqn_2.py
        → figures/SE1_nh_csp.png, figures/SE2_nh_ca_csp.png
  - SI Tables S2–S9 LaTeX: create_custom_selection_latex_tables.py

Prerequisites:
  - Pipeline outputs under --outputs-dir (master_alignment.csv per target, etc.)
  - confusion_matrix_per_system.csv (e.g. scripts/confusion_matrix_analysis.py)
  - Domain targets CSVs for ST7–ST9 (targets_BET_ET.csv, targets_TFIIH.csv, targets_ubiquitin.csv)
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

    scripts = [
        (
            "create_si_table_s1.py",
            [
                "--csp-csv",
                "data/CSP_UBQ_ph0.5_temp5C.csv",
                "--output",
                str(Path(figures) / "ST1_all_receptors.tex"),
            ],
        ),
        ("create_si_figs_s1_s9.py", ["--outputs-dir", outputs, "--figures-dir", figures, "--csv", csv_path]),
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
                str(Path(figures) / "SF10_ca_inclusive.png"),
            ],
        ),
        ("create_si_fig_s11.py", ["--outputs-dir", outputs]),
        (
            "create_si_fig_s12.py",
            ["--outputs-dir", outputs, "--figures-dir", figures],
        ),
        (
            "create_si_fig_s13.py",
            ["--outputs-dir", outputs, "--figures-dir", figures],
        ),
        (
            "create_si_fig_s14.py",
            [
                "--outputs-dir",
                outputs,
                "--output-image",
                str(Path(figures) / "SF14_1d_CSP_boxplot.png"),
            ],
        ),
        ("create_si_fig_s16.py", ["--outputs-dir", outputs, "--figures-dir", figures]),
        (
            "create_si_fig_s17.py",
            [
                "--outputs-dir",
                outputs,
                "--output",
                str(Path(figures) / "SF17_significance_threshold.png"),
            ],
        ),
    ]
    scripts.extend(
        [
            ("create_si_eqn_1.py", ["--output", str(Path(figures) / "SE1_nh_csp.png")]),
            ("create_si_eqn_2.py", ["--output", str(Path(figures) / "SE2_nh_ca_csp.png")]),
            (
                "create_si_fig_s18.py",
                [
                    "--outputs-dir",
                    outputs,
                    "--input",
                    confusion_csv,
                    "--output-image",
                    str(Path(figures) / "SF18_f1_vs_mcc.png"),
                ],
            ),
            ("create_custom_selection_latex_tables.py", ["--figures-dir", figures]),
        ]
    )

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
