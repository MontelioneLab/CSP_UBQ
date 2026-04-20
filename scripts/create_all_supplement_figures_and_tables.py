#!/usr/bin/env python3
"""
Wrapper script to run all scripts that create supplemental figures and tables.

SI numbering (current convention):
  - SI Table S1: create_si_table_s1.py (CSP_UBQ_ph0.5_temp5C.csv)
  - SI Tables S2–S9 / SI Figs S1–S8: create_si_figs_s1_s12.py
        (EC class / SCOPe fold / domain selections; hydrolases, transferases,
         all alpha, all beta, a+b, BET-ET, TFIIH, ubiquitin)
  - SI Table S10 + SI Fig. S9: create_si_st14_fig_s13_dissimilar_conditions.py
        (CSP_UBQ rows not in CSP_UBQ_ph0.5_temp5C.csv; dissimilar apo/holo conditions)
  - SI Fig. S10: create_si_fig_s13.py (CA-inclusive confusion histograms)
  - SI Fig. S11: create_si_fig_s14.py (CA-inclusive vs exclusive F1)
  - SI Fig. S12: create_si_fig_s15.py (N–N distance histograms)
  - SI Fig. S13: create_si_fig_s16.py (atom–atom distance histograms)
  - SI Fig. S14: create_si_fig_s17.py (PDB Advanced Search; placeholder)
  - SI Eqn. S1 / S2: create_si_eqn_1.py, create_si_eqn_2.py
  - SI Fig. S15: create_si_fig_s20.py (ideal N/H offsets)
  - SI Fig. S16: create_si_fig_s21.py (CSP significance threshold histogram)
  - SI Fig. S17: create_si_fig_s22.py (F1 vs MCC)
  - SI Tables S2–S9 LaTeX tables: create_custom_selection_latex_tables.py

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
        "--include-placeholders",
        action="store_true",
        help="Run placeholder scripts (SI Fig. S14 / create_si_fig_s17.py) that do not yet produce output.",
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
        ("create_si_figs_s1_s12.py", ["--outputs-dir", outputs, "--figures-dir", figures, "--csv", csv_path]),
        (
            "create_si_st14_fig_s13_dissimilar_conditions.py",
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
        ("create_si_fig_s13.py", ["--outputs-dir", outputs]),
        ("create_si_fig_s14.py", ["--outputs-dir", outputs]),
        ("create_si_fig_s15.py", ["--outputs-dir", outputs, "--figures-dir", figures]),
        ("create_si_fig_s16.py", ["--outputs-dir", outputs, "--figures-dir", figures]),
    ]
    if args.include_placeholders:
        scripts.append(("create_si_fig_s17.py", []))
    scripts.extend(
        [
            ("create_si_eqn_1.py", ["--output", str(Path(figures) / "SE1_nh_csp.png")]),
            ("create_si_eqn_2.py", ["--output", str(Path(figures) / "SE2_nh_ca_csp.png")]),
            ("create_si_fig_s20.py", ["--outputs-dir", outputs, "--figures-dir", figures]),
            (
                "create_si_fig_s21.py",
                ["--outputs-dir", outputs, "--output", str(Path(figures) / "SF16_significance_threshold.png")],
            ),
            ("create_si_fig_s22.py", ["--output-image", str(Path(figures) / "SF17_f1_vs_mcc.png")]),
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
