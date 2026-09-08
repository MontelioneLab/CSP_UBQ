#!/usr/bin/env python3
"""Run local RCI parameter sweeps against a server reference bundle."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RCI_SCRIPT = ROOT / "scripts" / "rci_v_1c_PyNMR-STAR.py"
SERVER_RESULTS_DIR = ROOT / "server_rci_results"
DEFAULT_INPUT = SERVER_RESULTS_DIR / "31122_3.str.shifty.txt"
DEFAULT_SERVER_RCI = SERVER_RESULTS_DIR / "31122_3.str.shifty.txt_BMRB.RCI.txt"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "comparisons" / "31122_workshop"


@dataclass(frozen=True)
class SeriesRecord:
    residue_number: int
    residue_name: str
    value: float


@dataclass(frozen=True)
class SweepConfig:
    incomplete_data_use: int
    end_corr: int

    @property
    def label(self) -> str:
        return f"inc{self.incomplete_data_use}_endcorr{self.end_corr}"

    @property
    def cli_args(self) -> list[str]:
        args = ["-mpl", f"-end_corr{self.end_corr}"]
        if self.incomplete_data_use == 0:
            args.append("-no_i")
        return args


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run local RCI calculations for a SHIFTY input across "
            "incomplete-data and terminal-correction settings, then compare "
            "them to a server RCI reference."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="SHIFTY input used for the local runs.",
    )
    parser.add_argument(
        "--server-rci",
        type=Path,
        default=DEFAULT_SERVER_RCI,
        help="Server RCI text file used as the reference trace.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where summary tables, plots, and per-run workspaces are stored.",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
        help="Python interpreter used to launch the legacy RCI script.",
    )
    parser.add_argument(
        "--n-term-window",
        type=int,
        default=15,
        help="Residue window used for the N-terminal focused RMSD metric and plot.",
    )
    parser.add_argument(
        "--skip-runs",
        action="store_true",
        help="Reuse previously generated local runs under the output directory.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip matplotlib plot generation.",
    )
    return parser.parse_args()


def parse_series_file(path: Path) -> dict[int, SeriesRecord]:
    records: dict[int, SeriesRecord] = {}
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 3:
            raise ValueError(f"{path}:{line_number}: expected 'resi score resn', got {raw_line!r}")
        residue_number = int(parts[0])
        value = float(parts[1])
        residue_name = parts[2]
        records[residue_number] = SeriesRecord(residue_number, residue_name, value)
    return records


def read_shifty_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.startswith("#NUM"):
            header = stripped.split()
            continue
        if stripped.startswith("#"):
            continue
        parts = stripped.split()
        if header is None or len(parts) != len(header):
            continue
        rows.append(dict(zip(header, parts)))
    if header is None:
        raise ValueError(f"No SHIFTY header found in {path}")
    return header, rows


def shifty_float(value: str) -> float | None:
    if value in {".", "?", "None", "NONE"}:
        return None
    return float(value)


def audit_shifty_file(path: Path) -> tuple[str, list[int]]:
    header, rows = read_shifty_rows(path)
    zero_counts = {column: 0 for column in header[2:]}
    zero_examples = {column: [] for column in header[2:]}
    gly_zero_ha: list[int] = []
    gly_rows: list[tuple[int, float | None, float | None]] = []
    suspicious_zero_notes: list[str] = []

    for row in rows:
        residue_number = int(row["#NUM"])
        residue_name = row["AA"]
        parsed_values = {column: shifty_float(row[column]) for column in header[2:]}
        for column in header[2:]:
            value = parsed_values[column]
            if value == 0.0:
                zero_counts[column] += 1
                if len(zero_examples[column]) < 8:
                    zero_examples[column].append(residue_number)
        if residue_name == "G":
            ha_value = parsed_values["HA"]
            cb_value = parsed_values["CB"]
            gly_rows.append((residue_number, ha_value, cb_value))
            if ha_value == 0.0:
                gly_zero_ha.append(residue_number)
        else:
            suspicious_columns: list[str] = []
            if parsed_values["HA"] == 0.0:
                suspicious_columns.append("HA")
            if parsed_values["CA"] == 0.0:
                suspicious_columns.append("CA")
            if residue_name != "P" and parsed_values["CB"] == 0.0:
                suspicious_columns.append("CB")
            if suspicious_columns:
                suspicious_zero_notes.append(
                    f"  residue {residue_number} {residue_name}: "
                    + ", ".join(f"{column}=0.00" for column in suspicious_columns)
                )

    lines = [
        f"Input audit for {path.name}",
        "",
        f"Residues parsed: {len(rows)}",
        f"Glycine residues parsed: {len(gly_rows)}",
        f"Glycine residues with HA=0.00: {gly_zero_ha if gly_zero_ha else 'none'}",
        "",
        "Zero-valued field counts:",
    ]
    for column in header[2:]:
        examples = zero_examples[column]
        example_text = ", ".join(str(item) for item in examples) if examples else "none"
        lines.append(
            f"  {column}: {zero_counts[column]} zero entries"
            f" (examples: {example_text})"
        )
    lines.extend(
        [
            "",
            "Glycine rows:",
        ]
    )
    for residue_number, ha_value, cb_value in gly_rows:
        lines.append(
            f"  residue {residue_number}: HA={ha_value:.4f} CB={cb_value:.4f}"
            if ha_value is not None and cb_value is not None
            else f"  residue {residue_number}: HA={ha_value} CB={cb_value}"
        )
    lines.extend(
        [
            "",
            "Suspicious non-glycine zeroes:",
        ]
    )
    if suspicious_zero_notes:
        lines.extend(suspicious_zero_notes)
    else:
        lines.append("  none")
    lines.extend(
        [
            "",
            "Interpretation:",
            "  Glycine HA values are the key bug check here.",
            "  Glycine CB=0.00 is expected because glycine has no CB.",
            "  Proline and the first residue can legitimately have N/HN zeros in SHIFTY-style exports.",
        ]
    )
    return "\n".join(lines) + "\n", gly_zero_ha


def residue_rmsd(
    left: dict[int, SeriesRecord],
    right: dict[int, SeriesRecord],
    residues: list[int],
) -> float | None:
    squared_diffs: list[float] = []
    for residue_number in residues:
        if residue_number not in left or residue_number not in right:
            continue
        diff = left[residue_number].value - right[residue_number].value
        squared_diffs.append(diff * diff)
    if not squared_diffs:
        return None
    return math.sqrt(sum(squared_diffs) / len(squared_diffs))


def compute_metrics(
    server_records: dict[int, SeriesRecord],
    local_records: dict[int, SeriesRecord],
    n_term_window: int,
) -> dict[str, float | int | str]:
    common_residues = sorted(set(server_records) & set(local_records))
    server_only = sorted(set(server_records) - set(local_records))
    local_only = sorted(set(local_records) - set(server_records))
    resname_mismatch_count = sum(
        1
        for residue_number in common_residues
        if server_records[residue_number].residue_name != local_records[residue_number].residue_name
    )

    squared_diffs: list[float] = []
    max_abs_diff = -1.0
    max_abs_diff_residue = -1
    for residue_number in common_residues:
        diff = abs(server_records[residue_number].value - local_records[residue_number].value)
        squared_diffs.append(diff * diff)
        if diff > max_abs_diff:
            max_abs_diff = diff
            max_abs_diff_residue = residue_number

    rmsd = math.sqrt(sum(squared_diffs) / len(squared_diffs)) if squared_diffs else None
    n_term_residues = common_residues[:n_term_window]
    n_term_rmsd = residue_rmsd(server_records, local_records, n_term_residues)

    return {
        "common_residues": len(common_residues),
        "server_only_count": len(server_only),
        "local_only_count": len(local_only),
        "resname_mismatch_count": resname_mismatch_count,
        "rmsd": rmsd,
        "n_term_rmsd": n_term_rmsd,
        "max_abs_diff": max_abs_diff if max_abs_diff_residue >= 0 else None,
        "max_abs_diff_residue": max_abs_diff_residue if max_abs_diff_residue >= 0 else None,
        "server_only_preview": ",".join(str(item) for item in server_only[:10]),
        "local_only_preview": ",".join(str(item) for item in local_only[:10]),
    }


def ensure_clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def run_local_rci(
    config: SweepConfig,
    input_path: Path,
    output_dir: Path,
    python_executable: Path,
) -> dict[str, object]:
    run_dir = output_dir / "runs" / config.label
    ensure_clean_dir(run_dir)
    local_input = run_dir / input_path.name
    shutil.copy2(input_path, local_input)

    mpl_config_dir = run_dir / ".mplconfig"
    mpl_config_dir.mkdir(exist_ok=True)

    command = [
        str(python_executable),
        str(RCI_SCRIPT),
        "-b",
        input_path.name,
        *config.cli_args,
    ]
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(mpl_config_dir)
    process = subprocess.run(
        command,
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    (run_dir / "command.log").write_text(" ".join(command) + "\n", encoding="utf-8")
    (run_dir / "stdout.log").write_text(process.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(process.stderr, encoding="utf-8")

    rci_path = run_dir / f"{input_path.name}.RCI.txt"
    if process.returncode != 0:
        raise RuntimeError(
            f"{config.label}: local RCI run failed with exit code {process.returncode}. "
            f"See {run_dir / 'stderr.log'}"
        )
    if not rci_path.is_file():
        raise FileNotFoundError(f"{config.label}: expected RCI output not found at {rci_path}")

    return {
        "run_dir": run_dir,
        "rci_path": rci_path,
        "return_code": process.returncode,
    }


def reuse_local_rci(config: SweepConfig, input_path: Path, output_dir: Path) -> dict[str, object]:
    run_dir = output_dir / "runs" / config.label
    rci_path = run_dir / f"{input_path.name}.RCI.txt"
    if not rci_path.is_file():
        raise FileNotFoundError(
            f"{config.label}: no reusable RCI output found at {rci_path}. "
            "Run without --skip-runs first."
        )
    return {
        "run_dir": run_dir,
        "rci_path": rci_path,
        "return_code": 0,
    }


def write_summary_csv(rows: list[dict[str, object]], path: Path) -> None:
    fieldnames = [
        "label",
        "incomplete_data_use",
        "end_corr",
        "rmsd",
        "n_term_rmsd",
        "max_abs_diff",
        "max_abs_diff_residue",
        "common_residues",
        "server_only_count",
        "local_only_count",
        "resname_mismatch_count",
        "server_only_preview",
        "local_only_preview",
        "rci_path",
        "run_dir",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary_markdown(
    target_id: str,
    rows: list[dict[str, object]],
    input_path: Path,
    server_rci_path: Path,
    audit_path: Path,
    output_path: Path,
) -> None:
    sorted_rows = sorted(rows, key=lambda row: (float(row["rmsd"]), float(row["n_term_rmsd"])))
    lines = [
        f"# RCI workshop sweep for {target_id}",
        "",
        f"- Input SHIFTY file: `{input_path}`",
        f"- Server RCI reference: `{server_rci_path}`",
        f"- Input audit: `{audit_path}`",
        "",
        "## Best matches by overall RMSD",
        "",
        "| Config | incomplete_data_use | end_corr | RMSD | N-term RMSD | Max abs diff |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in sorted_rows[:8]:
        lines.append(
            "| {label} | {incomplete_data_use} | {end_corr} | {rmsd:.6f} | {n_term_rmsd:.6f} | {max_abs_diff:.6f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Lower RMSD indicates a closer match to the server RCI trace.",
            "- `N-term RMSD` is computed on the first residues shared with the server trace, so it highlights the N-terminal behavior you called out.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_results(
    target_id: str,
    server_records: dict[int, SeriesRecord],
    result_rows: list[dict[str, object]],
    output_dir: Path,
    n_term_window: int,
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    server_x = sorted(server_records)
    server_y = [server_records[residue_number].value for residue_number in server_x]

    ranked_rows = sorted(result_rows, key=lambda row: float(row["rmsd"]))
    best_rows = ranked_rows[:4]

    generated: list[Path] = []
    end_corr_colors = {
        0: "#9a3412",
        1: "#c2410c",
        2: "#2563eb",
        3: "#7c3aed",
        4: "#db2777",
        5: "#059669",
        6: "#0f766e",
    }
    incomplete_styles = {
        0: ("solid", 1.35, 0.82),
        1: ((0, (5, 3)), 1.35, 0.65),
    }

    plt.close("all")
    fig, ax = plt.subplots(figsize=(13.0, 7.0), constrained_layout=True)
    fig.patch.set_facecolor("#f7f3eb")
    ax.set_facecolor("#fffdfa")
    ax.plot(server_x, server_y, color="#111827", linewidth=3.0, label="Server RCI", zorder=6)
    for row in ranked_rows:
        local_records = parse_series_file(Path(row["rci_path"]))
        local_x = sorted(local_records)
        local_y = [local_records[residue_number].value for residue_number in local_x]
        incomplete_value = int(row["incomplete_data_use"])
        end_corr = int(row["end_corr"])
        linestyle, linewidth, alpha = incomplete_styles[incomplete_value]
        ax.plot(
            local_x,
            local_y,
            color=end_corr_colors[end_corr],
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            label=(
                f"{row['label']} "
                f"(RMSD={float(row['rmsd']):.4f})"
            ),
            zorder=2,
        )
    ax.set_title(
        f"RCI comparison for {target_id}: server vs all local sweep traces",
        fontsize=16,
    )
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RCI")
    ax.grid(axis="y", color="#d6d3d1", linewidth=0.8)
    ax.legend(
        frameon=True,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        ncol=1,
    )
    all_overlay_path = output_dir / f"{target_id}_all_traces_overlay.png"
    fig.savefig(all_overlay_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    generated.append(all_overlay_path)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(12.0, 6.5), constrained_layout=True)
    fig.patch.set_facecolor("#f7f3eb")
    ax.set_facecolor("#fffdfa")
    ax.plot(server_x, server_y, color="#111827", linewidth=2.8, label="Server RCI", zorder=5)
    color_cycle = ["#1d4ed8", "#0f766e", "#b45309", "#be123c"]
    for index, row in enumerate(best_rows):
        local_records = parse_series_file(Path(row["rci_path"]))
        local_x = sorted(local_records)
        local_y = [local_records[residue_number].value for residue_number in local_x]
        ax.plot(
            local_x,
            local_y,
            label=f"{row['label']} (RMSD={float(row['rmsd']):.4f})",
            color=color_cycle[index % len(color_cycle)],
            linewidth=2.0,
            alpha=0.92,
        )
    ax.set_title(f"RCI comparison for {target_id}: server vs best local parameter sets", fontsize=16)
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RCI")
    ax.grid(axis="y", color="#d6d3d1", linewidth=0.8)
    ax.legend(frameon=True)
    overlay_path = output_dir / f"{target_id}_best_overlay.png"
    fig.savefig(overlay_path, dpi=180)
    plt.close(fig)
    generated.append(overlay_path)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(11.0, 5.8), constrained_layout=True)
    fig.patch.set_facecolor("#f7f3eb")
    ax.set_facecolor("#fffdfa")
    n_term_cutoff = min(server_x[:n_term_window]) if server_x else 1
    n_term_limit = server_x[min(n_term_window - 1, len(server_x) - 1)] if server_x else n_term_window
    ax.plot(server_x, server_y, color="#111827", linewidth=2.8, label="Server RCI", zorder=5)
    for index, row in enumerate(best_rows):
        local_records = parse_series_file(Path(row["rci_path"]))
        local_x = sorted(local_records)
        local_y = [local_records[residue_number].value for residue_number in local_x]
        ax.plot(
            local_x,
            local_y,
            label=row["label"],
            color=color_cycle[index % len(color_cycle)],
            linewidth=2.0,
            alpha=0.92,
        )
    ax.set_xlim(n_term_cutoff, n_term_limit)
    ax.set_title(f"N-terminal focus for {target_id} (first {n_term_window} residues)")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RCI")
    ax.grid(axis="y", color="#d6d3d1", linewidth=0.8)
    ax.legend(frameon=True, fontsize=9)
    n_term_path = output_dir / f"{target_id}_nterm_focus.png"
    fig.savefig(n_term_path, dpi=180)
    plt.close(fig)
    generated.append(n_term_path)

    plt.close("all")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    fig.patch.set_facecolor("#f7f3eb")
    metric_specs = [("rmsd", "Overall RMSD"), ("n_term_rmsd", f"N-term RMSD ({n_term_window} aa)")]
    for axis, (metric_key, title) in zip(axes, metric_specs):
        matrix = {0: {}, 1: {}}
        for row in result_rows:
            matrix[int(row["incomplete_data_use"])][int(row["end_corr"])] = float(row[metric_key])
        image = [
            [matrix[incomplete_value].get(end_corr, math.nan) for end_corr in range(7)]
            for incomplete_value in [0, 1]
        ]
        heatmap = axis.imshow(image, cmap="YlOrRd_r", aspect="auto")
        axis.set_title(title)
        axis.set_xlabel("end_corr")
        axis.set_ylabel("incomplete_data_use")
        axis.set_xticks(range(7), labels=[str(item) for item in range(7)])
        axis.set_yticks([0, 1], labels=["0", "1"])
        for row_index, incomplete_value in enumerate([0, 1]):
            for column_index, end_corr in enumerate(range(7)):
                value = matrix[incomplete_value].get(end_corr)
                if value is None:
                    continue
                axis.text(column_index, row_index, f"{value:.3f}", ha="center", va="center", fontsize=8)
        fig.colorbar(heatmap, ax=axis, shrink=0.84)
    heatmap_path = output_dir / f"{target_id}_rmsd_heatmaps.png"
    fig.savefig(heatmap_path, dpi=180)
    plt.close(fig)
    generated.append(heatmap_path)

    return generated


def infer_target_id(server_rci_path: Path, input_path: Path) -> str:
    match = re.match(r"^(\d+)", server_rci_path.name)
    if match:
        return match.group(1)
    match = re.match(r"^(\d+)", input_path.name)
    if match:
        return match.group(1)
    return "target"


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    server_rci_path = args.server_rci.resolve()
    output_dir = args.output_dir.resolve()
    python_executable = args.python_executable.resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not server_rci_path.is_file():
        raise FileNotFoundError(f"Server RCI file not found: {server_rci_path}")
    if not python_executable.is_file():
        raise FileNotFoundError(f"Python executable not found: {python_executable}")

    output_dir.mkdir(parents=True, exist_ok=True)
    mpl_config_dir = output_dir / ".mplconfig"
    mpl_config_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    target_id = infer_target_id(server_rci_path, input_path)

    audit_text, gly_zero_ha = audit_shifty_file(input_path)
    audit_path = output_dir / "shifty_input_audit.txt"
    audit_path.write_text(audit_text, encoding="utf-8")

    server_records = parse_series_file(server_rci_path)
    configs = [SweepConfig(incomplete_data_use=incomplete, end_corr=end_corr) for incomplete in (0, 1) for end_corr in range(7)]

    rows: list[dict[str, object]] = []
    for config in configs:
        run_payload = (
            reuse_local_rci(config, input_path, output_dir)
            if args.skip_runs
            else run_local_rci(config, input_path, output_dir, python_executable)
        )
        local_records = parse_series_file(Path(run_payload["rci_path"]))
        metrics = compute_metrics(server_records, local_records, args.n_term_window)
        rows.append(
            {
                "label": config.label,
                "incomplete_data_use": config.incomplete_data_use,
                "end_corr": config.end_corr,
                "rmsd": metrics["rmsd"],
                "n_term_rmsd": metrics["n_term_rmsd"],
                "max_abs_diff": metrics["max_abs_diff"],
                "max_abs_diff_residue": metrics["max_abs_diff_residue"],
                "common_residues": metrics["common_residues"],
                "server_only_count": metrics["server_only_count"],
                "local_only_count": metrics["local_only_count"],
                "resname_mismatch_count": metrics["resname_mismatch_count"],
                "server_only_preview": metrics["server_only_preview"],
                "local_only_preview": metrics["local_only_preview"],
                "rci_path": str(run_payload["rci_path"]),
                "run_dir": str(run_payload["run_dir"]),
            }
        )

    summary_csv_path = output_dir / "rci_parameter_sweep_summary.csv"
    write_summary_csv(rows, summary_csv_path)
    summary_md_path = output_dir / "README.md"
    write_summary_markdown(target_id, rows, input_path, server_rci_path, audit_path, summary_md_path)

    if not args.no_plots:
        plot_results(target_id, server_records, rows, output_dir, args.n_term_window)

    sorted_rows = sorted(rows, key=lambda row: (float(row["rmsd"]), float(row["n_term_rmsd"])))
    print(f"Wrote sweep outputs to {output_dir}")
    print(f"Input audit: {audit_path}")
    print(f"Summary CSV: {summary_csv_path}")
    print(f"Summary README: {summary_md_path}")
    print("Top 5 parameter sets by RMSD:")
    for row in sorted_rows[:5]:
        print(
            f"  {row['label']}: RMSD={float(row['rmsd']):.6f}, "
            f"N-term RMSD={float(row['n_term_rmsd']):.6f}, "
            f"max|diff|={float(row['max_abs_diff']):.6f}"
        )
    if gly_zero_ha:
        print(f"Warning: glycine HA=0.00 detected at residues {gly_zero_ha}")
    else:
        print("No glycine HA=0.00 values were found in the SHIFTY input.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
