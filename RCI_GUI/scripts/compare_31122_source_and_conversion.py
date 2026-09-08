#!/usr/bin/env python3
"""Compare 31122 source-file conversion and downstream local RCI outputs."""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from format_conversion import convert_file, list_input_entities


ROOT = Path(__file__).resolve().parents[1]
RCI_SCRIPT = ROOT / "scripts" / "rci_v_1c_PyNMR-STAR.py"
SERVER_RESULTS_DIR = ROOT / "server_rci_results"
DEFAULT_BMRB_INPUT = SERVER_RESULTS_DIR / "31122_3.str.txt"
DEFAULT_SERVER_SHIFTY = SERVER_RESULTS_DIR / "31122_3.str.shifty.txt"
DEFAULT_SERVER_RCI = SERVER_RESULTS_DIR / "31122_3.str.shifty.txt_BMRB.RCI.txt"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "comparisons" / "31122_workshop" / "source_compare"
SHIFTY_COLUMNS = ("HA", "CA", "CB", "CO", "N", "HN")


@dataclass(frozen=True)
class SeriesRecord:
    residue_number: int
    residue_name: str
    value: float


@dataclass(frozen=True)
class RunSpec:
    label: str
    input_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the 31122 BMRB source file to SHIFTY, compare that conversion "
            "against the server-reported SHIFTY file, and compare local RCI outputs "
            "from the source and SHIFTY inputs against the server RCI trace."
        )
    )
    parser.add_argument("--bmrb-input", type=Path, default=DEFAULT_BMRB_INPUT)
    parser.add_argument("--server-shifty", type=Path, default=DEFAULT_SERVER_SHIFTY)
    parser.add_argument("--server-rci", type=Path, default=DEFAULT_SERVER_RCI)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--python-executable", type=Path, default=Path(sys.executable))
    parser.add_argument("--entity-id", default=None, help="Optional entity id for NMR-STAR 3 input.")
    parser.add_argument("--value-tolerance", type=float, default=1e-6)
    parser.add_argument("--skip-runs", action="store_true")
    return parser.parse_args()


def ensure_clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


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


def read_shifty_table(path: Path) -> tuple[list[str], dict[int, dict[str, str]]]:
    header: list[str] | None = None
    rows: dict[int, dict[str, str]] = {}
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
        row = dict(zip(header, parts))
        rows[int(row["#NUM"])] = row
    if header is None:
        raise ValueError(f"No SHIFTY header found in {path}")
    return header, rows


def compare_shifty_tables(
    local_path: Path,
    server_path: Path,
    output_dir: Path,
    tolerance: float,
) -> tuple[Path, Path]:
    _header_local, local_rows = read_shifty_table(local_path)
    _header_server, server_rows = read_shifty_table(server_path)

    summary_lines = [
        f"Local converted SHIFTY: {local_path}",
        f"Server SHIFTY: {server_path}",
        "",
    ]
    all_residues = sorted(set(local_rows) | set(server_rows))
    diff_rows: list[dict[str, object]] = []

    missing_local = [residue for residue in all_residues if residue not in local_rows]
    missing_server = [residue for residue in all_residues if residue not in server_rows]
    summary_lines.append(f"Residues only in local conversion: {missing_server if missing_server else 'none'}")
    summary_lines.append(f"Residues only in server SHIFTY: {missing_local if missing_local else 'none'}")
    summary_lines.append("")

    per_column_counts = {column: 0 for column in ("AA", *SHIFTY_COLUMNS)}
    largest_diffs: dict[str, tuple[int, float] | None] = {column: None for column in SHIFTY_COLUMNS}

    for residue_number in sorted(set(local_rows) & set(server_rows)):
        local_row = local_rows[residue_number]
        server_row = server_rows[residue_number]
        if local_row["AA"] != server_row["AA"]:
            per_column_counts["AA"] += 1
            diff_rows.append(
                {
                    "residue_number": residue_number,
                    "column": "AA",
                    "local_value": local_row["AA"],
                    "server_value": server_row["AA"],
                    "difference": "",
                }
            )
        for column in SHIFTY_COLUMNS:
            local_value = float(local_row[column])
            server_value = float(server_row[column])
            difference = local_value - server_value
            if abs(difference) <= tolerance:
                continue
            per_column_counts[column] += 1
            if largest_diffs[column] is None or abs(difference) > abs(largest_diffs[column][1]):
                largest_diffs[column] = (residue_number, difference)
            diff_rows.append(
                {
                    "residue_number": residue_number,
                    "column": column,
                    "local_value": local_row[column],
                    "server_value": server_row[column],
                    "difference": difference,
                }
            )

    summary_lines.append("Per-column difference counts:")
    for column in ("AA", *SHIFTY_COLUMNS):
        summary_lines.append(f"  {column}: {per_column_counts[column]}")
    summary_lines.append("")
    summary_lines.append("Largest numeric deviations by column:")
    for column in SHIFTY_COLUMNS:
        largest = largest_diffs[column]
        if largest is None:
            summary_lines.append(f"  {column}: none")
        else:
            residue_number, difference = largest
            summary_lines.append(f"  {column}: residue {residue_number}, diff={difference:+.6f}")
    summary_lines.append("")

    diff_csv = output_dir / "shifty_conversion_differences.csv"
    with diff_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["residue_number", "column", "local_value", "server_value", "difference"],
        )
        writer.writeheader()
        writer.writerows(diff_rows)

    summary_path = output_dir / "shifty_conversion_summary.txt"
    summary_lines.append(f"Detailed per-field differences: {diff_csv}")
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    return summary_path, diff_csv


def compute_rmsd(reference: dict[int, SeriesRecord], probe: dict[int, SeriesRecord]) -> tuple[float | None, int]:
    common_residues = sorted(set(reference) & set(probe))
    squared_diffs: list[float] = []
    for residue_number in common_residues:
        diff = reference[residue_number].value - probe[residue_number].value
        squared_diffs.append(diff * diff)
    if not squared_diffs:
        return None, 0
    return math.sqrt(sum(squared_diffs) / len(squared_diffs)), len(common_residues)


def run_local_rci(
    spec: RunSpec,
    output_dir: Path,
    python_executable: Path,
    entity_id: str | None,
) -> Path:
    run_dir = output_dir / "rci_runs" / spec.label
    ensure_clean_dir(run_dir)
    local_input = run_dir / spec.input_path.name
    shutil.copy2(spec.input_path, local_input)
    mpl_config_dir = run_dir / ".mplconfig"
    mpl_config_dir.mkdir(exist_ok=True)

    command = [str(python_executable), str(RCI_SCRIPT), "-b", spec.input_path.name, "-mpl"]
    if entity_id:
        command.extend(["-entity", entity_id])
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
    if process.returncode != 0:
        raise RuntimeError(f"{spec.label}: RCI run failed; see {run_dir / 'stderr.log'}")
    rci_path = run_dir / f"{spec.input_path.name}.RCI.txt"
    if not rci_path.is_file():
        raise FileNotFoundError(f"{spec.label}: missing RCI output at {rci_path}")
    return rci_path


def reuse_local_rci(spec: RunSpec, output_dir: Path) -> Path:
    rci_path = output_dir / "rci_runs" / spec.label / f"{spec.input_path.name}.RCI.txt"
    if not rci_path.is_file():
        raise FileNotFoundError(f"{spec.label}: reusable output not found at {rci_path}")
    return rci_path


def plot_rci_overlay(
    output_path: Path,
    server_records: dict[int, SeriesRecord],
    run_records: dict[str, dict[int, SeriesRecord]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    server_x = sorted(server_records)
    server_y = [server_records[residue_number].value for residue_number in server_x]

    plt.close("all")
    fig, ax = plt.subplots(figsize=(11.5, 6.2), constrained_layout=True)
    fig.patch.set_facecolor("#f7f3eb")
    ax.set_facecolor("#fffdfa")
    ax.plot(server_x, server_y, color="#111827", linewidth=2.9, label="Server RCI reference")

    styles = {
        "bmrb_source": ("#1d4ed8", "Local RCI from 31122_3.str.txt"),
        "local_converted_shifty": ("#0f766e", "Local RCI from locally converted SHIFTY"),
        "server_shifty_input": ("#b45309", "Local RCI from server SHIFTY input"),
    }
    for label, records in run_records.items():
        x_values = sorted(records)
        y_values = [records[residue_number].value for residue_number in x_values]
        color, title = styles[label]
        ax.plot(x_values, y_values, color=color, linewidth=2.0, alpha=0.94, label=title)

    ax.set_title("31122 local RCI traces from different source files")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RCI")
    ax.grid(axis="y", color="#d6d3d1", linewidth=0.8)
    ax.legend(frameon=True, fontsize=9)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_rci_summary(
    output_path: Path,
    rows: list[dict[str, object]],
    rci_overlay_path: Path,
    shifty_summary_path: Path,
    diff_csv_path: Path,
) -> None:
    lines = [
        "# 31122 source and conversion comparison",
        "",
        f"- SHIFTY comparison summary: `{shifty_summary_path}`",
        f"- SHIFTY per-field diff CSV: `{diff_csv_path}`",
        f"- RCI overlay plot: `{rci_overlay_path}`",
        "",
        "| Source | RMSD to server RCI | Common residues | RCI file |",
        "|---|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['label']} | {float(row['rmsd']):.6f} | {int(row['common_residues'])} | `{row['rci_path']}` |"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    bmrb_input = args.bmrb_input.resolve()
    server_shifty = args.server_shifty.resolve()
    server_rci = args.server_rci.resolve()
    output_dir = args.output_dir.resolve()
    python_executable = args.python_executable.resolve()

    if not bmrb_input.is_file():
        raise FileNotFoundError(f"BMRB input not found: {bmrb_input}")
    if not server_shifty.is_file():
        raise FileNotFoundError(f"Server SHIFTY not found: {server_shifty}")
    if not server_rci.is_file():
        raise FileNotFoundError(f"Server RCI not found: {server_rci}")

    output_dir.mkdir(parents=True, exist_ok=True)
    mpl_config_dir = output_dir / ".mplconfig"
    mpl_config_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    entity_id = args.entity_id
    if entity_id is None:
        entity_infos = list_input_entities(bmrb_input)
        if len(entity_infos) == 1:
            entity_id = entity_infos[0].entity_id

    conversion_dir = output_dir / "conversion"
    conversion_dir.mkdir(parents=True, exist_ok=True)
    local_shifty = conversion_dir / "31122_3.local.shifty.txt"
    convert_file(bmrb_input, "SHIFTY", local_shifty, entity_id=entity_id)

    shifty_summary_path, diff_csv_path = compare_shifty_tables(
        local_shifty,
        server_shifty,
        output_dir,
        args.value_tolerance,
    )

    specs = [
        RunSpec("bmrb_source", bmrb_input),
        RunSpec("local_converted_shifty", local_shifty),
        RunSpec("server_shifty_input", server_shifty),
    ]
    server_records = parse_series_file(server_rci)
    summary_rows: list[dict[str, object]] = []
    run_records: dict[str, dict[int, SeriesRecord]] = {}
    for spec in specs:
        rci_path = reuse_local_rci(spec, output_dir) if args.skip_runs else run_local_rci(
            spec, output_dir, python_executable, entity_id if spec.input_path == bmrb_input else None
        )
        records = parse_series_file(rci_path)
        run_records[spec.label] = records
        rmsd, common_residues = compute_rmsd(server_records, records)
        summary_rows.append(
            {
                "label": spec.label,
                "rmsd": rmsd,
                "common_residues": common_residues,
                "rci_path": str(rci_path),
            }
        )

    overlay_path = output_dir / "31122_source_overlay.png"
    plot_rci_overlay(overlay_path, server_records, run_records)

    rci_summary_path = output_dir / "README.md"
    write_rci_summary(rci_summary_path, summary_rows, overlay_path, shifty_summary_path, diff_csv_path)

    print(f"Wrote conversion + source comparison outputs to {output_dir}")
    print(f"Local converted SHIFTY: {local_shifty}")
    print(f"SHIFTY summary: {shifty_summary_path}")
    print(f"RCI summary: {rci_summary_path}")
    for row in summary_rows:
        print(
            f"  {row['label']}: RMSD={float(row['rmsd']):.6f}, "
            f"common residues={int(row['common_residues'])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
