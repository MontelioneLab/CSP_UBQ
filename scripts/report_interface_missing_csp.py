#!/usr/bin/env python3
"""
Report missing holo HN shifts / missing CSPs at the binding-site interface.

Interface = union used by merge_csv.compute_classification / PyMOL viz:
  occluded OR CA-distance OR H-bond/charge/π OR any-atom < 2 Å.

Large-|ΔN| exclusions are reported separately and are not counted as missing CSP.
Primary fraction: n_missing_holo_hn / n_interface (from master_alignment.csv).

Examples:
  python scripts/report_interface_missing_csp.py --outputs-dir outputs
  python scripts/report_interface_missing_csp.py --targets-csv data/CSP_UBQ.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from .csp import _exceeds_max_abs_delta_n, _max_abs_delta_n_ppm
    from .target_resolution import load_and_resolve
except ImportError:
    from csp import _exceeds_max_abs_delta_n, _max_abs_delta_n_ppm
    from target_resolution import load_and_resolve

INTERFACE_FLAG_COLUMNS: Tuple[str, ...] = (
    "is_occluded_occlusion",
    "passes_filter_distance",
    "has_hbond_interaction",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "passes_sub_2A_filter_any_atom",
)

STRUCTURAL_INTERFACE_SPECS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("occlusion_analysis.csv", ("is_occluded",)),
    ("ca_distance_filter.csv", ("passes_filter",)),
    (
        "interaction_filter.csv",
        ("has_hbond", "has_charge_complement", "has_pi_contact"),
    ),
    ("any_atom_distance_filter.csv", ("passes_sub_2A_filter",)),
)

MANIFEST_FIELDS: Tuple[str, ...] = (
    "target_dir",
    "apo_bmrb",
    "holo_bmrb",
    "holo_pdb",
    "n_interface",
    "n_missing_holo_hn",
    "fraction_missing_holo_hn",
    "n_missing_csp",
    "fraction_missing_csp",
    "n_excluded_large_dn_interface",
    "n_interface_not_in_alignment",
)

DETAIL_FIELDS: Tuple[str, ...] = (
    "target_dir",
    "apo_bmrb",
    "holo_bmrb",
    "holo_pdb",
    "pdb_residue_number",
    "apo_resi",
    "apo_aa",
    "holo_resi",
    "holo_aa",
    "H_apo",
    "N_apo",
    "H_holo",
    "N_holo",
    "N_holo_original",
    "csp_A",
    "missing_holo_hn",
    "missing_csp",
    "excluded_large_dn",
    "is_occluded_occlusion",
    "passes_filter_distance",
    "has_hbond_interaction",
    "has_charge_complement_interaction",
    "has_pi_contact_interaction",
    "passes_sub_2A_filter_any_atom",
)


def _to_clean_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    return str(value).strip()


def to_bool(value: Any) -> bool:
    """Parse common bool-like CSV values (merge_csv / analyze_targets style)."""
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = _to_clean_str(value).lower()
    if not text:
        return False
    return text in ("true", "1", "yes", "y", "t")


def is_empty(value: Any) -> bool:
    return _to_clean_str(value) == ""


def _parse_float(value: Any) -> Optional[float]:
    text = _to_clean_str(value)
    if not text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def is_interface_row(row: Dict[str, str]) -> bool:
    return any(to_bool(row.get(col, "")) for col in INTERFACE_FLAG_COLUMNS)


def is_large_dn_row(row: Dict[str, str], *, threshold_ppm: Optional[float] = None) -> bool:
    n_apo = _parse_float(row.get("N_apo"))
    n_holo_raw = _parse_float(row.get("N_holo_original"))
    if n_holo_raw is None:
        n_holo_raw = _parse_float(row.get("N_holo"))
    return _exceeds_max_abs_delta_n(n_apo, n_holo_raw, threshold_ppm=threshold_ppm)


def is_missing_holo_hn(row: Dict[str, str]) -> bool:
    return is_empty(row.get("H_holo")) or is_empty(row.get("N_holo"))


def is_missing_csp(row: Dict[str, str], *, large_dn: bool) -> bool:
    if large_dn:
        return False
    return is_empty(row.get("csp_A"))


def _res_key(row: Dict[str, str], keys: Sequence[str]) -> str:
    for key in keys:
        text = _to_clean_str(row.get(key, ""))
        if text:
            return text
    return ""


def structural_interface_residue_numbers(target_dir: Path) -> Set[str]:
    """PDB residue numbers in the structural interface union (component CSVs)."""
    residues: Set[str] = set()
    for filename, flag_cols in STRUCTURAL_INTERFACE_SPECS:
        path = target_dir / filename
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if any(to_bool(row.get(col, "")) for col in flag_cols):
                    key = _res_key(
                        row,
                        ("residue_number", "resi", "pdb_residue_number", "resseq"),
                    )
                    if key:
                        residues.add(key)
    return residues


def discover_target_dirs(outputs_dir: Path) -> List[Path]:
    if not outputs_dir.is_dir():
        return []
    dirs = [
        p
        for p in sorted(outputs_dir.iterdir())
        if p.is_dir() and (p / "master_alignment.csv").is_file()
    ]
    return dirs


@dataclass
class TargetReport:
    target_dir: str
    apo_bmrb: str = ""
    holo_bmrb: str = ""
    holo_pdb: str = ""
    n_interface: int = 0
    n_missing_holo_hn: int = 0
    n_missing_csp: int = 0
    n_excluded_large_dn_interface: int = 0
    n_interface_not_in_alignment: int = 0
    detail_rows: List[Dict[str, str]] = field(default_factory=list)

    @property
    def fraction_missing_holo_hn(self) -> float:
        if self.n_interface <= 0:
            return 0.0
        return self.n_missing_holo_hn / self.n_interface

    @property
    def fraction_missing_csp(self) -> float:
        if self.n_interface <= 0:
            return 0.0
        return self.n_missing_csp / self.n_interface

    def manifest_row(self) -> Dict[str, str]:
        return {
            "target_dir": self.target_dir,
            "apo_bmrb": self.apo_bmrb,
            "holo_bmrb": self.holo_bmrb,
            "holo_pdb": self.holo_pdb,
            "n_interface": str(self.n_interface),
            "n_missing_holo_hn": str(self.n_missing_holo_hn),
            "fraction_missing_holo_hn": f"{self.fraction_missing_holo_hn:.6f}",
            "n_missing_csp": str(self.n_missing_csp),
            "fraction_missing_csp": f"{self.fraction_missing_csp:.6f}",
            "n_excluded_large_dn_interface": str(self.n_excluded_large_dn_interface),
            "n_interface_not_in_alignment": str(self.n_interface_not_in_alignment),
        }


def analyze_target(
    target_dir: Path,
    *,
    threshold_ppm: Optional[float] = None,
) -> TargetReport:
    alignment_path = target_dir / "master_alignment.csv"
    report = TargetReport(target_dir=target_dir.name)

    with alignment_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if rows:
        first = rows[0]
        report.apo_bmrb = _to_clean_str(first.get("apo_bmrb"))
        report.holo_bmrb = _to_clean_str(first.get("holo_bmrb"))
        report.holo_pdb = _to_clean_str(first.get("holo_pdb"))

    aligned_pdb_resis: Set[str] = set()
    for row in rows:
        key = _to_clean_str(row.get("pdb_residue_number"))
        if key:
            aligned_pdb_resis.add(key)

        if not is_interface_row(row):
            continue

        report.n_interface += 1
        large_dn = is_large_dn_row(row, threshold_ppm=threshold_ppm)
        if large_dn:
            report.n_excluded_large_dn_interface += 1

        missing_holo = is_missing_holo_hn(row)
        missing_csp = is_missing_csp(row, large_dn=large_dn)

        if missing_holo:
            report.n_missing_holo_hn += 1
        if missing_csp:
            report.n_missing_csp += 1

        if missing_holo or missing_csp:
            detail = {
                "target_dir": target_dir.name,
                "apo_bmrb": report.apo_bmrb,
                "holo_bmrb": report.holo_bmrb,
                "holo_pdb": report.holo_pdb,
                "pdb_residue_number": _to_clean_str(row.get("pdb_residue_number")),
                "apo_resi": _to_clean_str(row.get("apo_resi")),
                "apo_aa": _to_clean_str(row.get("apo_aa")),
                "holo_resi": _to_clean_str(row.get("holo_resi")),
                "holo_aa": _to_clean_str(row.get("holo_aa")),
                "H_apo": _to_clean_str(row.get("H_apo")),
                "N_apo": _to_clean_str(row.get("N_apo")),
                "H_holo": _to_clean_str(row.get("H_holo")),
                "N_holo": _to_clean_str(row.get("N_holo")),
                "N_holo_original": _to_clean_str(row.get("N_holo_original")),
                "csp_A": _to_clean_str(row.get("csp_A")),
                "missing_holo_hn": "1" if missing_holo else "0",
                "missing_csp": "1" if missing_csp else "0",
                "excluded_large_dn": "1" if large_dn else "0",
            }
            for col in INTERFACE_FLAG_COLUMNS:
                detail[col] = "1" if to_bool(row.get(col, "")) else "0"
            report.detail_rows.append(detail)

    structural = structural_interface_residue_numbers(target_dir)
    if structural:
        report.n_interface_not_in_alignment = len(structural - aligned_pdb_resis)

    return report


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Report fraction of missing holo HN shifts and missing CSPs at the "
            "binding-site interface (large-|ΔN| exclusions ignored for missing CSP)."
        )
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing per-target subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=None,
        help="Optional CSV with holo_pdb (+ apo_bmrb) to select targets.",
    )
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=None,
        help=(
            "Directory for the dataset manifest "
            "(default: <outputs-dir>/summary_statistics)."
        ),
    )
    parser.add_argument(
        "--write-per-target",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write per-target interface_missing_csp.csv (default: on).",
    )
    parser.add_argument(
        "--threshold-ppm",
        type=float,
        default=None,
        help=(
            "Large-|ΔN| cutoff in ppm "
            f"(default: config Compute.max_abs_delta_n_ppm = {_max_abs_delta_n_ppm():.1f})."
        ),
    )
    return parser.parse_args(argv)


def resolve_targets(outputs_dir: Path, targets_csv: Optional[Path]) -> List[Path]:
    if targets_csv is not None:
        resolved = load_and_resolve(targets_csv, outputs_dir)
        out: List[Path] = []
        for path in resolved:
            if (path / "master_alignment.csv").is_file():
                out.append(path)
            else:
                print(
                    f"[WARN] Skipping {path.name}: missing master_alignment.csv",
                    file=sys.stderr,
                )
        return out
    return discover_target_dirs(outputs_dir)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    outputs_dir = args.outputs_dir
    summary_dir = args.summary_dir or (outputs_dir / "summary_statistics")

    if not outputs_dir.is_dir():
        print(f"[ERROR] Outputs directory not found: {outputs_dir}", file=sys.stderr)
        return 1

    target_dirs = resolve_targets(outputs_dir, args.targets_csv)
    if not target_dirs:
        print("[WARN] No targets with master_alignment.csv found.", file=sys.stderr)
        write_csv(summary_dir / "interface_missing_csp_per_target.csv", MANIFEST_FIELDS, [])
        return 0

    reports: List[TargetReport] = []
    for target_dir in target_dirs:
        try:
            report = analyze_target(target_dir, threshold_ppm=args.threshold_ppm)
        except Exception as exc:  # noqa: BLE001 — keep scanning other targets
            print(f"[WARN] Failed on {target_dir.name}: {exc}", file=sys.stderr)
            continue
        reports.append(report)
        if args.write_per_target:
            write_csv(
                target_dir / "interface_missing_csp.csv",
                DETAIL_FIELDS,
                report.detail_rows,
            )

    manifest_path = summary_dir / "interface_missing_csp_per_target.csv"
    write_csv(manifest_path, MANIFEST_FIELDS, (r.manifest_row() for r in reports))

    n_targets = len(reports)
    sum_iface = sum(r.n_interface for r in reports)
    sum_holo = sum(r.n_missing_holo_hn for r in reports)
    sum_csp = sum(r.n_missing_csp for r in reports)
    sum_ldn = sum(r.n_excluded_large_dn_interface for r in reports)
    with_iface = [r for r in reports if r.n_interface > 0]
    mean_frac_holo = (
        sum(r.fraction_missing_holo_hn for r in with_iface) / len(with_iface)
        if with_iface
        else 0.0
    )
    mean_frac_csp = (
        sum(r.fraction_missing_csp for r in with_iface) / len(with_iface)
        if with_iface
        else 0.0
    )
    overall_frac_holo = (sum_holo / sum_iface) if sum_iface else 0.0
    overall_frac_csp = (sum_csp / sum_iface) if sum_iface else 0.0

    print(f"Targets analyzed: {n_targets}")
    print(f"Interface residues (sum): {sum_iface}")
    print(
        f"Missing holo HN at interface: {sum_holo} "
        f"(overall fraction {overall_frac_holo:.4f}; "
        f"mean per-target {mean_frac_holo:.4f})"
    )
    print(
        f"Missing CSP at interface (excl. large-|ΔN|): {sum_csp} "
        f"(overall fraction {overall_frac_csp:.4f}; "
        f"mean per-target {mean_frac_csp:.4f})"
    )
    print(f"Large-|ΔN| at interface (ignored for missing CSP): {sum_ldn}")
    print(f"Manifest: {manifest_path}")
    if args.write_per_target:
        print("Per-target detail: <outputs>/<TARGET>/interface_missing_csp.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
