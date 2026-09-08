#!/usr/bin/env python3
"""
Compute N/H CSPs and HSQC figures for accepted apo–apo BMRB pairs.

For each row in apo_apo_pairs.csv:
  - parse both STAR files
  - compute N/H CSPs with grid referencing (query = apo, match = second state)
  - write csp_table.csv, csp_sequence_alignment.txt, offset_grid_* artifacts
  - write hsqc_scatter.png and csp_distribution.png (apo-sequence x ticks)

Output directory per pair:
  {out_root}/{query_apo_bmrb}_{match_apo_bmrb}/
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.align import align_global  # noqa: E402
from scripts.apo_holo_exp_conditions import resolve_star_path  # noqa: E402
from scripts.bmrb_io import parse_sequence_and_shifts_from_saveframes  # noqa: E402
from scripts.config import Paths, Referencing  # noqa: E402
from scripts.csp import CSPResult, compute_csp_multiple_saveframes  # noqa: E402
from scripts.HSQC_visualize import plot_hsqc_variants  # noqa: E402
from scripts.visualize import plot_csp_distribution  # noqa: E402

# Force non-interactive matplotlib before any plotting
import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)

# Matplotlib is not thread-safe (mathtext / Agg); serialize CSP grid + plots.
_MPL_LOCK = threading.Lock()


def resolve_pair_star(
    bmrb_id: str,
    csv_path: str,
    *,
    shifts_dir: Path,
    cs_dir: Path,
) -> Optional[Path]:
    """Resolve STAR path from CSV value, shifts/, or CS_Lists/."""
    if csv_path:
        p = Path(csv_path)
        if not p.is_absolute():
            p = _REPO_ROOT / p
        if p.is_file() and p.stat().st_size > 0:
            return p
    for base in (shifts_dir, cs_dir):
        got = resolve_star_path(bmrb_id, base)
        if got is not None:
            return got
    return None


def write_sequence_alignment(
    out_path: Path,
    *,
    query_id: str,
    match_id: str,
    query_seqs: List[Tuple[Any, ...]],
    match_seqs: List[Tuple[Any, ...]],
) -> Optional[float]:
    """Write csp_sequence_alignment.txt; return best alignment score."""
    best_score = float("-inf")
    best_aligned_apo = None
    best_aligned_holo = None
    best_mapping = None

    for apo_seq, *_rest in query_seqs:
        for holo_seq, *_rest2 in match_seqs:
            aligned_apo, aligned_holo, mapping, score = align_global(apo_seq, holo_seq)
            if score > best_score:
                best_score = score
                best_aligned_apo = aligned_apo
                best_aligned_holo = aligned_holo
                best_mapping = mapping

    if not best_aligned_apo or not best_aligned_holo:
        return None

    matches = sum(
        1 for a, b in zip(best_aligned_apo, best_aligned_holo) if a == b and a != "-"
    )
    gaps_apo = sum(1 for c in best_aligned_apo if c == "-")
    gaps_holo = sum(1 for c in best_aligned_holo if c == "-")
    mismatches = len(best_aligned_apo) - matches - gaps_apo - gaps_holo
    total_aligned = len(best_aligned_apo) - gaps_apo - gaps_holo
    identity = (matches / total_aligned * 100) if total_aligned > 0 else 0.0
    match_string = "".join(
        "|" if a == b and a != "-" else " "
        for a, b in zip(best_aligned_apo, best_aligned_holo)
    )

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Sequence Alignment: Apo vs Apo (match) Chemical Shift Lists\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Query apo BMRB ID: {query_id}\n")
        f.write(f"Match apo BMRB ID: {match_id}\n")
        f.write(f"Alignment Score: {best_score:.2f}\n")
        f.write(f"Number of Mapped Pairs: {len(best_mapping) if best_mapping else 0}\n\n")
        f.write("Alignment Statistics:\n")
        f.write(f"  Matches:    {matches:4d}\n")
        f.write(f"  Mismatches: {mismatches:4d}\n")
        f.write(f"  Gaps (query): {gaps_apo:4d}\n")
        f.write(f"  Gaps (match): {gaps_holo:4d}\n")
        f.write(f"  Identity:   {identity:.1f}%\n\n")
        f.write("=" * 80 + "\n\n")

        line_width = 80
        for i in range(0, len(best_aligned_apo), line_width):
            chunk_apo = best_aligned_apo[i : i + line_width]
            chunk_holo = best_aligned_holo[i : i + line_width]
            chunk_match = match_string[i : i + line_width]
            start_pos = i + 1
            end_pos = min(i + line_width, len(best_aligned_apo))
            f.write(f"Query {start_pos:4d}-{end_pos:4d}: {chunk_apo}\n")
            f.write(f"      {' ' * 11} {chunk_match}\n")
            f.write(f"Match {start_pos:4d}-{end_pos:4d}: {chunk_holo}\n")
            f.write("\n")

    return best_score


def write_csp_table(
    out_path: Path,
    results: Sequence[CSPResult],
    *,
    query_id: str,
    match_id: str,
) -> None:
    """Write slim N/H csp_table.csv (no SASA/occlusion columns)."""
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "query_apo_bmrb",
                "match_apo_bmrb",
                "apo_resi",
                "apo_aa",
                "match_resi",
                "match_aa",
                "H_apo",
                "N_apo",
                "H_match",
                "N_match",
                "H_match_original",
                "N_match_original",
                "H_offset",
                "N_offset",
                "dH",
                "dN",
                "csp_A",
                "csp_z",
                "significant",
                "significant_1sd",
                "significant_2sd",
            ]
        )
        for r in results:
            w.writerow(
                [
                    query_id,
                    match_id,
                    r.apo_index,
                    r.apo_aa,
                    r.holo_index,
                    r.holo_aa,
                    f"{r.H_apo:.4f}" if r.H_apo is not None else "",
                    f"{r.N_apo:.4f}" if r.N_apo is not None else "",
                    f"{r.H_holo:.4f}" if r.H_holo is not None else "",
                    f"{r.N_holo:.4f}" if r.N_holo is not None else "",
                    f"{r.H_holo_original:.4f}" if r.H_holo_original is not None else "",
                    f"{r.N_holo_original:.4f}" if r.N_holo_original is not None else "",
                    f"{r.H_offset:.4f}" if r.H_offset is not None else "",
                    f"{r.N_offset:.4f}" if r.N_offset is not None else "",
                    f"{r.dH:.4f}" if r.dH is not None else "",
                    f"{r.dN:.4f}" if r.dN is not None else "",
                    f"{r.csp_A:.4f}" if r.csp_A is not None else "",
                    f"{r.z_score:.4f}" if r.z_score is not None else "",
                    int(bool(r.significant)) if r.significant is not None else "",
                    int(bool(r.significant_1sd)) if r.significant_1sd is not None else "",
                    int(bool(r.significant_2sd)) if r.significant_2sd is not None else "",
                ]
            )


def process_pair(
    row: Dict[str, str],
    *,
    out_root: Path,
    shifts_dir: Path,
    cs_dir: Path,
    ref_method: str,
) -> Dict[str, Any]:
    query_id = (row.get("query_apo_bmrb") or "").strip()
    match_id = (row.get("match_apo_bmrb") or "").strip()
    pair_label = f"{query_id}_{match_id}"
    result: Dict[str, Any] = {
        "pair": pair_label,
        "ok": False,
        "error": "",
        "n_csp": 0,
    }
    if not query_id or not match_id:
        result["error"] = "missing_ids"
        return result

    pair_dir = out_root / pair_label
    pair_dir.mkdir(parents=True, exist_ok=True)

    query_star = resolve_pair_star(
        query_id,
        row.get("query_star_path") or "",
        shifts_dir=shifts_dir,
        cs_dir=cs_dir,
    )
    match_star = resolve_pair_star(
        match_id,
        row.get("match_star_path") or "",
        shifts_dir=shifts_dir,
        cs_dir=cs_dir,
    )
    if query_star is None:
        result["error"] = "query_star_missing"
        return result
    if match_star is None:
        result["error"] = "match_star_missing"
        return result

    try:
        query_seqs = parse_sequence_and_shifts_from_saveframes(str(query_star))
        match_seqs = parse_sequence_and_shifts_from_saveframes(str(match_star))
    except Exception as exc:
        result["error"] = f"parse_failed:{exc}"
        return result

    if not query_seqs:
        result["error"] = "query_no_hn_sequences"
        return result
    if not match_seqs:
        result["error"] = "match_no_hn_sequences"
        return result

    write_sequence_alignment(
        pair_dir / "csp_sequence_alignment.txt",
        query_id=query_id,
        match_id=match_id,
        query_seqs=query_seqs,
        match_seqs=match_seqs,
    )

    try:
        # Grid heatmap + HSQC/distribution plots all use matplotlib.
        with _MPL_LOCK:
            results = compute_csp_multiple_saveframes(
                query_seqs,
                match_seqs,
                query_id,
                match_id,
                "",  # no holo PDB
                None,
                referencing_method=ref_method,
                grid_params=None,
                target_id=pair_label,
                # CSP joins output_root/target_id → pair_dir for offset_grid_*
                output_root=str(out_root),
            )
            if not results:
                result["error"] = "no_csp_results"
                return result

            write_csp_table(
                pair_dir / "csp_table.csv",
                results,
                query_id=query_id,
                match_id=match_id,
            )

            plot_hsqc_variants(
                results,
                str(pair_dir / "hsqc_scatter.png"),
                title=f"{query_id} vs {match_id} HSQC comparison",
                apo_label=query_id,
                holo_label=match_id,
                binding_results=None,
            )
            plot_csp_distribution(
                results,
                str(pair_dir / "csp_distribution.png"),
                title=f"{query_id} vs {match_id} CSP along query apo sequence",
            )
    except Exception as exc:
        result["error"] = f"csp_or_plot_failed:{exc}"
        return result

    n_valid = sum(1 for r in results if r.csp_A is not None)
    result["n_csp"] = n_valid
    result["ok"] = True
    return result


def load_pairs(
    pairs_csv: Path,
    id_filter: Optional[set],
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(pairs_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            q = (row.get("query_apo_bmrb") or "").strip()
            m = (row.get("match_apo_bmrb") or "").strip()
            label = f"{q}_{m}"
            if id_filter and label not in id_filter and q not in id_filter and m not in id_filter:
                continue
            rows.append(row)
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    cfg = Paths()
    default_out = Path(cfg.outputs_dir) / "apo_apo_matches"
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--pairs",
        type=Path,
        default=default_out / "apo_apo_pairs.csv",
        help="Accepted apo–apo pairs CSV",
    )
    ap.add_argument(
        "--out-root",
        type=Path,
        default=default_out,
        help="Root for per-pair subdirs (default: outputs/apo_apo_matches)",
    )
    ap.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(cfg.cs_cache_dir),
        help="STAR cache (CS_Lists)",
    )
    ap.add_argument(
        "--shifts-dir",
        type=Path,
        default=None,
        help="Preferred STAR dir (default: {out-root}/shifts)",
    )
    ap.add_argument(
        "--ids",
        type=str,
        default="",
        help="Comma-separated pair labels (query_match) or BMRB IDs to filter",
    )
    ap.add_argument("--workers", type=int, default=1, help="Parallel pair workers")
    args = ap.parse_args(list(argv) if argv is not None else None)

    id_filter: Optional[set] = None
    if args.ids.strip():
        id_filter = {x.strip() for x in args.ids.split(",") if x.strip()}

    pairs = load_pairs(args.pairs, id_filter)
    if not pairs:
        print("No pairs to process.", file=sys.stderr)
        return 1

    shifts_dir = args.shifts_dir or (args.out_root / "shifts")
    ref_method = Referencing().method
    args.out_root.mkdir(parents=True, exist_ok=True)

    print(f"Processing {len(pairs)} apo–apo pairs (referencing={ref_method})")

    summaries: List[Dict[str, Any]] = []

    def _run(row: Dict[str, str]) -> Dict[str, Any]:
        return process_pair(
            row,
            out_root=args.out_root,
            shifts_dir=shifts_dir,
            cs_dir=args.cs_dir,
            ref_method=ref_method,
        )

    if args.workers <= 1:
        for row in pairs:
            q = (row.get("query_apo_bmrb") or "").strip()
            m = (row.get("match_apo_bmrb") or "").strip()
            print(f"[PAIR] {q}_{m} ...", flush=True)
            summary = _run(row)
            summaries.append(summary)
            status = "ok" if summary["ok"] else f"FAIL ({summary['error']})"
            print(f"  → {status} n_csp={summary['n_csp']}", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_run, row): row for row in pairs}
            for fut in as_completed(futs):
                row = futs[fut]
                q = (row.get("query_apo_bmrb") or "").strip()
                m = (row.get("match_apo_bmrb") or "").strip()
                try:
                    summary = fut.result()
                except Exception as exc:
                    summary = {
                        "pair": f"{q}_{m}",
                        "ok": False,
                        "error": f"exception:{exc}",
                        "n_csp": 0,
                    }
                summaries.append(summary)
                status = "ok" if summary["ok"] else f"FAIL ({summary['error']})"
                print(f"[PAIR] {summary['pair']} → {status} n_csp={summary['n_csp']}", flush=True)

    n_ok = sum(1 for s in summaries if s["ok"])
    n_fail = len(summaries) - n_ok
    print(f"Done: {n_ok} ok, {n_fail} failed out of {len(summaries)}")

    # Write run summary next to pairs CSV
    summary_path = args.out_root / "apo_apo_csp_run_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pair", "ok", "n_csp", "error"])
        w.writeheader()
        for s in sorted(summaries, key=lambda x: x["pair"]):
            w.writerow(
                {
                    "pair": s["pair"],
                    "ok": s["ok"],
                    "n_csp": s["n_csp"],
                    "error": s["error"],
                }
            )
    print(f"Wrote {summary_path}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
