#!/usr/bin/env python3
"""
Split a CSP_UBQ-style table into domain-like vs full-length protein targets:

- Receptor-chain sequences come from holo/apo NMR-STAR (see bmrb_io) with PDB fallback.
- UniProt accession: prefer **RCSB polymer-entity deposition** (same entity ID as receptor chain via
  SIFS) when available — **no blastp** needed; otherwise fall back to **blastp** and same-species ranking.
- Full-length UniProt sequence is fetched from rest.uniprot.org.
- If ``uniprot_seq_length > LENGTH_RATIO * bmrb_holo_seq_length`` → ``domain``, else ``full_length``.

Incremental runs resume from ``{basename}_bifurcation_all_rows.csv`` unless ``--force-recompute``.
After each processed row all CSVs plus the summary JSON are rewritten.

"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

_repo_root = Path(__file__).resolve().parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

try:
    from .annotate_csp_csv_metadata import _load_receptor_chain_map
except ImportError:
    from scripts.annotate_csp_csv_metadata import _load_receptor_chain_map  # noqa: WPS433

try:
    from .apo_holo_exp_conditions import ensure_star_path
    from .bmrb_io import receptor_sequence_from_assigned_shifts
    from .config import paths
    from .rcsb_io import fetch_pdb, parse_pdb_sequences
    from .rcsb_data_io import fetch_holo_polymer_context, normalize_pdb_id
    from .uniprot_blast import BlastHit, best_coverage_fraction, run_blastp
    from .uniprot_io import fetch_fasta_sequence, fetch_taxonomy_id, normalize_uniprot_accession
except ImportError:
    from scripts.apo_holo_exp_conditions import ensure_star_path
    from scripts.bmrb_io import receptor_sequence_from_assigned_shifts
    from scripts.config import paths
    from scripts.rcsb_io import fetch_pdb, parse_pdb_sequences
    from scripts.rcsb_data_io import fetch_holo_polymer_context, normalize_pdb_id
    from scripts.uniprot_blast import BlastHit, best_coverage_fraction, run_blastp
    from scripts.uniprot_io import fetch_fasta_sequence, fetch_taxonomy_id, normalize_uniprot_accession


LENGTH_RATIO_DEFAULT = 1.5


def _short_err(reason: Optional[str], max_len: int = 80) -> str:
    if not reason:
        return ""
    reason = reason.replace("\n", " ")
    if len(reason) <= max_len:
        return reason
    return reason[: max_len - 3] + "..."


def row_cache_key(apo_bmrb: str, holo_bmrb: str, holo_pdb: str) -> Tuple[str, str, str]:
    """Stable tuple for incremental cache rows (aligned with receptor map join keys)."""
    return (
        apo_bmrb.strip(),
        holo_bmrb.strip(),
        (holo_pdb or "").strip().lower(),
    )


def load_bifurcation_cache(csv_path: Path) -> Tuple[Dict[Tuple[str, str, str], Dict[str, str]], List[str]]:
    """
    Load ``*_bifurcation_all_rows.csv`` into keyed rows. Last occurrence wins duplicate keys.

    Returns (cache dict, merged fieldnames from file).
    """
    if not csv_path.is_file():
        return {}, []
    merged_fields: List[str] = []
    result: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        merged_fields = list(reader.fieldnames or [])
        for row in reader:
            apo = (row.get("apo_bmrb") or "").strip()
            holo_b = (row.get("holo_bmrb") or "").strip()
            hpdb = (row.get("holo_pdb") or "").strip().lower()
            rk = row_cache_key(apo, holo_b, hpdb)
            flat = {str(k): ("" if row.get(k) is None else str(row.get(k))) for k in merged_fields}
            result[rk] = flat
    return result, merged_fields


def _merge_output_fieldnames(primary: List[str], cache_fields: List[str]) -> List[str]:
    out = list(primary)
    for c in cache_fields:
        if c and c not in out:
            out.append(c)
    return out


def _tally_enriched(rows: List[Dict[str, Any]]) -> Tuple[Dict[str, int], int]:
    """Return (count dict, domains+full_length classified)."""
    unresolved = domains = ful = 0
    for r in rows:
        if r.get("error_reason"):
            unresolved += 1
        elif (r.get("bifurcation_label") or "") == "domain":
            domains += 1
        elif (r.get("bifurcation_label") or "") == "full_length":
            ful += 1
        else:
            unresolved += 1
    return {"rows": len(rows), "domains": domains, "full_length": ful, "unresolved": unresolved}, domains + ful


def write_bifurcation_bundle(
    enriched: List[Dict[str, Any]],
    *,
    base_fieldnames: List[str],
    basename: str,
    out_dir: Path,
    started_at_perf: float,
    skipped_cached: int,
) -> Dict[str, Any]:
    """Rewrite all bifurcation outputs from the current ``enriched`` list."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_domains = out_dir / f"{basename}_domains.csv"
    out_full = out_dir / f"{basename}_full_length.csv"
    out_unresolved = out_dir / f"{basename}_bifurcation_unresolved.csv"
    out_all = out_dir / f"{basename}_bifurcation_all_rows.csv"
    summary_path = out_dir / f"{basename}_bifurcation_summary.json"

    def _filtered(path: Path, predicate: Callable[[Dict[str, Any]], bool]) -> None:
        subset = [r for r in enriched if predicate(r)]
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=base_fieldnames, extrasaction="ignore")
            w.writeheader()
            for row in subset:
                w.writerow({k: row.get(k, "") for k in base_fieldnames})

    _filtered(
        out_domains,
        lambda r: ((r.get("bifurcation_label") == "domain") and not r.get("error_reason")),
    )
    _filtered(
        out_full,
        lambda r: ((r.get("bifurcation_label") == "full_length") and not r.get("error_reason")),
    )
    _filtered(out_unresolved, lambda r: bool(r.get("error_reason")))

    with open(out_all, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=base_fieldnames, extrasaction="ignore")
        wr.writeheader()
        for row in enriched:
            wr.writerow({k: row.get(k, "") for k in base_fieldnames})

    counts, classified = _tally_enriched(enriched)
    elapsed = time.perf_counter() - started_at_perf
    summary_payload = {
        "counts": dict(counts),
        "elapsed_seconds": round(elapsed, 3),
        "skipped_cached": skipped_cached,
        "outputs": {
            "domains_csv": str(out_domains.resolve()),
            "full_length_csv": str(out_full.resolve()),
            "unresolved_csv": str(out_unresolved.resolve()),
            "all_rows_csv": str(out_all.resolve()),
            "summary_json": str(summary_path.resolve()),
        },
    }
    with open(summary_path, "w", encoding="utf-8") as sf:
        json.dump(summary_payload, sf, indent=2)

    out = dict(counts)
    out["skipped_cached_run"] = skipped_cached
    out["elapsed_seconds"] = round(elapsed, 3)
    out["classified_ok"] = classified
    out["outputs"] = summary_payload["outputs"]
    return out


def _pdb_fallback_sequence(pdb_raw: Optional[str], chain: str, pdb_cache_dir: str) -> Tuple[Optional[str], bool]:
    """Return PDB chain sequence (uppercase AA) when path exists locally or can be fetched."""
    pid = normalize_pdb_id(pdb_raw or "")
    ch = (chain or "").strip()
    if not pid or not ch:
        return None, False
    try:
        ppath = fetch_pdb(pid, cache_dir=pdb_cache_dir, force=False)
        chains = parse_pdb_sequences(ppath)
        ck = None
        for cand in chains:
            if cand.upper() == ch.upper():
                ck = cand
                break
        if ck is None:
            ck = ch if ch in chains else None
        if ck is None:
            return None, False
        return chains.get(ck) or chains.get(ch), True
    except Exception:
        return None, False


def _extract_sequence_bmrb_then_pdb(
    bmrb_id: str,
    receptor_chain_id: str,
    cs_dir: Path,
    do_fetch_bmrb: bool,
    pdb_fallback_id: Optional[str],
    pdb_cache_dir: str,
) -> Tuple[Optional[str], str, bool]:
    star = ensure_star_path(bmrb_id, cs_dir, do_fetch_bmrb)
    if star is None:
        if pdb_fallback_id:
            seq, ok = _pdb_fallback_sequence(pdb_fallback_id, receptor_chain_id, pdb_cache_dir)
            if seq and ok:
                return seq, "pdb_fallback_only", True
        return None, "none", False

    seq, src = receptor_sequence_from_assigned_shifts(str(star), receptor_chain_id)
    if seq:
        return seq, src, False

    if pdb_fallback_id:
        pseq, ok = _pdb_fallback_sequence(pdb_fallback_id, receptor_chain_id, pdb_cache_dir)
        if pseq and ok:
            return pseq, "pdb_fallback_only", True
    return None, "none", False


def _pick_uniprot_accession(
    hits: List[BlastHit],
    holo_tax_id: Optional[int],
    deposited: List[str],
    query_len: int,
    *,
    tax_cache: Dict[str, Optional[int]],
    min_pident: float,
    min_query_cov: float,
) -> Tuple[Optional[BlastHit], str]:
    if not hits:
        return None, "no_blast_hits"

    dep = {normalize_uniprot_accession(u) for u in deposited}
    short_dep = {d.split("-")[0] for d in dep}

    def passes(h: BlastHit) -> bool:
        cov = best_coverage_fraction(h, query_len)
        return h.pident >= min_pident and cov >= min_query_cov

    ranked = sorted(hits, key=lambda h: (-h.bitscore, -h.pident, h.subject_accession))

    passing = [h for h in ranked if passes(h)]
    pool = passing if passing else ranked

    for h in pool:
        acc_h = normalize_uniprot_accession(h.subject_accession)
        base = acc_h.split("-")[0]
        if acc_h in dep or base in short_dep:
            return h, "matched_rcsb_reference"

    if holo_tax_id is not None:
        for h in pool:
            acc_u = normalize_uniprot_accession(h.subject_accession)
            if acc_u not in tax_cache:
                tax_cache[acc_u] = fetch_taxonomy_id(acc_u)
                time.sleep(0.03)
            tid = tax_cache.get(acc_u)
            if tid == holo_tax_id:
                return h, "taxonomy_match_holo"

    return pool[0], "top_blast_hit"


def _classify_row(
    bmrb_holo_len: int,
    uniprot_len: int,
    ratio: float,
) -> str:
    if uniprot_len > ratio * bmrb_holo_len:
        return "domain"
    return "full_length"


def process_dataset(
    *,
    input_csv: Path,
    receptor_csv: Path,
    out_dir: Path,
    basename: str,
    cs_dir: Path,
    pdb_cache_dir: str,
    fetch_bmrb: bool,
    blast_db: str,
    blast_remote: bool,
    min_pident: float,
    min_query_cov: float,
    ratio: float,
    sleep_between_uniprot: float,
    dry_run_rest: bool,
    quiet: bool = False,
    force_recompute: bool = False,
) -> Dict[str, Any]:
    t0 = time.perf_counter()
    receptor_map = _load_receptor_chain_map(receptor_csv)

    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames_in = reader.fieldnames or []
        rows = list(reader)

    n_input = len(rows)

    extra_cols = [
        "receptor_chain_id",
        "bmrb_holo_seq_source",
        "bmrb_holo_seq_from_pdb_fallback",
        "bmrb_apo_seq_source",
        "bmrb_apo_seq_from_pdb_fallback",
        "bmrb_holo_seq",
        "bmrb_apo_seq",
        "bmrb_holo_seq_length",
        "bmrb_apo_seq_length",
        "pdb_polymer_entity_id",
        "uniprot_resolution_source",
        "uniprot_accession",
        "uniprot_seq",
        "uniprot_seq_length",
        "blast_bitscore",
        "blast_identity",
        "blast_rank_reason",
        "rcsb_ncbi_taxonomy_id",
        "rcsb_reference_uniprots_json",
        "bifurcation_label",
        "length_ratio_used",
        "classification_rule",
        "error_reason",
    ]

    cached_by_key: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    cache_fieldnames_from_disk: List[str] = []

    out_all_path = out_dir / f"{basename}_bifurcation_all_rows.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    if out_all_path.exists() and not force_recompute:
        cached_by_key, cache_fieldnames_from_disk = load_bifurcation_cache(out_all_path)

    base_fieldnames_in = list(fieldnames_in) + [c for c in extra_cols if c not in fieldnames_in]
    base_fieldnames = _merge_output_fieldnames(base_fieldnames_in, cache_fieldnames_from_disk)

    if not quiet:
        blast_mode = f"blastp -db {blast_db}" + (" -remote" if blast_remote else " (local)")
        print("Domain / full-length bifurcation", file=sys.stdout)
        print(f"  Input table:       {input_csv.resolve()}", file=sys.stdout)
        print(f"  Receptor chains:   {receptor_csv.resolve()} ({len(receptor_map)} mappings)", file=sys.stdout)
        print(f"  Rows (input):      {n_input}", file=sys.stdout)
        print(f"  CS cache:          {cs_dir.resolve()}", file=sys.stdout)
        print(f"  PDB cache:         {pdb_cache_dir}", file=sys.stdout)
        print(f"  Output dir:        {out_dir.resolve()}", file=sys.stdout)
        print(f"  Basename:          {basename}", file=sys.stdout)
        print(f"  Classify rule:     UniProt length > {ratio} × bmrb_holo_seq_length → domain", file=sys.stdout)
        print(f"  UniProt pref:      RCSB polymer entity (SIFS) → then {blast_mode} if missing", file=sys.stdout)
        print(
            f"  Blast hit filters: min %identity {min_pident}, min query coverage {min_query_cov:.2f}",
            file=sys.stdout,
        )
        print(f"  Fetch BMRB:        {fetch_bmrb}", file=sys.stdout)
        print(f"  Force recompute:   {force_recompute}", file=sys.stdout)
        if not force_recompute and cached_by_key:
            print(f"  Resume cache rows: {len(cached_by_key)} (merge from {out_all_path.name})", file=sys.stdout)
        if dry_run_rest:
            print("  UniProt fetch:    DRY-RUN (using holo length as UniProt length)", file=sys.stdout)
        print("", file=sys.stdout)

    enriched: List[Dict[str, Any]] = []
    tax_cache: Dict[str, Optional[int]] = {}

    run_skipped_so_far = 0
    bundle_result: Optional[Dict[str, Any]] = None

    for i, row in enumerate(rows, start=1):
        apo_bmrb = (row.get("apo_bmrb") or "").strip()
        holo_bmrb = (row.get("holo_bmrb") or "").strip()
        holo_pdb = (row.get("holo_pdb") or "").strip().lower()
        apo_pdb = (row.get("apo_pdb") or "").strip()

        rk = row_cache_key(apo_bmrb, holo_bmrb, holo_pdb)
        reused = False

        if rk in cached_by_key and not force_recompute:
            out_row = copy.deepcopy(dict(cached_by_key[rk]))
            for fk in base_fieldnames:
                out_row.setdefault(fk, "")
            enriched.append(out_row)
            reused = True
            run_skipped_so_far += 1

        if reused:
            if not quiet:
                ch_disp = ((out_row.get("receptor_chain_id") or row.get("")) or "").upper() or "-"
                lbl = ((out_row.get("bifurcation_label")) or "").upper()
                errs = _short_err(out_row.get("error_reason"))
                status = (
                    f"CACHED (resume)  {errs}" if errs else f"CACHED (resume)  {lbl or '?'}"
                )
                print(
                    f"[{i}/{n_input}] ... holo_bmrb={holo_bmrb} holo_pdb={holo_pdb or '-'} chain={ch_disp} | {status}",
                    file=sys.stdout,
                )

            bundle_result = write_bifurcation_bundle(
                enriched,
                base_fieldnames=base_fieldnames,
                basename=basename,
                out_dir=out_dir,
                started_at_perf=t0,
                skipped_cached=run_skipped_so_far,
            )
            continue

        out_row = {k: row.get(k, "") for k in fieldnames_in}
        chain = receptor_map.get((apo_bmrb, holo_bmrb, holo_pdb)) if holo_pdb else None
        err: Optional[str] = None
        if not chain:
            err = "no_receptor_chain_mapping"
        out_row["receptor_chain_id"] = chain or ""

        holo_seq: Optional[str] = None
        holo_src = ""
        holo_pdb_fb = False
        apo_seq: Optional[str] = None
        apo_src = ""
        apo_pdb_fb = False
        holo_ctx = None
        acc = ""
        uni_seq = ""
        uni_len = 0
        blast_bs = ""
        blast_id = ""
        rank_reason = ""
        rcsb_tax = ""
        rcsb_refs = "[]"
        pdb_entity_id_col = ""
        uniprot_source = ""

        rule = f"uniprot_len > {ratio} * bmrb_holo_seq_length => domain"

        if not err:
            assert chain is not None
            holo_seq, holo_src, holo_pdb_fb = _extract_sequence_bmrb_then_pdb(
                holo_bmrb, chain, cs_dir, fetch_bmrb, row.get("holo_pdb") or None, pdb_cache_dir,
            )

            apo_seq, apo_src, apo_pdb_fb = _extract_sequence_bmrb_then_pdb(
                apo_bmrb, chain, cs_dir, fetch_bmrb, apo_pdb or None, pdb_cache_dir,
            )

            if not holo_seq:
                err = err or "empty_bmrb_holo_sequence"

        if not err and holo_seq:
            pdb_for_meta = normalize_pdb_id(row.get("holo_pdb") or "")
            if pdb_for_meta:
                holo_ctx = fetch_holo_polymer_context(pdb_for_meta, chain)
                rcsb_tax = holo_ctx.ncbi_taxonomy_id if holo_ctx and holo_ctx.ncbi_taxonomy_id else ""
                rcsb_refs = json.dumps((holo_ctx.uniprot_accessions if holo_ctx else []))
                pdb_entity_id_col = holo_ctx.entity_id or "" if holo_ctx else ""

            direct_acc_raw: Optional[str] = None
            if holo_ctx and getattr(holo_ctx, "preferred_uniprot_accession", None):
                direct_acc_raw = holo_ctx.preferred_uniprot_accession
            elif holo_ctx and holo_ctx.uniprot_accessions:
                direct_acc_raw = holo_ctx.uniprot_accessions[0]

            direct_acc_norm = normalize_uniprot_accession(direct_acc_raw or "")

            if direct_acc_norm:
                acc = direct_acc_norm
                rank_reason = "rcsb_polymer_entity_direct"
                blast_bs = "n/a"
                blast_id = "n/a"
                uniprot_source = "rcsb"
            else:
                blast_hits: List[BlastHit] = []
                uniprot_source = "blastp"
                try:
                    blast_hits = run_blastp(
                        holo_seq,
                        blast_db=blast_db,
                        blast_remote=blast_remote,
                        num_threads=min(8, os.cpu_count() or 4),
                        max_target_seqs=30,
                        log_command=not quiet,
                    )
                except Exception as ex:
                    err = f"blastp_failed:{ex}"

                if not err and not blast_hits:
                    err = "no_blast_hits"

                if blast_hits and not err:
                    query_len_b = len(holo_seq or "")
                    deposited_uni = holo_ctx.uniprot_accessions if holo_ctx else []
                    bh, rank_reason = _pick_uniprot_accession(
                        blast_hits,
                        holo_ctx.ncbi_taxonomy_id if holo_ctx else None,
                        deposited_uni,
                        query_len_b,
                        tax_cache=tax_cache,
                        min_pident=min_pident,
                        min_query_cov=min_query_cov,
                    )
                    if bh:
                        blast_bs = str(bh.bitscore)
                        blast_id = str(bh.pident)
                        acc = bh.subject_accession
                    else:
                        err = err or "no_uniprot_candidate"

            if acc and not err:
                if dry_run_rest:
                    uni_seq = "(dry-run)"
                    uni_len = len(holo_seq)
                else:
                    try:
                        uni_seq, uni_len = fetch_fasta_sequence(acc)
                        time.sleep(sleep_between_uniprot)
                    except Exception as ex:
                        err = f"uniprot_fetch_failed:{ex}"

        bh_len = len(holo_seq) if holo_seq else 0
        bas_len = len(apo_seq) if apo_seq else 0

        out_row["bmrb_holo_seq_source"] = holo_src
        out_row["bmrb_holo_seq_from_pdb_fallback"] = "1" if holo_pdb_fb else "0"
        out_row["bmrb_apo_seq_source"] = apo_src
        out_row["bmrb_apo_seq_from_pdb_fallback"] = "1" if apo_pdb_fb else "0"
        out_row["bmrb_holo_seq"] = holo_seq or ""
        out_row["bmrb_apo_seq"] = apo_seq or ""
        out_row["bmrb_holo_seq_length"] = str(bh_len)
        out_row["bmrb_apo_seq_length"] = str(bas_len)
        out_row["pdb_polymer_entity_id"] = pdb_entity_id_col
        out_row["uniprot_resolution_source"] = uniprot_source
        out_row["uniprot_accession"] = acc
        out_row["uniprot_seq"] = uni_seq if not dry_run_rest else "(dry-run)"
        out_row["uniprot_seq_length"] = str(uni_len)
        out_row["blast_bitscore"] = blast_bs
        out_row["blast_identity"] = blast_id
        out_row["blast_rank_reason"] = rank_reason
        out_row["rcsb_ncbi_taxonomy_id"] = str(rcsb_tax) if str(rcsb_tax) != "" else ""
        out_row["rcsb_reference_uniprots_json"] = rcsb_refs
        out_row["length_ratio_used"] = str(ratio)
        out_row["classification_rule"] = rule

        if err:
            out_row["bifurcation_label"] = ""
            out_row["error_reason"] = err
        elif bh_len <= 0 or uni_len <= 0:
            out_row["error_reason"] = err or "zero_length_sequences"
            out_row["bifurcation_label"] = ""
        else:
            out_row["error_reason"] = ""
            out_row["bifurcation_label"] = _classify_row(bh_len, uni_len, ratio)

        for c in extra_cols:
            out_row.setdefault(c, "")
        enriched.append(out_row)

        if not quiet:
            hp = holo_pdb or "-"
            ap = apo_pdb or "-"
            ch = (chain or "-").upper()
            errs = _short_err(out_row.get("error_reason"))
            if errs:
                status = f"COMPUTED UNRESOLVED  {errs}"
            else:
                lbl = (out_row["bifurcation_label"] or "—").upper()
                src = uniprot_source or "-"
                status = (
                    f"{lbl:13}"
                    f" src={src:6} UniProt:{acc}/L_uni={uni_len}/L_bmrb={bh_len}/reason:{rank_reason or '-'}"
                )
            print(
                f"[{i}/{n_input}] COMPUTE apo_bmrb={apo_bmrb} holo_bmrb={holo_bmrb} holo_pdb={hp} apo_pdb={ap} chain={ch} | {status}",
                file=sys.stdout,
            )

        bundle_result = write_bifurcation_bundle(
            enriched,
            base_fieldnames=base_fieldnames,
            basename=basename,
            out_dir=out_dir,
            started_at_perf=t0,
            skipped_cached=run_skipped_so_far,
        )

    if bundle_result is None:
        bundle_result = write_bifurcation_bundle(
            enriched,
            base_fieldnames=base_fieldnames,
            basename=basename,
            out_dir=out_dir,
            started_at_perf=t0,
            skipped_cached=run_skipped_so_far,
        )
    final = bundle_result

    if not quiet:
        print("", file=sys.stdout)
        print("Done.", file=sys.stdout)
        print(f"  Elapsed (wall): ~{final.get('elapsed_seconds', 0):.1f}s", file=sys.stdout)
        print(
            f"  Tallies:  domain={final.get('domains', 0)}, full_length={final.get('full_length', 0)}, "
            f"unresolved={final.get('unresolved', 0)} "
            f"(cached rows reused this run: {run_skipped_so_far})",
            file=sys.stdout,
        )
        print(f"  Wrote all outputs under basename «{basename}» (domains, full_length, unresolved, all_rows, JSON).")
    return final


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Bifurcate CSP_UBQ CSV into domain vs full_length.")
    p.add_argument(
        "--input",
        type=Path,
        default=Path(paths.filtered_input_csv),
        help=f"Filtered CSP CSV (default: {paths.filtered_input_csv})",
    )
    p.add_argument(
        "--receptor-chains",
        type=Path,
        default=Path(paths.receptor_chain_csv),
        help=f"Receptor chains CSV (default: {paths.receptor_chain_csv})",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path(paths.data_dir),
        help=f"Directory for bifurcation outputs (default: {paths.data_dir})",
    )
    p.add_argument(
        "--basename",
        default="CSP_UBQ_ph0.5_temp5C",
        help="Stem for output filenames (default: CSP_UBQ_ph0.5_temp5C)",
    )
    p.add_argument(
        "--cs-dir",
        type=Path,
        default=Path(paths.cs_cache_dir),
        help="BMRB NMR-STAR cache directory",
    )
    p.add_argument(
        "--pdb-cache-dir",
        default=paths.pdb_cache_dir,
        help="Cached PDB coordinates directory",
    )
    p.add_argument("--no-fetch-bmrb", action="store_true", help="Do not download missing BMRB STR files.")
    p.add_argument("--blast-db", default="swissprot", help="blastp -db (default swissprot; use nr, uniprot, etc.).")
    p.add_argument(
        "--blast-remote",
        action="store_true",
        help="Pass -remote to blastp (requires network; uses NCBI servers).",
    )
    p.add_argument("--min-pident", type=float, default=25.0, help="Minimum %% identity filter for accepting hits.")
    p.add_argument(
        "--min-query-coverage",
        type=float,
        default=0.5,
        help="Minimum blast query coverage (HSP extent / probe length)",
    )
    p.add_argument(
        "--length-ratio",
        type=float,
        default=LENGTH_RATIO_DEFAULT,
        help="Classify domain if uniprot length > ratio * bmrb holo length",
    )
    p.add_argument(
        "--sleep-uniprot",
        type=float,
        default=0.08,
        help="Seconds between UniProt REST calls (gentle pacing)",
    )
    p.add_argument(
        "--dry-run-rest",
        action="store_true",
        help="Still run STAR/PDB/blast parsing but stub UniProt sequence fetch (local testing).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress banner, per-row progress, and 'Done.' block; still prints the JSON summary at the end.",
    )
    p.add_argument(
        "--force-recompute",
        action="store_true",
        help=(
            "When set, ignore <basename>_bifurcation_all_rows.csv and recompute each row "
            "(default: merge/cache from that file for incremental/resumable runs)."
        ),
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    args = build_arg_parser().parse_args(argv)
    final = process_dataset(
        input_csv=args.input,
        receptor_csv=args.receptor_chains,
        out_dir=args.out_dir,
        basename=args.basename,
        cs_dir=args.cs_dir.resolve(),
        pdb_cache_dir=os.path.abspath(args.pdb_cache_dir),
        fetch_bmrb=not args.no_fetch_bmrb,
        blast_db=args.blast_db,
        blast_remote=args.blast_remote,
        min_pident=args.min_pident,
        min_query_cov=args.min_query_coverage,
        ratio=args.length_ratio,
        sleep_between_uniprot=args.sleep_uniprot,
        dry_run_rest=args.dry_run_rest,
        quiet=args.quiet,
        force_recompute=args.force_recompute,
    )
    if args.quiet:
        print(json.dumps(final, indent=2, sort_keys=True), file=sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
