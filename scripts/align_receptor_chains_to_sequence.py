#!/usr/bin/env python3
"""
Align receptor chains from PDB files to a user-provided sequence.

Loops through PDB files in a directory (default: PDB_FILES), extracts the
receptor chain (longest chain per PDB) or all chains if requested, performs
global sequence alignment to the user-provided query sequence using
Biopython + BLOSUM62, and writes a CSV table of results sorted by E-value
(best first).

Requires: Biopython (pip install biopython).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

try:
    from .config import paths
except Exception:
    import sys

    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, _root)
    from scripts.config import paths

from scripts import rcsb_io


def _clean_sequence(seq: str) -> str:
    """Normalize to 20 canonical one-letter amino acids; drop other characters."""
    if not seq:
        return ""
    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.strip().upper().replace(" ", "")
    return "".join(ch for ch in seq if ch in allowed)


def _approximate_e_value(score: float, m: int, n: int) -> float:
    """Approximate E-value: E = K * m * n * exp(-lambda * S) for BLOSUM62."""
    if m <= 0 or n <= 0:
        return float("nan")
    lambda_ = 0.318
    K = 0.13
    exponent = -lambda_ * score
    if exponent < -700:
        return 0.0
    if exponent > 700:
        return float("inf")
    return K * m * n * math.exp(exponent)


def _align_and_score(
    query_seq: str,
    chain_seq: str,
) -> Tuple[float, int, float]:
    """Global alignment (BLOSUM62); returns (score, overlap_length, identity_fraction)."""
    if not query_seq or not chain_seq:
        return float("nan"), 0, float("nan")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11.0
    aligner.extend_gap_score = -1.0
    alignments = list(aligner.align(query_seq, chain_seq))
    if not alignments:
        return float("nan"), 0, float("nan")
    best = alignments[0]
    score = float(best.score)
    # best[0] and best[1] are the two aligned sequences (with gaps)
    aligned_q = str(best[0])
    aligned_c = str(best[1])
    overlap = 0
    identities = 0
    for a, b in zip(aligned_q, aligned_c):
        if a == "-" or b == "-":
            continue
        overlap += 1
        if a == b:
            identities += 1
    identity_frac = identities / overlap if overlap else float("nan")
    return score, overlap, identity_frac


def _longest_chain(chain_sequences: Dict[str, str]) -> str:
    """Return chain ID with longest sequence; arbitrary tie-break."""
    if not chain_sequences:
        return ""
    return max(chain_sequences.keys(), key=lambda c: len(chain_sequences[c]))


def run(
    query_sequence: str,
    pdb_dir: Path,
    output_csv: Path,
    all_chains: bool = False,
) -> List[dict]:
    """Align receptor (or all) chains from PDBs to query; return rows, sorted by E-value."""
    query_clean = _clean_sequence(query_sequence)
    if not query_clean:
        raise ValueError("Query sequence is empty or has no valid amino acids.")

    pdb_dir = Path(pdb_dir)
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    rows: List[dict] = []

    for pdb_path in pdb_files:
        pdb_id = pdb_path.stem
        try:
            chain_seqs = rcsb_io.parse_pdb_sequences(str(pdb_path))
        except Exception:
            continue
        if not chain_seqs:
            continue

        if all_chains:
            chains_to_use = list(chain_seqs.items())
        else:
            rec_chain = _longest_chain(chain_seqs)
            if not rec_chain:
                continue
            chains_to_use = [(rec_chain, chain_seqs[rec_chain])]

        for chain_id, chain_seq in chains_to_use:
            chain_clean = _clean_sequence(chain_seq)
            if not chain_clean:
                continue
            score, overlap, identity_frac = _align_and_score(query_clean, chain_clean)
            e_val = _approximate_e_value(score, len(query_clean), len(chain_clean))
            rows.append({
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "query_seq_length": len(query_clean),
                "chain_seq_length": len(chain_clean),
                "alignment_score": score,
                "overlap_length": overlap,
                "identity_fraction": identity_frac,
                "approx_e_value": e_val,
            })

    # Sort by E-value ascending (best first); NaN/Inf last
    def _e_sort_key(r):
        v = r["approx_e_value"]
        if math.isnan(v):
            return (2, 0.0)
        if math.isinf(v) and v > 0:
            return (2, 0.0)
        return (0, v)

    rows.sort(key=_e_sort_key)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id", "chain_id", "query_seq_length", "chain_seq_length",
        "alignment_score", "overlap_length", "identity_fraction", "approx_e_value",
    ]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            out = {k: r[k] for k in fieldnames}
            out["alignment_score"] = f"{r['alignment_score']:.3f}"
            out["identity_fraction"] = f"{r['identity_fraction']:.4f}" if not math.isnan(r["identity_fraction"]) else ""
            ev = r["approx_e_value"]
            out["approx_e_value"] = f"{ev:.3e}" if math.isfinite(ev) else ""
            w.writerow(out)

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Align receptor chains from PDB files to a user-provided sequence; output CSV sorted by E-value.",
    )
    parser.add_argument(
        "sequence",
        type=str,
        help="Query sequence (one-letter amino acids).",
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=Path(paths.pdb_cache_dir),
        help="Directory containing PDB files (default: PDB_FILES).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("receptor_chain_alignments.csv"),
        help="Output CSV path (default: receptor_chain_alignments.csv).",
    )
    parser.add_argument(
        "--all-chains",
        action="store_true",
        help="Align every chain in each PDB; default is longest chain only (receptor).",
    )
    args = parser.parse_args()

    rows = run(
        query_sequence=args.sequence,
        pdb_dir=args.pdb_dir,
        output_csv=args.output,
        all_chains=args.all_chains,
    )
    print(f"Wrote {len(rows)} alignment(s) to {args.output} (sorted by E-value).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
