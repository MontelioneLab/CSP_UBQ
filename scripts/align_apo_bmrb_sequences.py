"""
Align apo BMRB sequences against the apo BMRB 17769 sequence.

Reads CSP_UBQ.csv from the repository root, identifies the row with apo_bmrb == 17769
and uses its match_seq as the reference sequence. All other rows with a non-empty
apo_bmrb and match_seq are globally aligned against this reference using Biopython
pairwise2 + BLOSUM62. For each comparison, the script reports:

- alignment score
- overlap length (non-gap positions in the aligned pair)
- identity fraction over the overlap
- an approximate E-value (Karlin–Altschul style)

Results are written to apo_bmrb_17769_alignments.csv in the repository root.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import warnings
from dataclasses import dataclass
from typing import Iterable, List, Tuple

try:
    from Bio import pairwise2
    from Bio.Align import substitution_matrices
except ImportError as exc:
    raise ImportError(
        "Biopython is required for sequence alignment. Install with: pip install biopython"
    ) from exc

from scripts import rcsb_io


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_CSV = os.path.join(ROOT_DIR, "CSP_UBQ.csv")
OUTPUT_CSV = os.path.join(ROOT_DIR, "apo_bmrb_17769_alignments.csv")
REF_APO_BMRB_ID = "17769"


@dataclass
class AlignmentResult:
    ref_apo_bmrb: str
    ref_seq_length: int
    query_apo_bmrb: str
    query_row_index: int
    query_seq_length: int
    alignment_score: float
    aligned_length: int
    overlap_length: int
    identity_count: int
    identity_fraction: float
    approx_e_value: float


def _allowed_sequence_alphabet() -> set[str]:
    """
    Keep one-letter symbols consistent with parse_pdb_sequences in rcsb_io.py.

    parse_pdb_sequences maps standard residues to the 20 canonical amino acids and
    non-standard residues to X, so we accept the same alphabet here.
    """
    _ = rcsb_io  # ensure module import is intentional and not removed
    return set("ACDEFGHIKLMNPQRSTVWYX")


def _clean_sequence(seq: str) -> str:
    """
    Normalize a sequence string to a standard one-letter amino acid alphabet.

    - Strips whitespace
    - Uppercases letters
    - Keeps only the canonical amino acid alphabet used by parse_pdb_sequences
      plus X for unknown/non-standard residues
    """
    if not seq:
        return ""
    allowed = _allowed_sequence_alphabet()
    seq = "".join(seq.strip().upper().split())
    return "".join(ch for ch in seq if ch in allowed)


def _approximate_e_value(score: float, m: int, n: int) -> float:
    """
    Approximate a BLAST-style E-value from an alignment score using
    Karlin–Altschul statistics:

        E = K * m * n * exp(-lambda * S)

    where:
    - S is the alignment score
    - m, n are the sequence lengths
    - lambda, K are empirical constants for BLOSUM62 (approximate).
    """
    if m <= 0 or n <= 0:
        return float("nan")
    # Approximate parameters often used for BLOSUM62
    lambda_ = 0.318
    K = 0.13
    exponent = -lambda_ * score
    # Guard against extreme underflow/overflow
    if exponent < -700:
        return 0.0
    if exponent > 700:
        return float("inf")
    return K * m * n * math.exp(exponent)


def align_and_score(
    ref_seq: str,
    query_seq: str,
) -> Tuple[float, int, int, int, float]:
    """
    Run a global alignment (Needleman–Wunsch) and compute:

    - alignment_score
    - aligned_length
    - overlap_length (non-gap positions)
    - identity_count over overlap
    - identity_fraction over the overlap
    """
    if not ref_seq or not query_seq:
        return float("nan"), 0, 0, 0, float("nan")

    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -11.0
    gap_extend = -1.0

    alignments = pairwise2.align.globalds(
        ref_seq,
        query_seq,
        matrix,
        gap_open,
        gap_extend,
        one_alignment_only=True,
    )
    if not alignments:
        return float("nan"), 0, 0, 0, float("nan")

    best = alignments[0]
    aligned_ref = best.seqA
    aligned_query = best.seqB
    score = float(best.score)
    aligned_length = len(aligned_ref)

    overlap = 0
    identities = 0
    for a, b in zip(aligned_ref, aligned_query):
        if a == "-" or b == "-":
            continue
        overlap += 1
        if a == b:
            identities += 1

    if overlap == 0:
        identity_fraction = float("nan")
    else:
        identity_fraction = identities / overlap

    return score, aligned_length, overlap, identities, identity_fraction


def _iter_apo_rows(path: str) -> Iterable[Tuple[int, dict]]:
    """
    Yield (1-based row_index, row_dict) for rows in CSP_UBQ.csv that have
    a non-empty apo_bmrb field.
    """
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=2):  # header is line 1
            apo_id = (row.get("apo_bmrb") or "").strip()
            if not apo_id:
                continue
            yield idx, row


def _find_reference_sequence(path: str, ref_apo_id: str) -> Tuple[str, int]:
    """
    Look up the reference sequence and its row index for apo_bmrb == ref_apo_id.
    """
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=2):
            apo_id = (row.get("apo_bmrb") or "").strip()
            if apo_id != ref_apo_id:
                continue
            seq_raw = row.get("match_seq") or ""
            seq = _clean_sequence(seq_raw)
            if not seq:
                raise ValueError(
                    f"Reference row for apo_bmrb {ref_apo_id} has empty/invalid match_seq"
                )
            return seq, idx
    raise ValueError(f"No row found in {path} with apo_bmrb == {ref_apo_id}")


def run_alignments(
    input_csv: str = INPUT_CSV,
    output_csv: str = OUTPUT_CSV,
    ref_apo_id: str = REF_APO_BMRB_ID,
) -> List[AlignmentResult]:
    """
    Main worker: run alignments of all apo sequences against the reference.
    """
    ref_seq, ref_row_idx = _find_reference_sequence(input_csv, ref_apo_id)
    ref_len = len(ref_seq)

    results: List[AlignmentResult] = []
    warned_missing_seq = 0
    for row_idx, row in _iter_apo_rows(input_csv):
        apo_id = (row.get("apo_bmrb") or "").strip()
        if apo_id == ref_apo_id:
            # Skip self
            continue
        seq_raw = row.get("match_seq") or ""
        seq = _clean_sequence(seq_raw)
        if not seq:
            warned_missing_seq += 1
            warnings.warn(
                f"Skipping row {row_idx} (apo_bmrb={apo_id}): empty/invalid match_seq after cleaning.",
                RuntimeWarning,
                stacklevel=2,
            )
            continue

        score, aligned_length, overlap, identity_count, identity_fraction = align_and_score(
            ref_seq, seq
        )
        e_value = _approximate_e_value(score, ref_len, len(seq))

        results.append(
            AlignmentResult(
                ref_apo_bmrb=ref_apo_id,
                ref_seq_length=ref_len,
                query_apo_bmrb=apo_id,
                query_row_index=row_idx,
                query_seq_length=len(seq),
                alignment_score=score,
                aligned_length=aligned_length,
                overlap_length=overlap,
                identity_count=identity_count,
                identity_fraction=identity_fraction,
                approx_e_value=e_value,
            )
        )

    # Write CSV
    fieldnames = [
        "ref_apo_bmrb",
        "ref_seq_length",
        "query_apo_bmrb",
        "query_row_index",
        "query_seq_length",
        "alignment_score",
        "aligned_length",
        "overlap_length",
        "identity_count",
        "identity_fraction",
        "approx_e_value",
    ]
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(
                {
                    "ref_apo_bmrb": r.ref_apo_bmrb,
                    "ref_seq_length": r.ref_seq_length,
                    "query_apo_bmrb": r.query_apo_bmrb,
                    "query_row_index": r.query_row_index,
                    "query_seq_length": r.query_seq_length,
                    "alignment_score": f"{r.alignment_score:.3f}",
                    "aligned_length": r.aligned_length,
                    "overlap_length": r.overlap_length,
                    "identity_count": r.identity_count,
                    "identity_fraction": (
                        f"{r.identity_fraction:.4f}"
                        if not math.isnan(r.identity_fraction)
                        else ""
                    ),
                    "approx_e_value": (
                        f"{r.approx_e_value:.3e}"
                        if math.isfinite(r.approx_e_value)
                        else ""
                    ),
                }
            )

    if warned_missing_seq:
        print(f"Skipped {warned_missing_seq} rows with apo_bmrb set but empty/invalid match_seq.")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Align all apo_bmrb sequences in CSP_UBQ.csv to the apo_bmrb 17769 reference "
            "and report overlap and approximate E-values."
        )
    )
    parser.add_argument("--input", default=INPUT_CSV, help="Input CSV path (default: CSP_UBQ.csv)")
    parser.add_argument(
        "--output",
        default=OUTPUT_CSV,
        help="Output CSV path (default: apo_bmrb_17769_alignments.csv)",
    )
    parser.add_argument(
        "--ref-apo-id",
        default=REF_APO_BMRB_ID,
        help="Reference apo_bmrb ID (default: 17769)",
    )
    parser.add_argument(
        "--summary-count",
        type=int,
        default=5,
        help="Number of alignment rows to summarize to stdout (default: 5)",
    )
    args = parser.parse_args()

    print(f"Input CSV: {args.input}")
    print(f"Output CSV: {args.output}")
    print(f"Reference apo_bmrb: {args.ref_apo_id}")
    results = run_alignments(input_csv=args.input, output_csv=args.output, ref_apo_id=args.ref_apo_id)
    print(f"Wrote {len(results)} alignment rows to {args.output}")

    preview_count = max(0, min(args.summary_count, len(results)))
    if preview_count:
        print("Alignment summary preview:")
        for result in results[:preview_count]:
            print(
                f"  row={result.query_row_index}, apo_bmrb={result.query_apo_bmrb}, "
                f"score={result.alignment_score:.2f}, overlap={result.overlap_length}, "
                f"identity={result.identity_fraction:.3f}, e={result.approx_e_value:.3e}"
            )


if __name__ == "__main__":
    main()

