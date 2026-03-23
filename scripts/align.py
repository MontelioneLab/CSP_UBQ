"""
Global alignment utilities to map apo and holo sequences.

Prefers Biopython's PairwiseAligner when available; otherwise falls back to a
simple Needleman–Wunsch implementation with linear gap penalties.
"""

from __future__ import annotations

from typing import List, Tuple

try:
    from Bio import Align
    _HAS_BIOPYTHON = True
except Exception:  # Biopython not installed
    _HAS_BIOPYTHON = False

try:
    from .config import alignment as cfg
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import alignment as cfg


def align_global(apo_seq: str, holo_seq: str) -> Tuple[str, str, List[Tuple[int, int]], float]:
    """Align sequences and return alignment strings and index mapping.

    Returns
    -------
    (aligned_apo, aligned_holo, mapping, score)
        mapping: list of (apo_index, holo_index) 1-based indices for matched (non-gap) columns
    """
    # Biopython raises ValueError for empty sequences; handle explicitly
    if not apo_seq or not holo_seq:
        if not apo_seq and not holo_seq:
            return "", "", [], 0.0
        if not apo_seq:
            gaps = "-" * len(holo_seq)
            return gaps, holo_seq, [], len(holo_seq) * cfg.gap_open_penalty
        gaps = "-" * len(apo_seq)
        return apo_seq, gaps, [], len(apo_seq) * cfg.gap_open_penalty
    if _HAS_BIOPYTHON:
        return _align_biopython(apo_seq, holo_seq)
    return _align_nw(apo_seq, holo_seq)


def _align_biopython(a: str, b: str) -> Tuple[str, str, List[Tuple[int, int]]]:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = cfg.match_score
    aligner.mismatch_score = cfg.mismatch_penalty
    aligner.open_gap_score = cfg.gap_open_penalty
    aligner.extend_gap_score = cfg.gap_extend_penalty
    # Use iterator and get first alignment (all optimal alignments have the same score)
    # This avoids materializing millions of alignments when sequences are very different
    alignment_iter = aligner.align(a, b)
    best = next(alignment_iter)
    aligned_a = str(best.aligned[0])  # coordinates; need string sequences
    # Use aligner to get aligned strings via format
    aligned_str = best.format()
    # Parse the format output more carefully
    lines = [ln.rstrip("\n") for ln in aligned_str.splitlines() if ln.strip()]
    
    # Extract sequences from format output
    # Format is: query line, match line, target line (repeated)
    apo_aln: List[str] = []
    holo_aln: List[str] = []
    
    for i in range(0, len(lines), 3):
        if i + 2 < len(lines):
            # Extract just the sequence part (after the colon and position info)
            query_line = lines[i]
            target_line = lines[i + 2]
            
            # Find the sequence part after the colon
            if ':' in query_line:
                query_seq = query_line.split(':', 1)[1].strip()
            else:
                query_seq = query_line.strip()
                
            if ':' in target_line:
                target_seq = target_line.split(':', 1)[1].strip()
            else:
                target_seq = target_line.strip()
            
            apo_aln.append(query_seq)
            holo_aln.append(target_seq)
    
    aligned_apo = "".join(apo_aln)
    aligned_holo = "".join(holo_aln)
    
    
    # Clean the alignment strings to remove any formatting artifacts
    aligned_apo = _clean_alignment_sequence(aligned_apo)
    aligned_holo = _clean_alignment_sequence(aligned_holo)
    
    mapping = _build_mapping_from_alignment(aligned_apo, aligned_holo)
    
    # Debug output for Biopython alignment
    import os as _os
    if (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[ALIGN] Using Biopython PairwiseAligner")
        print(f"[ALIGN] Alignment score: {best.score}")
    
    return aligned_apo, aligned_holo, mapping, best.score


def _align_nw(a: str, b: str) -> Tuple[str, str, List[Tuple[int, int]], float]:
    # Simple Needleman–Wunsch with linear gap penalties
    match = cfg.match_score
    mismatch = cfg.mismatch_penalty
    gap = cfg.gap_open_penalty  # linear gap cost

    na, nb = len(a), len(b)
    score = [[0] * (nb + 1) for _ in range(na + 1)]
    trace = [[0] * (nb + 1) for _ in range(na + 1)]  # 0 diag, 1 up, 2 left
    for i in range(1, na + 1):
        score[i][0] = i * gap
        trace[i][0] = 1
    for j in range(1, nb + 1):
        score[0][j] = j * gap
        trace[0][j] = 2
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            s = match if a[i-1] == b[j-1] else mismatch
            diag = score[i-1][j-1] + s
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            best, move = max((diag, 0), (up, 1), (left, 2), key=lambda x: x[0])
            score[i][j] = best
            trace[i][j] = move
    # traceback
    i, j = na, nb
    a_aln: List[str] = []
    b_aln: List[str] = []
    while i > 0 or j > 0:
        mv = trace[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and mv == 0:
            a_aln.append(a[i-1])
            b_aln.append(b[j-1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or mv == 1):
            a_aln.append(a[i-1])
            b_aln.append("-")
            i -= 1
        else:
            a_aln.append("-")
            b_aln.append(b[j-1])
            j -= 1
    aligned_apo = "".join(reversed(a_aln))
    aligned_holo = "".join(reversed(b_aln))
    
    # Clean the alignment strings to remove any formatting artifacts
    aligned_apo = _clean_alignment_sequence(aligned_apo)
    aligned_holo = _clean_alignment_sequence(aligned_holo)
    
    mapping = _build_mapping_from_alignment(aligned_apo, aligned_holo)
    
    # Debug output for Needleman-Wunsch alignment
    import os as _os
    if (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[ALIGN] Using custom Needleman-Wunsch implementation")
        print(f"[ALIGN] Alignment score: {score[na][nb]}")
    
    return aligned_apo, aligned_holo, mapping, score[na][nb]


def _clean_alignment_sequence(seq: str) -> str:
    """Remove any non A-Z (uppercase) characters or '-' from the sequence."""
    cleaned = "".join(c for c in seq if (c == '-' or ('A' <= c <= 'Z')))
    return cleaned


def _build_mapping_from_alignment(aligned_apo: str, aligned_holo: str) -> List[Tuple[int, int]]:
    mapping: List[Tuple[int, int]] = []
    i_apo = 0
    i_holo = 0
    for ca, ch in zip(aligned_apo, aligned_holo):
        if ca != "-":
            i_apo += 1
        if ch != "-":
            i_holo += 1
        if ca != "-" and ch != "-":
            mapping.append((i_apo, i_holo))
    # Debug summary
    if (len(aligned_apo) > 0):
        # Avoid import cycle by local env check
        import os as _os
        if (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[ALIGN] alignment length {len(aligned_apo)}; mapped pairs {len(mapping)}")
            _print_neat_alignment(aligned_apo, aligned_holo)
    return mapping


def _print_neat_alignment(seq1: str, seq2: str, line_width: int = 80) -> None:
    """Print a neatly formatted alignment with line numbers and match indicators."""
    match_string = ''.join('|' if a == b and a != '-' else ' ' for a, b in zip(seq1, seq2))
    
    print("\n[ALIGN] Sequence Alignment:")
    print("=" * (line_width + 20))
    
    # Split into chunks for better readability
    for i in range(0, len(seq1), line_width):
        chunk1 = seq1[i:i+line_width]
        chunk2 = seq2[i:i+line_width]
        chunk_match = match_string[i:i+line_width]
        
        # Add line numbers
        start_pos = i + 1
        end_pos = min(i + line_width, len(seq1))
        
        print(f"[ALIGN] Apo  {start_pos:4d}-{end_pos:4d}: {chunk1}")
        print(f"[ALIGN]      {' ' * 10} {chunk_match}")
        print(f"[ALIGN] Holo {start_pos:4d}-{end_pos:4d}: {chunk2}")
        print()
    
    # Calculate and display statistics
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
    gaps1 = sum(1 for c in seq1 if c == '-')
    gaps2 = sum(1 for c in seq2 if c == '-')
    mismatches = len(seq1) - matches - gaps1 - gaps2
    
    print(f"[ALIGN] Statistics:")
    print(f"[ALIGN]   Matches:    {matches:4d}")
    print(f"[ALIGN]   Mismatches: {mismatches:4d}")
    print(f"[ALIGN]   Gaps (apo): {gaps1:4d}")
    print(f"[ALIGN]   Gaps (holo):{gaps2:4d}")
    print(f"[ALIGN]   Identity:   {matches/(len(seq1)-gaps1-gaps2)*100:.1f}%")
    print("=" * (line_width + 20))


