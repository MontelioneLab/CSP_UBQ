"""
Receptor-centric multi-sequence alignment for pipeline visualization.

Uses the existing global pairwise aligner (``align_global`` from ``align.py``) and
a star-style merge on the longest (reference) sequence: gap columns on the
reference are merged by taking the column-wise maximum at each residue boundary,
then partner rows are lifted into the merged reference with block-wise padding.

This is **visualization-grade** (not MAFFT / Clustal quality) but deterministic
and free of external aligner binaries.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List, Optional, Tuple

try:
    from .align import align_global
except ImportError:
    from scripts.align import align_global


def row_key(apo_bmrb: str, holo_bmrb: str, holo_pdb: str) -> Tuple[str, str, str]:
    return (
        (apo_bmrb or "").strip(),
        (holo_bmrb or "").strip(),
        (holo_pdb or "").strip().lower(),
    )


def load_bifurcation_row(
    data_dir: Path,
    basename: str,
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
) -> Optional[Dict[str, str]]:
    """
    Return the first CSV row matching ``row_key`` searching ``{basename}_domains.csv``
    then ``{basename}_full_length.csv`` under ``data_dir``.
    """
    key = row_key(apo_bmrb, holo_bmrb, holo_pdb)
    for suffix in ("_domains.csv", "_full_length.csv"):
        path = data_dir / f"{basename}{suffix}"
        if not path.is_file():
            continue
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                rk = row_key(
                    row.get("apo_bmrb") or "",
                    row.get("holo_bmrb") or "",
                    row.get("holo_pdb") or "",
                )
                if rk == key:
                    return {k: ("" if v is None else str(v)) for k, v in row.items()}
    return None


def _verify_center(C: str, center: str) -> None:
    core = "".join(c for c in C if c != "-").upper()
    if core != center.upper():
        raise ValueError("aligned center fragment does not decode to reference")


def _gap_profile(C: str, center: str) -> List[int]:
    """
    Return ``gaps[k]`` = gap count before residue ``k`` for k in ``0..L-1``,
    and ``gaps[L]`` = trailing gaps after the last residue.
    """
    L = len(center)
    gaps = [0] * (L + 1)
    i = 0
    while i < len(C) and C[i] == "-":
        gaps[0] += 1
        i += 1
    for r in range(L):
        if i >= len(C) or C[i].upper() != center[r].upper():
            raise ValueError("truncated gapped center string or residue mismatch")
        i += 1
        if r < L - 1:
            g = 0
            while i < len(C) and C[i] == "-":
                g += 1
                i += 1
            gaps[r + 1] = g
        else:
            tr = 0
            while i < len(C) and C[i] == "-":
                tr += 1
                i += 1
            gaps[L] = tr
    if i != len(C):
        raise ValueError("extra characters in gapped center string")
    return gaps


def merge_reference_gaps(center: str, aligned_centers: List[str]) -> str:
    """Merge several gapped reference strings (all decode to ``center``) by max gaps per boundary."""
    L = len(center)
    if not aligned_centers:
        return center
    for C in aligned_centers:
        _verify_center(C, center)

    merged_gaps = [0] * (L + 1)
    for C in aligned_centers:
        g = _gap_profile(C, center)
        for k in range(L + 1):
            merged_gaps[k] = max(merged_gaps[k], g[k])

    parts: List[str] = []
    for r in range(L):
        parts.append("-" * merged_gaps[r])
        parts.append(center[r])
    parts.append("-" * merged_gaps[L])
    return "".join(parts)


def _residue_block_spans(C: str, center: str) -> List[Tuple[int, int]]:
    """Half-open [start,end) per residue: leading gaps for that residue + the residue letter."""
    L = len(center)
    spans: List[Tuple[int, int]] = []
    i = 0
    for r in range(L):
        start = i
        while i < len(C) and C[i] == "-":
            i += 1
        if i >= len(C) or C[i].upper() != center[r].upper():
            raise ValueError("premature end or residue mismatch in gapped reference")
        i += 1
        end = i
        spans.append((start, end))
    while i < len(C) and C[i] == "-":
        i += 1
    if i != len(C):
        raise ValueError("trailing non-gap content in gapped reference")
    return spans


def _lift_block(seg_m: str, seg_c: str, seg_p: str) -> str:
    """Pad / lift partner segment to match wider reference segment ``seg_m``."""
    if len(seg_c) != len(seg_p):
        raise ValueError("pair seg length mismatch")
    if len(seg_m) < len(seg_c):
        raise ValueError("reference block shorter than pair block")
    out: List[str] = []
    ii, jj = 0, 0
    while ii < len(seg_m):
        mch = seg_m[ii]
        cch = seg_c[jj]
        pch = seg_p[jj]
        if mch == "-":
            if cch == "-":
                out.append(pch)
                ii += 1
                jj += 1
            else:
                out.append("-")
                ii += 1
        else:
            if cch == "-":
                out.append(pch)
                ii += 1
                jj += 1
            else:
                if mch.upper() != cch.upper():
                    raise ValueError("block residue mismatch")
                out.append(pch)
                ii += 1
                jj += 1
    if jj != len(seg_c):
        raise ValueError("lift_block did not consume pair block")
    return "".join(out)


def lift_partner_to_master(master: str, pair_c: str, pair_p: str, center: str) -> str:
    """Map ``pair_p`` (aligned to ``pair_c``) onto columns of ``master`` (merged ref)."""
    _verify_center(pair_c, center)
    if len(pair_c) != len(pair_p):
        raise ValueError("pair_c / pair_p length mismatch")
    spans_m = _residue_block_spans(master, center)
    spans_c = _residue_block_spans(pair_c, center)
    out_chars: List[str] = []
    for (sm, em), (sc, ec) in zip(spans_m, spans_c):
        seg_m = master[sm:em]
        seg_c = pair_c[sc:ec]
        seg_p = pair_p[sc:ec]
        out_chars.append(_lift_block(seg_m, seg_c, seg_p))
    return "".join(out_chars)


def star_msa_strings(
    labeled_seqs: List[Tuple[str, str]],
) -> List[Tuple[str, str]]:
    """
    Build equal-length aligned strings from (label, ungapped_sequence) pairs.

    Skips sequences that are empty after strip. Requires at least two non-empty
    sequences. The longest sequence is the star center.
    """
    cleaned: List[Tuple[str, str]] = []
    for lab, seq in labeled_seqs:
        s = (seq or "").strip().upper().replace("\n", "")
        if not s:
            continue
        cleaned.append((lab, s))
    if len(cleaned) < 2:
        raise ValueError("need at least two non-empty sequences for MSA")
    center_idx = max(range(len(cleaned)), key=lambda i: (len(cleaned[i][1]), -i))
    center_label, center = cleaned[center_idx]

    aligned_centers: List[str] = []
    pair_data: List[Tuple[str, str, str]] = []
    for lab, seq in cleaned:
        ca, pb, _, _ = align_global(center, seq)
        aligned_centers.append(ca)
        pair_data.append((lab, ca, pb))

    master = merge_reference_gaps(center, aligned_centers)

    out: List[Tuple[str, str]] = []
    for i, (lab, ca, pb) in enumerate(pair_data):
        if i == center_idx:
            row_s = master
        else:
            row_s = lift_partner_to_master(master, ca, pb, center)
        if len(row_s) != len(master):
            raise RuntimeError(f"row width {len(row_s)} != master {len(master)}")
        out.append((lab, row_s))
    return out


def write_msa_png(
    aligned_rows: List[Tuple[str, str]],
    out_path: Path,
    *,
    title: str = "",
    wrap_width: int = 90,
    dpi: int = 120,
) -> None:
    """Render aligned (label, gapped_seq) rows to a monospace PNG; long lines wrap in blocks."""
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mono = font_manager.FontProperties(family="monospace", size=8)
    lines: List[str] = []
    for lab, gseq in aligned_rows:
        block_start = 0
        while block_start < len(gseq):
            chunk = gseq[block_start : block_start + wrap_width]
            prefix = f"{lab:16s} " if block_start == 0 else f"{'':16s} "
            lines.append(prefix + chunk)
            block_start += wrap_width

    fig_h = max(2.0, 0.18 * len(lines) + (1.2 if title else 0.4))
    fig_w = min(18.0, max(8.0, wrap_width * 0.055 + 1.5))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10, family="sans-serif", loc="left")
    nlines = max(len(lines), 1)
    for i, ln in enumerate(lines):
        ax.text(0.0, nlines - 1 - i, ln, fontproperties=mono, va="top", fontsize=8)
    ax.set_xlim(0, wrap_width + 24)
    ax.set_ylim(-0.5, nlines + 0.6)
    ax.invert_yaxis()
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def try_write_receptor_alignment_png(
    *,
    data_dir: Path,
    bifurcation_basename: str,
    apo_bmrb: str,
    holo_bmrb: str,
    holo_pdb: str,
    receptor_chain: Optional[str],
    seq_apo_pdb: Optional[str],
    seq_holo_pdb: Optional[str],
    out_png: Path,
) -> Tuple[bool, str]:
    """
    Load bifurcation row; collect up to five sequences; write PNG. Returns (ok, message).
    """
    row = load_bifurcation_row(data_dir, bifurcation_basename, apo_bmrb, holo_bmrb, holo_pdb)
    if not row:
        return False, "no bifurcation row in domains/full_length CSV"
    uni = (row.get("uniprot_seq") or "").strip()
    if not uni or uni == "(dry-run)":
        return False, "empty or dry-run uniprot_seq in bifurcation row"
    apo_bmrb_seq = (row.get("bmrb_apo_seq") or "").strip()
    holo_bmrb_seq = (row.get("bmrb_holo_seq") or "").strip()
    if not apo_bmrb_seq or not holo_bmrb_seq:
        return False, "missing bmrb_apo_seq or bmrb_holo_seq in bifurcation row"

    ch_bif = (row.get("receptor_chain_id") or "").strip().lower()
    ch_pipe = (receptor_chain or "").strip().lower()
    if ch_bif and ch_pipe and ch_bif != ch_pipe:
        # non-fatal
        pass

    labeled: List[Tuple[str, str]] = []
    if seq_apo_pdb:
        labeled.append(("apo_PDB", seq_apo_pdb))
    if seq_holo_pdb:
        labeled.append(("holo_PDB", seq_holo_pdb))
    labeled.append(("apo_BMRB", apo_bmrb_seq))
    labeled.append(("holo_BMRB", holo_bmrb_seq))
    acc = (row.get("uniprot_accession") or "").strip()
    uni_label = f"UniProt({acc})" if acc else "UniProt"
    labeled.append((uni_label, uni))

    if len([x for x in labeled if x[1]]) < 2:
        return False, "fewer than two sequences after filtering"

    try:
        aligned = star_msa_strings(labeled)
    except Exception as ex:
        return False, f"MSA failed: {ex}"

    note = "progressive star MSA (union gaps on longest ref); visualization-grade"
    title = f"{holo_pdb.upper()}  |  {note}"
    if not seq_apo_pdb:
        title += "  |  apo PDB missing"
    write_msa_png(aligned, out_png, title=title)
    return True, f"wrote {out_png}"


__all__ = [
    "load_bifurcation_row",
    "star_msa_strings",
    "merge_reference_gaps",
    "lift_partner_to_master",
    "write_msa_png",
    "try_write_receptor_alignment_png",
    "row_key",
]
