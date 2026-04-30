"""
Run blastp against a protein sequence and parse tabular output.
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import List, Optional

try:
    from .uniprot_io import normalize_uniprot_accession
except ImportError:
    from scripts.uniprot_io import normalize_uniprot_accession
@dataclass
class BlastHit:
    """One row from blastp ``-outfmt 6`` (default classic fields)."""

    subject_accession: str  # normalized UniProt-like accession
    pident: float
    length: int  # alignment length
    bitscore: float
    qstart: int
    qend: int
    sstart: int
    send: int
    subject_title: str


def _parse_blast_tab_line(line: str) -> Optional[BlastHit]:
    line = line.strip()
    if not line or line.startswith("#") or line.startswith("Query=") or line.startswith("Field"):
        return None
    parts = line.split("\t")
    if len(parts) < 12:
        parts = line.split()
    if len(parts) < 12:
        return None
    try:
        sseqid = parts[1]
        pident = float(parts[2])
        length = int(parts[3])
        # parts[4-9] mismatches, gaps, qstart, qend, sstart, send
        qstart, qend, sstart, send = (
            int(parts[6]),
            int(parts[7]),
            int(parts[8]),
            int(parts[9]),
        )
        bitscore = float(parts[11])
        title = " ".join(parts[12:]) if len(parts) > 12 else ""
        acc = normalize_uniprot_accession(sseqid)
        if not acc:
            return None
        return BlastHit(
            subject_accession=acc,
            pident=pident,
            length=length,
            bitscore=bitscore,
            qstart=qstart,
            qend=qend,
            sstart=sstart,
            send=send,
            subject_title=title,
        )
    except (ValueError, IndexError):
        return None


DEFAULT_OUTFMT = (
    "qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
)


def run_blastp(
    query_sequence: str,
    *,
    blast_db: str = "nr",
    blast_remote: bool = False,
    num_threads: int = 2,
    max_target_seqs: int = 25,
    expect: float = 10.0,
    blastp_executable: str = "blastp",
    extra_args: Optional[List[str]] = None,
    log_command: bool = True,
) -> List[BlastHit]:
    """Run ``blastp`` with tabular output and return parsed hits (best-effort).

    When ``log_command`` is true, prints the full ``blastp`` argv to stdout (shell-quoted)
    immediately before the subprocess runs.
    """
    seq = (query_sequence or "").strip().upper().replace(" ", "").replace("\n", "")
    if len(seq) < 5:
        return []

    extra_args = list(extra_args or [])
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False, encoding="ascii") as fq:
        fq.write(">query\n")
        fq.write(seq + "\n")
        qpath = fq.name

    try:
        # NCBI blastp rejects -num_threads together with -remote (prints USAGE, exit 1).
        cmd: List[str] = [
            blastp_executable,
            "-query",
            qpath,
            "-db",
            blast_db,
            "-outfmt",
            f"6 {DEFAULT_OUTFMT}",
            "-max_target_seqs",
            str(max_target_seqs),
            "-evalue",
            str(expect),
        ]
        if not blast_remote:
            cmd.extend(["-num_threads", str(num_threads)])
        if blast_remote:
            cmd.append("-remote")
        cmd.extend(extra_args)

        if log_command:
            print(f"blastp: {shlex.join(cmd)}", file=sys.stdout, flush=True)

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,
            check=False,
            env={**os.environ.copy(), "PATH": os.environ.get("PATH", "")},
        )
        stdout = proc.stdout or ""
        stderr = proc.stderr or ""
        if proc.returncode != 0:
            raise RuntimeError(f"blastp failed (code {proc.returncode}): {stderr[:500]}")

        hits: List[BlastHit] = []
        for line in stdout.splitlines():
            h = _parse_blast_tab_line(line)
            if h:
                hits.append(h)
        if log_command:
            print(
                f"blastp: finished (exit 0), {len(hits)} hit(s) parsed from tabular output",
                file=sys.stdout,
                flush=True,
            )
        return hits
    finally:
        try:
            os.unlink(qpath)
        except OSError:
            pass


def best_coverage_fraction(hit: BlastHit, query_len: int) -> float:
    """Fraction of query covered by HSP (coarse)."""
    if query_len <= 0:
        return 0.0
    return max(0.0, min(1.0, (hit.qend - hit.qstart + 1) / query_len))
