"""
Minimal UniProt REST helpers (sequence and organism taxonomy).
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import requests

try:
    from .config import network
except Exception:
    import os as _os, sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import network

_VALID_AA = frozenset("ACDEFGHIKLMNPQRSTVWYUXOBZ")


def _parse_uniprot_fasta_body(text: str) -> str:
    """Join all sequence lines after the first FASTA header (UniProt uses fixed-width lines)."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return ""
    i = 0
    if lines[0].startswith(">"):
        i = 1
    chunks: list[str] = []
    for ln in lines[i:]:
        if ln.startswith(">"):
            break
        chunks.append(ln)
    raw = "".join(chunks).replace("*", "").upper()
    return "".join(a for a in raw if a in _VALID_AA)


def normalize_uniprot_accession(raw: str) -> str:
    """Strip version and isolate accession (handles sp|P12345|NAME)."""
    if not raw:
        return ""
    s = raw.strip()
    if "|" in s:
        parts = [p for p in s.split("|") if p]
        # typical: sp|P12345|NAME
        if len(parts) >= 2 and _looks_like_acc(parts[1]):
            return parts[1].split(".")[0].upper()
        return parts[0].split(".")[0].upper()
    return s.split()[0].split(".")[0].upper()


def _looks_like_acc(tok: str) -> bool:
    t = tok.upper()
    return len(t) >= 4 and (t[0].isalpha() or t.startswith("A0A"))


def fetch_fasta_sequence(accession: str) -> Tuple[str, int]:
    """Return (one-letter sequence, length) from UniProt FASTA."""
    acc = normalize_uniprot_accession(accession)
    if not acc:
        raise ValueError("Empty UniProt accession")
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
    if resp.status_code != 200:
        raise RuntimeError(f"UniProt FASTA {acc}: HTTP {resp.status_code}")
    seq = _parse_uniprot_fasta_body(resp.text)
    if not seq:
        seq = "".join(c for c in resp.text.upper() if c.isalpha() and c in _VALID_AA)
        if seq:
            return seq, len(seq.replace("U", "C"))
        raise RuntimeError(f"No sequence in UniProt FASTA response for {acc}")
    return seq, len(seq)


def fetch_taxonomy_id(accession: str) -> Optional[int]:
    """NCBI taxonomy id for the UniProt entry (organism.taxonId), if present."""
    acc = normalize_uniprot_accession(accession)
    if not acc:
        return None
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
    if resp.status_code != 200:
        return None
    try:
        data = resp.json()
        org = data.get("organism") or {}
        tid = org.get("taxonId")
        if tid is not None:
            return int(tid)
    except (ValueError, KeyError, TypeError):
        pass
    return None


def fetch_entry_summary(accession: str) -> Dict[str, object]:
    """JSON summary for ranking (taxonomy id)."""
    acc = normalize_uniprot_accession(accession)
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
    if resp.status_code != 200:
        return {}
    try:
        data = resp.json()
        org = data.get("organism") or {}
        tid = org.get("taxonId")
        return {"taxonomy_id": int(tid) if tid is not None else None, "raw": data}
    except (ValueError, KeyError, TypeError):
        return {}
