"""
RCSB PDB downloader and simple parser utilities.
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple, Optional

import requests

try:
    from .config import network, paths
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import network, paths


class PDBDownloadError(RuntimeError):
    pass


def fetch_pdb(pdb_id: str, cache_dir: Optional[str] = None, force: bool = False) -> str:
    pdb_id = pdb_id.strip()
    if not pdb_id:
        raise ValueError("Empty pdb_id")
    cache_dir = cache_dir or paths.pdb_cache_dir
    os.makedirs(cache_dir, exist_ok=True)
    out_path = os.path.join(cache_dir, f"{pdb_id}.pdb")
    if not force and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PDB] Cache hit {out_path}")
        return out_path
    url = network.rcsb_pdb_url_template.format(pdb_id=pdb_id)
    last_err: Optional[Exception] = None
    for _ in range(network.retries):
        try:
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PDB] GET {url}")
            resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
            if resp.status_code == 200 and resp.content:
                with open(out_path, "wb") as f:
                    f.write(resp.content)
                if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                    print(f"[PDB] OK {pdb_id} ({len(resp.content)} bytes)")
                return out_path
            last_err = PDBDownloadError(f"HTTP {resp.status_code}")
        except Exception as e:
            last_err = e
    raise PDBDownloadError(f"Failed to download PDB {pdb_id}: {last_err}")


def parse_pdb_sequences(pdb_path: str) -> Dict[str, str]:
    """Extract per-chain sequences (one-letter) from ATOM records.

    Uses the standard residue name to one-letter map. Non-standard residues are "X".
    """
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    chain_to_seq: Dict[str, List[Tuple[int, str]]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue
            resname = line[17:20].strip().upper()
            chain_id = (line[21] or " ").strip() or "A"
            resseq = int(line[22:26])
            aa = three_to_one.get(resname, "X")
            chain_to_seq.setdefault(chain_id, []).append((resseq, aa))
    chain_sequences: Dict[str, str] = {}
    for ch, items in chain_to_seq.items():
        # keep first occurrence per residue number
        seen = set()
        seq_list: List[str] = []
        for resseq, aa in items:
            if resseq in seen:
                continue
            seen.add(resseq)
            seq_list.append(aa)
        chain_sequences[ch] = "".join(seq_list)
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PDB] Parsed chains: {', '.join(sorted(chain_sequences.keys())) or 'none'}")
    return chain_sequences



