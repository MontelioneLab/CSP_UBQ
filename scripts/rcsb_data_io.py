"""
Lightweight RCSB Data API access: chain → entity, taxonomy, and UniProt cross-refs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import requests

try:
    from .config import network
except Exception:
    import os as _os, sys as _sys

    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import network

RCSB_REST_PI = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance"
RCSB_REST_PE = "https://data.rcsb.org/rest/v1/core/polymer_entity"


@dataclass
class HoloPolymerContext:
    pdb_id: str
    chain_id: str
    entity_id: Optional[str]
    ncbi_taxonomy_id: Optional[int]
    uniprot_accessions: List[str]
    preferred_uniprot_accession: Optional[str]


def normalize_pdb_id(pdb: str) -> str:
    return (pdb or "").strip().upper()


def polymer_entity_id_for_chain(pdb_id: str, chain_id: str) -> Optional[str]:
    """Resolve PDB entity id (e.g. \"1\") for auth/label chain id (e.g. A)."""
    pid = normalize_pdb_id(pdb_id)
    ch_raw = (chain_id or "").strip()
    if not pid or not ch_raw:
        return None
    ch_path = ch_raw.upper()
    url = f"{RCSB_REST_PI}/{pid}/{ch_path}"
    resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
    if resp.status_code != 200:
        return None
    try:
        data = resp.json()
        identifiers = data.get("polymer_entity_instance") or data
        ident = identifiers.get("rcsb_polymer_entity_instance_container_identifiers") or {}
        eid = ident.get("entity_id")
        return str(eid) if eid is not None else None
    except (KeyError, TypeError, ValueError):
        return None


def fetch_holo_polymer_context(pdb_id: str, chain_id: str) -> HoloPolymerContext:
    """
    UniProt accessions and taxonomy for the holo polymer entity that contains ``chain_id``.
    """
    pid = normalize_pdb_id(pdb_id)
    ch_norm = ((chain_id or "").strip() or "A").upper()
    tax: Optional[int] = None
    eid = polymer_entity_id_for_chain(pid, chain_id)

    if not eid:
        return HoloPolymerContext(pid, ch_norm, None, None, [], None)

    url = f"{RCSB_REST_PE}/{pid}/{eid}"
    resp = requests.get(url, timeout=(network.connect_timeout, network.read_timeout))
    if resp.status_code != 200:
        return HoloPolymerContext(pid, ch_norm, eid, None, [], None)

    preferred: Optional[str] = None
    uni_final: List[str] = []

    try:
        data = resp.json()
        pdata = data.get("polymer_entity") or data
        cid = pdata.get("rcsb_polymer_entity_container_identifiers") or {}
        refs = cid.get("reference_sequence_identifiers") or []
        ranked: List[Tuple[float, str]] = []
        for r in refs:
            if (r.get("database_name") or "").strip() != "UniProt":
                continue
            acc = (r.get("database_accession") or "").strip()
            if not acc:
                continue
            raw_cov = r.get("reference_sequence_coverage") or r.get("entity_sequence_coverage")
            try:
                cov = float(raw_cov) if raw_cov is not None and str(raw_cov) not in ("", ".", "?") else 0.0
            except (TypeError, ValueError):
                cov = 0.0
            base = acc.split(".")[0]
            ranked.append((cov, base))
        ranked.sort(key=lambda x: (-x[0], x[1]))

        for _cov, acc in ranked:
            if acc and acc not in uni_final:
                uni_final.append(acc)
        for r in refs:
            if (r.get("database_name") or "").strip() == "UniProt":
                acc = (r.get("database_accession") or "").strip()
                if acc:
                    base_acc = acc.split(".")[0]
                    if base_acc not in uni_final:
                        uni_final.append(base_acc)
        for u in cid.get("uniprot_ids") or []:
            s = str(u).strip().split(".")[0]
            if s and s not in uni_final:
                uni_final.append(s)

        if ranked:
            preferred = ranked[0][1]
        elif uni_final:
            preferred = uni_final[0]

        orgs = pdata.get("rcsb_entity_source_organism") or []
        if orgs and isinstance(orgs, list):
            t = orgs[0].get("ncbi_taxonomy_id")
            if t is not None:
                tax = int(t)
    except (KeyError, TypeError, ValueError, IndexError):
        uni_final = []
        preferred = None

    return HoloPolymerContext(
        pdb_id=pid,
        chain_id=ch_norm,
        entity_id=eid,
        ncbi_taxonomy_id=tax,
        uniprot_accessions=uni_final,
        preferred_uniprot_accession=preferred,
    )

