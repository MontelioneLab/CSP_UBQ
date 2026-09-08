"""Secondary-structure helpers (DSSP H/E/C labels)."""

from __future__ import annotations

import os
import tempfile
from typing import Dict, Optional

# DSSP 8-state codes collapsed to helix / sheet / coil-loop.
DSSP_REDUCE = {
    "H": "H",
    "G": "H",
    "I": "H",
    "B": "E",
    "E": "E",
    "T": "C",
    "S": "C",
    "-": "C",
    "C": "C",
}

SS_NAME = {"H": "helix", "E": "sheet", "C": "loop"}


def _write_heavy_atom_pdb(pdb_path: str, out_path: str) -> bool:
    """
    Write model-0 heavy atoms only to ``out_path``.

    mkdssp 4.x often segfaults on NMR ensembles with hydrogens; a stripped
    single-model heavy-atom PDB is much more reliable.
    """
    try:
        from Bio.PDB import PDBIO, PDBParser, Select  # type: ignore
    except Exception:
        return False

    class _HeavySelect(Select):  # type: ignore[misc]
        def accept_model(self, model):  # noqa: ANN001
            return model.id == 0

        def accept_atom(self, atom):  # noqa: ANN001
            try:
                if atom.element and str(atom.element).upper() == "H":
                    return False
            except Exception:
                pass
            name = atom.get_name().strip().upper()
            return not name.startswith("H")

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ss_model", pdb_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_path, _HeavySelect())
        return os.path.isfile(out_path) and os.path.getsize(out_path) > 0
    except Exception:
        return False


def _dssp_map_from_pdb(pdb_path: str, receptor_chain: Optional[str]) -> Dict[int, str]:
    """Run BioPython DSSP on ``pdb_path`` and return {resnum -> H/E/C}."""
    ss_map: Dict[int, str] = {}
    try:
        from Bio.PDB import PDBParser  # type: ignore
        from Bio.PDB.DSSP import DSSP  # type: ignore
    except Exception:
        return ss_map

    chain_filter = (receptor_chain or "").strip()
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("ss_model", pdb_path)
        model = structure[0]
        dssp = DSSP(model, pdb_path, dssp="mkdssp")
        for key in dssp.keys():
            chain_id = key[0]
            res_id = key[1]
            if chain_filter and chain_id != chain_filter:
                continue
            try:
                resseq = int(res_id[1])
            except Exception:
                continue
            ss_raw = str(dssp[key][2]).strip() if dssp[key] is not None else "C"
            ss_map[resseq] = DSSP_REDUCE.get(ss_raw, "C")
    except Exception:
        return {}
    return ss_map


def compute_dssp_secondary_structure(
    pdb_path: str, receptor_chain: Optional[str] = None
) -> Dict[int, str]:
    """
    Compute per-residue secondary structure with DSSP.

    Returns ``{residue_number -> 'H'|'E'|'C'}`` for the requested chain
    (or all chains when ``receptor_chain`` is empty/None). Missing or failed
    DSSP runs return an empty dict.

    Always runs mkdssp on a temporary model-0 heavy-atom rewrite so NMR
    structures with hydrogens do not crash mkdssp 4.x.
    """
    if not pdb_path or not os.path.exists(pdb_path):
        return {}

    try:
        fd, tmp_path = tempfile.mkstemp(suffix=".pdb", prefix="dssp_heavy_")
        os.close(fd)
    except Exception:
        return {}

    try:
        if not _write_heavy_atom_pdb(pdb_path, tmp_path):
            return {}
        ss_map = _dssp_map_from_pdb(tmp_path, receptor_chain)
        if ss_map:
            return ss_map
        # Chain ID mismatch: retry without chain filter.
        if (receptor_chain or "").strip():
            return _dssp_map_from_pdb(tmp_path, None)
        return {}
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
