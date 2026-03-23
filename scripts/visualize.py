"""
Visualization utilities: PyMOL coloring script and simple plots.
"""

from __future__ import annotations

import os
import csv

from typing import List, Optional, Dict, Set, Sequence, Tuple

try:
    import matplotlib
    # Set non-interactive backend to avoid GUI issues in threaded environments
    matplotlib.use('Agg')  # Use Agg backend (no GUI)
    import matplotlib.pyplot as plt
    _HAS_PLT = True
except Exception:
    _HAS_PLT = False

try:
    from .csp import CSPResult, compute_atom_deltas_with_offset
    from .config import classification_colors, hex_to_rgb01 as _config_hex_to_rgb01
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.csp import CSPResult, compute_atom_deltas_with_offset
    from scripts.config import classification_colors, hex_to_rgb01 as _config_hex_to_rgb01


AA_THREE_TO_ONE: Dict[str, str] = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def _three_letter_to_one(residue_name: str) -> str:
    return AA_THREE_TO_ONE.get(residue_name.upper(), 'X')


def _get_classification_colors() -> Dict[str, str]:
    return {
        "TP": classification_colors.TP,
        "FP": classification_colors.FP,
        "TN": classification_colors.TN,
        "FN": classification_colors.FN,
    }


def _load_metrics_from_confusion_csv(out_png: str) -> Optional[Dict[str, float]]:
    """Load precision/recall/F1/MCC for this system from outputs/confusion_matrix_per_system.csv."""
    system_id = os.path.basename(os.path.dirname(os.path.abspath(out_png))).strip()
    if not system_id:
        return None

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv_path = os.path.join(project_root, "outputs", "confusion_matrix_per_system.csv")
    if not os.path.exists(csv_path):
        return None

    try:
        with open(csv_path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row_system_id = (row.get("system_id") or "").strip()
                if row_system_id.lower() != system_id.lower():
                    continue
                return {
                    "precision": float(row.get("precision") or 0.0),
                    "recall": float(row.get("recall") or 0.0),
                    "f1_score": float(row.get("f1_score") or 0.0),
                    "mcc": float(row.get("mcc") or 0.0),
                }
    except Exception:
        return None

    return None


def create_three_row_residue_ticks(residue_indices):
    """
    Create three rows of residue index ticks:
    - Row 1: One's place (always shown)
    - Row 2: Ten's place (only when it changes)
    - Row 3: Hundred's place (only when it changes)
    
    Args:
        residue_indices: List of residue index numbers
        
    Returns:
        tuple: (ones_ticks, tens_ticks, hundreds_ticks) where each is a list of (position, label) tuples
    """
    ones_ticks = []
    tens_ticks = []
    hundreds_ticks = []
    
    # Track previous values to detect changes
    prev_tens = None
    prev_hundreds = None
    
    for i, residue_idx in enumerate(residue_indices):
        # Always add ones place - use actual residue index as position
        ones_ticks.append((residue_idx, str(residue_idx % 10)))
        
        # Add tens place only when it changes
        tens_digit = (residue_idx // 10) % 10
        if tens_digit != prev_tens and tens_digit > 0:
            tens_ticks.append((residue_idx, str(tens_digit)))
            prev_tens = tens_digit
        
        # Add hundreds place only when it changes
        hundreds_digit = (residue_idx // 100) % 10
        if hundreds_digit != prev_hundreds and hundreds_digit > 0:
            hundreds_ticks.append((residue_idx, str(hundreds_digit)))
            prev_hundreds = hundreds_digit
    
    return ones_ticks, tens_ticks, hundreds_ticks


def _infer_binding_dataset(residue_info: Sequence[Dict]) -> str:
    """
    Detect which analysis produced the residue info.

    Returns:
        'sasa' for occlusion results,
        'interaction' for interaction filter results,
        'distance' for CA distance filter,
        'union' for combined interaction/SASA/CA annotations,
        'unknown' otherwise.
    """
    for info in residue_info:
        if (
            'has_sasa_occlusion' in info
            or 'has_ca_distance' in info
        ):
            return 'union'
        if 'is_occluded' in info:
            return 'sasa'
        if (
            'interaction_category' in info
            or 'has_hbond' in info
            or 'has_charge_complement' in info
            or 'has_pi_contact' in info
        ):
            return 'interaction'
        if 'passes_filter' in info:
            return 'distance'
    return 'unknown'


def _build_binding_lookup(binding_results: Optional[Dict]) -> Tuple[Dict[int, bool], str]:
    """
    Create a lookup {residue_number -> bool} indicating binding/contact residues.
    Returns the lookup and the detected dataset type.
    """
    lookup: Dict[int, bool] = {}
    if not binding_results:
        return lookup, 'unknown'

    residue_info = binding_results.get('residue_info') or []
    dataset_type = binding_results.get('dataset_type') or _infer_binding_dataset(residue_info)

    for info in residue_info:
        res_num = info.get('residue_number')
        if res_num is None:
            continue
        if dataset_type == 'sasa':
            is_binding = bool(info.get('is_occluded'))
        elif dataset_type == 'union':
            is_binding = bool(
                info.get('has_hbond')
                or info.get('has_charge_complement')
                or info.get('has_pi_contact')
                or info.get('has_ca_distance')
                or info.get('has_sasa_occlusion')
                or info.get('has_any_atom_sub_2A')
            )
        elif dataset_type == 'interaction':
            is_binding = bool(
                info.get('has_hbond')
                or info.get('has_charge_complement')
                or info.get('has_pi_contact')
            )
        elif dataset_type == 'distance':
            is_binding = bool(info.get('passes_filter'))
        else:
            is_binding = bool(info.get('is_binding'))
        lookup[int(res_num)] = is_binding

    return lookup, dataset_type


def _parse_pdb_secondary_structure_records(pdb_path: str, receptor_chain: Optional[str]) -> Dict[int, str]:
    """
    Parse HELIX/SHEET records from a PDB file.
    Returns a map {residue_number -> 'H'|'E'} for the requested chain.
    """
    ss_map: Dict[int, str] = {}
    if not pdb_path or not os.path.exists(pdb_path):
        return ss_map

    chain_filter = (receptor_chain or "").strip()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("HELIX"):
                    start_chain = line[19].strip() if len(line) > 19 else ""
                    end_chain = line[31].strip() if len(line) > 31 else ""
                    if chain_filter and start_chain != chain_filter and end_chain != chain_filter:
                        continue
                    start_res = line[21:25].strip()
                    end_res = line[33:37].strip()
                    if not start_res or not end_res:
                        continue
                    try:
                        start_i = int(start_res)
                        end_i = int(end_res)
                    except ValueError:
                        continue
                    for resi in range(min(start_i, end_i), max(start_i, end_i) + 1):
                        ss_map[resi] = "H"
                elif line.startswith("SHEET"):
                    start_chain = line[21].strip() if len(line) > 21 else ""
                    end_chain = line[32].strip() if len(line) > 32 else ""
                    if chain_filter and start_chain != chain_filter and end_chain != chain_filter:
                        continue
                    start_res = line[22:26].strip()
                    end_res = line[33:37].strip()
                    if not start_res or not end_res:
                        continue
                    try:
                        start_i = int(start_res)
                        end_i = int(end_res)
                    except ValueError:
                        continue
                    for resi in range(min(start_i, end_i), max(start_i, end_i) + 1):
                        # Keep helix labels if present; otherwise annotate as strand.
                        ss_map[resi] = ss_map.get(resi, "E")
    except Exception:
        return {}
    return ss_map


def _compute_dssp_secondary_structure(pdb_path: str, receptor_chain: Optional[str]) -> Dict[int, str]:
    """
    Compute residue secondary structure with DSSP and return {residue_number -> H/E/C}.
    """
    ss_map: Dict[int, str] = {}
    if not pdb_path or not os.path.exists(pdb_path):
        return ss_map
    try:
        from Bio.PDB import PDBParser  # type: ignore
        from Bio.PDB.DSSP import DSSP  # type: ignore
    except Exception:
        return ss_map

    reduce_ss = {"H": "H", "G": "H", "I": "H", "B": "E", "E": "E", "T": "C", "S": "C", "-": "C", "C": "C"}
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
            ss_map[resseq] = reduce_ss.get(ss_raw, "C")
    except Exception:
        return {}
    return ss_map


def _resolve_secondary_structure_by_holo_index(
    results_hn: List[CSPResult],
    binding_results: Optional[Dict],
    residue_indices: List[int],
    *,
    pdb_id: Optional[str] = None,
    structure_pdb_path: Optional[str] = None,
    receptor_chain: Optional[str] = None,
) -> List[str]:
    """
    Resolve per-holo-index secondary structure labels (H/E/C) for the panel x-axis.
    """
    if not residue_indices:
        return []

    chain_id = (receptor_chain or (binding_results or {}).get("receptor_chain") or "").strip()
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates: List[str] = []

    if structure_pdb_path:
        candidates.append(structure_pdb_path)
    if pdb_id:
        candidates.append(os.path.join(project_root, "PDB_FILES", f"{str(pdb_id).lower()}.pdb"))
    # Last-resort guess for pipeline-generated PDB near output image.
    out_png = None
    if isinstance(binding_results, dict):
        out_png = binding_results.get("_plot_out_png")
    if out_png:
        out_dir = os.path.dirname(os.path.abspath(str(out_png)))
        if pdb_id:
            candidates.append(os.path.join(out_dir, f"{str(pdb_id).lower()}_csp.pdb"))

    pdb_path = next((c for c in candidates if c and os.path.exists(c)), None)
    if not pdb_path:
        # If no structure file, return all-coil so subplot still renders deterministically.
        return ["C"] * len(residue_indices)

    ss_by_resnum = _parse_pdb_secondary_structure_records(pdb_path, chain_id)
    if not ss_by_resnum:
        ss_by_resnum = _compute_dssp_secondary_structure(pdb_path, chain_id)

    position_map = create_sequence_alignment_map_from_results(results_hn, binding_results or {})
    labels: List[str] = []
    for holo_index in residue_indices:
        pdb_resnum = position_map.get(holo_index)
        if pdb_resnum is None:
            labels.append("C")
        else:
            labels.append(ss_by_resnum.get(int(pdb_resnum), "C"))
    return labels


def _draw_secondary_structure_track(
    ax,
    residue_indices: List[int],
    ss_labels: List[str],
    *,
    inside_axis: bool = False,
) -> None:
    """
    Draw a compact secondary-structure track aligned with residue x-positions.
    """
    if not residue_indices or not ss_labels:
        return

    from matplotlib import patches as mpatches  # type: ignore

    reduced: List[str] = []
    for label in ss_labels:
        s = (label or "C").strip().upper()
        if s in ("H", "G", "I"):
            reduced.append("H")
        elif s in ("E", "B"):
            reduced.append("E")
        else:
            reduced.append("C")

    # Fill any index gaps as loops so the SS track is continuous.
    index_to_ss = {idx: ss for idx, ss in zip(residue_indices, reduced)}
    full_indices = list(range(min(residue_indices), max(residue_indices) + 1))
    reduced_full = [index_to_ss.get(idx, "C") for idx in full_indices]

    blocks: List[List[object]] = []
    prev_ss: Optional[str] = None
    for x, ss in zip(full_indices, reduced_full):
        if prev_ss is None or ss != prev_ss:
            blocks.append([ss, x, x])
        else:
            blocks[-1][2] = x
        prev_ss = ss

    if inside_axis:
        y_min, y_max = ax.get_ylim()
        if y_max <= y_min:
            return
        y_span = y_max - y_min
        width = max(0.06, y_span * 0.075)
        band_bottom = y_min - (width * 1.7)
        band_top = band_bottom + width
        ax.set_ylim(band_bottom - (width * 0.2), y_max)
        y_mid = (band_bottom + band_top) / 2.0
        zorder = 5
    else:
        width = 1.0
        y_mid = width / 2.0
        zorder = 3

    fc_helix = "firebrick"
    fc_sheet = "#F0E442"  # Yellow
    fc_coil = "#0072B2"   # Blue (loops / missing SS)
    ec = "none"
    helix_arc_width = 3.5
    coil_thickness = width / 4.0
    sheet_thickness = width * (2.0 / 3.0)

    for blk_idx, block in enumerate(blocks):
        ss_type = str(block[0])
        start_x = int(block[1])
        end_x = int(block[2])
        length = end_x - start_x + 1
        start = start_x - 0.5

        if ss_type == "H":
            arc_len = 0.5
            st_theta = 0
            en_theta = 180
            t = start
            while t <= end_x + 0.5:
                origin = (t + 0.25, y_mid)
                arc = mpatches.Arc(
                    origin,
                    arc_len,
                    width,
                    linewidth=helix_arc_width,
                    theta1=st_theta - 1,
                    theta2=en_theta + 1,
                    edgecolor=fc_helix,
                    zorder=zorder,
                )
                ax.add_patch(arc)
                st_theta += 180
                en_theta += 180
                t += arc_len
        elif ss_type == "E":
            arrow = mpatches.FancyArrow(
                start,
                y_mid,
                float(length),
                0.0,
                length_includes_head=True,
                head_length=min(float(length) / 4.0, 1.0),
                head_width=width - 0.001,
                width=sheet_thickness,
                facecolor=fc_sheet,
                edgecolor=ec,
                linewidth=1.0,
                zorder=zorder,
            )
            ax.add_patch(arrow)
        else:
            seg_start = start
            seg_len = float(length)
            if blk_idx > 0 and str(blocks[blk_idx - 1][0]) in ("H", "E"):
                seg_start -= 0.1
                seg_len += 0.1
            if blk_idx + 1 < len(blocks) and str(blocks[blk_idx + 1][0]) in ("H", "E"):
                seg_len += 0.1
            rect = mpatches.Rectangle(
                (seg_start, y_mid - coil_thickness / 2.0),
                seg_len,
                coil_thickness,
                linewidth=1.0,
                edgecolor=ec,
                facecolor=fc_coil,
                zorder=zorder,
            )
            ax.add_patch(rect)

    if not inside_axis:
        ax.set_ylim([-0.05, 1.15])
        ax.set_xlim(min(residue_indices) - 0.5, max(residue_indices) + 0.5)
        ax.set_aspect(0.5)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)


def _extract_csp_sequence(results: Sequence[CSPResult]) -> Tuple[str, List[int]]:
    sequence_chars: List[str] = []
    positions: List[int] = []
    for r in results:
        aa = (r.holo_aa or '').strip().upper()
        if not aa or aa == 'P':
            continue
        sequence_chars.append(aa)
        positions.append(r.holo_index)
    return ''.join(sequence_chars), positions


def _extract_binding_sequence(binding_results: Optional[Dict], receptor_chain: Optional[str]) -> Tuple[str, List[int]]:
    residue_info = (binding_results or {}).get('residue_info', [])
    sequence_chars: List[str] = []
    residue_numbers: List[int] = []
    for info in residue_info:
        chain_id = (info.get('chain_id') or '').strip()
        if chain_id and receptor_chain and chain_id != receptor_chain:
            continue
        res_name = (info.get('residue_name') or '').strip()
        res_num = info.get('residue_number')
        if not res_name or res_num is None:
            continue
        letter = _three_letter_to_one(res_name)
        if not letter:
            continue
        sequence_chars.append(letter)
        try:
            residue_numbers.append(int(res_num))
        except Exception:
            continue
    return ''.join(sequence_chars), residue_numbers


def _extract_pdb_sequence(pdb_path: str, receptor_chain: Optional[str]) -> Tuple[str, List[int]]:
    if not pdb_path or not os.path.exists(pdb_path):
        return '', []
    sequence_chars: List[str] = []
    residue_numbers: List[int] = []
    seen_keys: Set[Tuple[str, int, str]] = set()
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if not line.startswith('ATOM'):
                continue
            chain_id = (line[21].strip() or '')
            if receptor_chain and chain_id and chain_id != receptor_chain:
                continue
            res_name = line[17:20].strip()
            try:
                res_num = int(line[22:26])
            except ValueError:
                continue
            insertion_code = line[26].strip()
            key = (chain_id, res_num, insertion_code)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            letter = _three_letter_to_one(res_name)
            if not letter:
                continue
            sequence_chars.append(letter)
            residue_numbers.append(res_num)
    return ''.join(sequence_chars), residue_numbers


def _build_index_map(seq_a: str, seq_b: str, len_a: int, len_b: int) -> Dict[int, int]:
    if not seq_a or not seq_b:
        return {}
    try:
        from .align import align_global
    except Exception:
        import os as _os, sys as _sys
        _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
        from scripts.align import align_global
    _, _, mapping, _ = align_global(seq_a, seq_b)
    index_map: Dict[int, int] = {}
    for pos_a, pos_b in mapping:
        idx_a = pos_a - 1
        idx_b = pos_b - 1
        if 0 <= idx_a < len_a and 0 <= idx_b < len_b:
            index_map[idx_a] = idx_b
    return index_map


def _create_csp_binding_pdb_alignment(
    results: Sequence[CSPResult],
    binding_results: Optional[Dict],
    pdb_path: str,
    receptor_chain: Optional[str],
) -> Dict[int, int]:
    csp_seq, csp_positions = _extract_csp_sequence(results)
    binding_seq, binding_positions = _extract_binding_sequence(binding_results, receptor_chain)
    pdb_seq, pdb_positions = _extract_pdb_sequence(pdb_path, receptor_chain)

    if not csp_seq or not binding_seq:
        return create_sequence_alignment_map_from_results(results, binding_results or {})

    csp_to_binding = _build_index_map(csp_seq, binding_seq, len(csp_positions), len(binding_positions))
    if not csp_to_binding:
        return create_sequence_alignment_map_from_results(results, binding_results or {})

    binding_to_pdb: Dict[int, int] = {}
    if pdb_seq:
        binding_to_pdb = _build_index_map(binding_seq, pdb_seq, len(binding_positions), len(pdb_positions))

    position_map: Dict[int, int] = {}
    for csp_idx, binding_idx in csp_to_binding.items():
        if not (0 <= csp_idx < len(csp_positions)) or not (0 <= binding_idx < len(binding_positions)):
            continue
        pdb_residue_number = None
        if binding_to_pdb:
            pdb_idx = binding_to_pdb.get(binding_idx)
            if pdb_idx is not None and 0 <= pdb_idx < len(pdb_positions):
                pdb_residue_number = pdb_positions[pdb_idx]
        residue_number = pdb_residue_number if pdb_residue_number is not None else binding_positions[binding_idx]
        if residue_number is not None:
            position_map[csp_positions[csp_idx]] = residue_number

    if not position_map:
        return create_sequence_alignment_map_from_results(results, binding_results or {})
    return position_map


def write_pymol_color_script(
    results: List[CSPResult], 
    pdb_id: str, 
    out_path: str,
    sasa_results: Optional[dict] = None,
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Minimal script that loads the pdb from current directory; users may open it in PyMOL.
    Colors the peptide chain (smaller chain) cyan.
    
    Args:
        results: List of CSPResult objects (not used but kept for compatibility)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        sasa_results: Optional results dictionary from compute_sasa_occlusion (for chain info)
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from sasa_results)
        output_dir: Output directory path (used to construct PDB file paths)
    """
    import os
    
    # Get chains from parameters, fallback to sasa_results if not provided
    if receptor_chain is None and sasa_results:
        receptor_chain = sasa_results.get('receptor_chain')
    if ligand_chain is None and sasa_results:
        ligand_chain = sasa_results.get('ligand_chain')
    
    # Construct PDB path relative to project root
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        # Ensure it ends with / and use string concatenation (os.path.join doesn't handle ./ well)
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_csp.pdb"
        # Convert relative path to absolute for file operations
        if pdb_path.startswith('./'):
            # Get project root: go up from out_path until we find the directory containing 'outputs'
            current = os.path.dirname(os.path.abspath(out_path))
            while current != os.path.dirname(current):  # Stop at filesystem root
                if os.path.basename(current) == 'outputs':
                    project_root = os.path.dirname(current)
                    break
                current = os.path.dirname(current)
            else:
                # Fallback: assume project root is 2 levels up from out_path
                project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
            pdb_path_abs = os.path.join(project_root, pdb_path[2:])
        else:
            pdb_path_abs = pdb_path
    else:
        # Fallback to old behavior
        pdb_path = f"./outputs/{pdb_id}/{pdb_id}_csp.pdb"
        # Get project root: go up from out_path until we find the directory containing 'outputs'
        current = os.path.dirname(os.path.abspath(out_path))
        while current != os.path.dirname(current):  # Stop at filesystem root
            if os.path.basename(current) == 'outputs':
                project_root = os.path.dirname(current)
                break
            current = os.path.dirname(current)
        else:
            # Fallback: assume project root is 2 levels up from out_path
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_path_abs = os.path.join(project_root, pdb_path[2:])
    
    # Create the PDB file if it doesn't exist
    if results and not os.path.exists(pdb_path_abs):
        create_pdb_with_csp_bfactors(results, pdb_id, pdb_path_abs, sasa_results=sasa_results)
    
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {pdb_path}, {pdb_id}_structure")
    lines.append(f"hide everything, {pdb_id}_structure")
    lines.append(f"show cartoon, {pdb_id}_structure")
    lines.append(f"spectrum b, blue_white_red, {pdb_id}_structure")
    
    # Color ligand chain cyan
    if ligand_chain:
        lines.append(f"color cyan, {pdb_id}_structure and chain {ligand_chain}")
    
    lines.append(f"set cartoon_transparency, 0.2, {pdb_id}_structure")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_pymol_color_csp_mask_script(
    results: List[CSPResult],
    pdb_id: str,
    out_path: str,
    sasa_results: Optional[dict] = None,
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    significance_field: str = 'significant',
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate PyMOL script with three-color masking:
    - Peptide chain: cyan
    - Protein residues with significant CSPs: red
    - Protein residues with insignificant or no CSPs: gray30
    
    Args:
        results: List of CSPResult objects from CSP analysis
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        sasa_results: Optional results dictionary from compute_sasa_occlusion (for chain info and alignment)
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from sasa_results)
        significance_field: CSPResult attribute used to determine binary significance (default: 'significant')
        output_dir: Output directory path (used to construct PDB file paths)
    """
    import os
    # Get chains from parameters, fallback to sasa_results if not provided
    if receptor_chain is None and sasa_results:
        receptor_chain = sasa_results.get('receptor_chain', 'A')
    if receptor_chain is None:
        receptor_chain = 'A'
    if ligand_chain is None and sasa_results:
        ligand_chain = sasa_results.get('ligand_chain', 'B')
    if ligand_chain is None:
        ligand_chain = 'B'
    
    # Get significant residue numbers using sequence alignment
    significant_residues: Set[int] = set()
    if sasa_results:
        significant_residues = _get_significant_residue_numbers(results, sasa_results, significance_field)
    
    # Construct PDB path relative to project root
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        # Ensure it ends with / and use string concatenation (os.path.join doesn't handle ./ well)
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_csp.pdb"
        # Convert relative path to absolute for file operations
        if pdb_path.startswith('./'):
            # Get project root: go up from out_path until we find the directory containing 'outputs'
            current = os.path.dirname(os.path.abspath(out_path))
            while current != os.path.dirname(current):  # Stop at filesystem root
                if os.path.basename(current) == 'outputs':
                    project_root = os.path.dirname(current)
                    break
                current = os.path.dirname(current)
            else:
                # Fallback: assume project root is 2 levels up from out_path
                project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
            pdb_path_abs = os.path.join(project_root, pdb_path[2:])
        else:
            pdb_path_abs = pdb_path
    else:
        pdb_path = f"./outputs/{pdb_id}/{pdb_id}_csp.pdb"
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_path_abs = os.path.join(project_root, pdb_path[2:])
    
    # Create the PDB file if it doesn't exist
    if results and not os.path.exists(pdb_path_abs):
        create_pdb_with_csp_bfactors(results, pdb_id, pdb_path_abs, sasa_results=sasa_results)
    
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {pdb_path}, {pdb_id}_structure")
    lines.append(f"hide everything, {pdb_id}_structure")
    lines.append(f"show cartoon, {pdb_id}_structure")
    
    # Color peptide chain cyan
    lines.append(f"color cyan, {pdb_id}_structure and chain {ligand_chain}")
    
    # Color all protein chain residues gray30 (default)
    lines.append(f"color gray30, {pdb_id}_structure and chain {receptor_chain}")
    
    # Override with red for significant CSP residues
    if significant_residues:
        for res_num in sorted(significant_residues):
            lines.append(f"color red, {pdb_id}_structure and chain {receptor_chain} and resi {res_num}")
    
    lines.append(f"set cartoon_transparency, 0.2, {pdb_id}_structure")
    
    # Add comments for clarity
    lines.append(f"# Color scheme:")
    lines.append(f"# Ligand chain ({ligand_chain}): cyan")
    lines.append(f"# Receptor chain ({receptor_chain}): gray30 (non-significant/no CSP), red (significant CSP via {significance_field})")
    lines.append(f"# Significant residues: {len(significant_residues)}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def create_pdb_with_delta_sasa_bfactors(sasa_results: dict, pdb_id: str, out_path: str) -> str:
    """
    Create a new PDB file with B-factor column updated with delta SASA values.
    Uses medoid model selection for multi-model PDB files.
    
    Args:
        sasa_results: Results dictionary from compute_sasa_occlusion
        pdb_id: PDB identifier for the structure
        out_path: Path to save the modified PDB file
        
    Returns:
        Path to the created PDB file
    """
    import os
    import tempfile
    import shutil
    
    # Import the medoid selection function from sasa_analysis
    try:
        from .sasa_analysis import find_medoid_model_from_pdb
    except Exception:
        import os as _os, sys as _sys
        _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
        from scripts.sasa_analysis import find_medoid_model_from_pdb
    
    # Input PDB file path
    input_pdb_path = f"./PDB_FILES/{pdb_id}.pdb"
    
    if not os.path.exists(input_pdb_path):
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb_path}")
    
    # Find medoid model if it's a multi-model PDB
    medoid_pdb_path = find_medoid_model_from_pdb(input_pdb_path)
    print(f"[DEBUG] Using PDB file: {medoid_pdb_path}")
    
    # Get chain information
    receptor_chain = sasa_results.get('receptor_chain', 'A')
    ligand_chain = sasa_results.get('ligand_chain', 'B')
    residue_info = sasa_results.get('residue_info', [])
    
    print(f"[DEBUG] Creating PDB with delta SASA B-factors:")
    print(f"[DEBUG] Input PDB: {input_pdb_path}")
    print(f"[DEBUG] Medoid PDB: {medoid_pdb_path}")
    print(f"[DEBUG] Output PDB: {out_path}")
    print(f"[DEBUG] Receptor chain: {receptor_chain}")
    print(f"[DEBUG] Processing {len(residue_info)} residues")
    
    # Create dictionary for quick lookup of delta SASA values
    delta_sasa_lookup = {}
    for info in residue_info:
        res_num = info['residue_number']
        delta_sasa = info['delta_sasa']
        delta_sasa_lookup[res_num] = delta_sasa
        if len(delta_sasa_lookup) <= 10:  # Debug first 10
            print(f"[DEBUG] Residue {res_num}: delta_sasa = {delta_sasa:.3f} Å²")
    
    # Read medoid PDB and modify B-factors
    modified_lines = []
    atoms_modified = 0
    
    with open(medoid_pdb_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract chain ID and residue number
                chain_id = line[21:22].strip()
                res_num = int(line[22:26].strip())
                
                # Check if this atom belongs to the protein chain and has delta SASA data
                if chain_id == receptor_chain and res_num in delta_sasa_lookup:
                    # Update B-factor column (columns 60-66)
                    delta_sasa = delta_sasa_lookup[res_num]
                    # Format B-factor to 6.2f (6 characters, 2 decimal places)
                    bfactor_str = f"{delta_sasa:6.2f}"
                    
                    # Replace B-factor in the line
                    modified_line = line[:60] + bfactor_str + line[66:]
                    modified_lines.append(modified_line)
                    atoms_modified += 1
                    
                    if atoms_modified <= 5:  # Debug first 5 atoms
                        atom_name = line[12:16].strip()
                        print(f"[DEBUG] Modified atom {atom_name} in residue {res_num}: B-factor = {delta_sasa:.2f}")
                else:
                    # Keep original line for non-protein atoms or residues without data
                    modified_lines.append(line)
            else:
                # Keep non-ATOM lines as-is
                modified_lines.append(line)
    
    # Write modified PDB file
    with open(out_path, 'w') as f:
        f.writelines(modified_lines)
    
    print(f"[DEBUG] PDB file created: {out_path}")
    print(f"[DEBUG] Total atoms modified: {atoms_modified}")
    print(f"[DEBUG] Delta SASA range: {min(delta_sasa_lookup.values()):.2f} to {max(delta_sasa_lookup.values()):.2f} Å²")
    
    # Clean up temporary files if we created them
    if medoid_pdb_path != input_pdb_path:
        try:
            temp_dir = os.path.dirname(medoid_pdb_path)
            shutil.rmtree(temp_dir)
            print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
        except Exception as e:
            print(f"[DEBUG] Warning: Could not clean up temporary files: {e}")
    
    return out_path


def write_pymol_delta_sasa_script(
    sasa_results: dict, 
    pdb_id: str, 
    out_path: str,
    binary_mode: bool = False,
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate PyMOL script to visualize delta SASA values as backbone coloring.
    Creates a modified PDB file with B-factors updated and loads that instead.
    
    Args:
        sasa_results: Results dictionary from compute_sasa_occlusion
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        binary_mode: If True, color occluded residues red and others blue instead of a heatmap.
    """
    import os
    
    # Get chain information
    receptor_chain = sasa_results.get('receptor_chain', 'A')
    ligand_chain = sasa_results.get('ligand_chain', 'B')
    residue_info = sasa_results.get('residue_info', [])
    
    print(f"[DEBUG] PyMOL Delta SASA Heatmap Generation:")
    print(f"[DEBUG] Receptor chain: {receptor_chain}")
    print(f"[DEBUG] Ligand chain: {ligand_chain}")
    print(f"[DEBUG] Number of residues in residue_info: {len(residue_info)}")
    
    # Create output directory if it doesn't exist
    script_output_dir = os.path.dirname(out_path)
    if script_output_dir:
        os.makedirs(script_output_dir, exist_ok=True)
    
    # Use provided output_dir or derive from script location
    pdb_output_dir = output_dir or script_output_dir
    # If output_dir is provided, it should be relative to project root (e.g., "./outputs/2mur_1/")
    # For file operations, we need the actual directory path
    if output_dir and output_dir.startswith('./'):
        # Convert relative path to absolute for file operations
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_output_dir_abs = os.path.join(project_root, output_dir[2:].rstrip('/'))
    else:
        pdb_output_dir_abs = pdb_output_dir
    
    # Create modified PDB file with delta SASA B-factors
    modified_pdb_path = os.path.join(pdb_output_dir_abs, f"{pdb_id}_delta_sasa.pdb")
    create_pdb_with_delta_sasa_bfactors(sasa_results, pdb_id, modified_pdb_path)
    
    # Generate PyMOL script - use path relative to project root
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_delta_sasa.pdb"
        if not pdb_path.startswith('./') and not pdb_path.startswith('/'):
            pdb_path = './' + pdb_path
    else:
        pdb_path = modified_pdb_path
    
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {pdb_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    # Color peptide chain in green
    lines.append(f"color green, structure and chain {ligand_chain}")
    
    # Calculate delta SASA range for spectrum
    delta_sasa_values = [info['delta_sasa'] for info in residue_info]
    occluded_residues = [
        info['residue_number']
        for info in residue_info
        if info.get('is_occluded') and info.get('residue_number') is not None
    ]
    
    if binary_mode:
        print(f"[DEBUG] Applying binary SASA coloring (occluded -> red)")
        lines.append("# Binary SASA coloring: red=occluded, blue=non-occluded")
        lines.append(f"color blue, structure and chain {receptor_chain}")
        for res_num in sorted(occluded_residues):
            lines.append(f"color red, structure and chain {receptor_chain} and resi {res_num}")
        lines.append(f"# Occluded residues colored red: {len(occluded_residues)}")
    else:
        if delta_sasa_values:
            min_delta = min(delta_sasa_values)
            max_delta = max(delta_sasa_values)
            print(f"[DEBUG] Delta SASA range: {min_delta:.2f} to {max_delta:.2f} Å²")
            print(f"[DEBUG] Applying spectrum coloring with min={min_delta:.2f}, max={max_delta:.2f}")
            
            lines.append(f"# Delta SASA range: {min_delta:.2f} to {max_delta:.2f} Å²")
            
            # Apply heatmap coloring: Blue (low) -> White -> Red (high delta SASA)
            lines.append(f"spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum={min_delta:.2f}, maximum={max_delta:.2f}")
        else:
            print("[DEBUG] No delta SASA values found, using auto-scaling")
            lines.append(f"spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum=0, maximum=auto")
    
    # Set transparency and visual enhancements
    lines.append("set cartoon_transparency, 0.2, structure")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_ring_mode, 1")
    
    # Add colorbar/legend information
    lines.append("# Colorbar information:")
    lines.append("# Blue regions: Low delta SASA (minimal occlusion)")
    lines.append("# White regions: Medium delta SASA (moderate occlusion)")
    lines.append("# Red regions: High delta SASA (strong occlusion)")
    
    # Add title and summary
    n_occluded = sasa_results.get('n_occluded', 0)
    fraction_occluded = sasa_results.get('fraction_occluded', 0.0)
    threshold = sasa_results.get('sasa_threshold', 5.0)
    avg_percent_burial = sasa_results.get('avg_percent_burial', 0.0)
    
    lines.append(f"# Delta SASA Heatmap Analysis Summary:")
    lines.append(f"# Receptor chain: {receptor_chain}")
    lines.append(f"# Ligand chain: {ligand_chain}")
    lines.append(f"# Occluded residues: {n_occluded}")
    lines.append(f"# Fraction occluded: {fraction_occluded:.2%}")
    lines.append(f"# SASA threshold: {threshold} Å²")
    lines.append(f"# Average percent burial: {avg_percent_burial:.1f}%")
    lines.append(f"# Heatmap: Blue (low ΔSASA) -> White -> Red (high ΔSASA)")
    lines.append(f"# Backbone colored by delta SASA values (protein-peptide occlusion)")
    lines.append(f"# Modified PDB file: {modified_pdb_path}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"[DEBUG] PyMOL script saved to: {out_path}")
    print(f"[DEBUG] Modified PDB file: {modified_pdb_path}")
    print(f"[DEBUG] Total PyMOL commands generated: {len(lines)}")
    if delta_sasa_values:
        min_delta = min(delta_sasa_values)
        max_delta = max(delta_sasa_values)
        print(f"[DEBUG] Key commands:")
        print(f"[DEBUG]   - Load: load {modified_pdb_path}, structure")
        print(f"[DEBUG]   - Spectrum: spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum={min_delta:.2f}, maximum={max_delta:.2f}")
        print(f"[DEBUG]   - B-factors pre-loaded from modified PDB file")


def create_pdb_with_csp_bfactors(
    results: List[CSPResult], 
    pdb_id: str, 
    out_path: str,
    sasa_results: Optional[dict] = None,
    binding_results: Optional[dict] = None,
) -> str:
    """
    Create a new PDB file with B-factor column updated with CSP values.
    Uses medoid model selection for multi-model PDB files.
    Uses sequence alignment to map CSP sequential positions to PDB residue numbers.
    Sets B-factors to zero for the peptide chain (smaller chain).
    
    Args:
        results: List of CSPResult objects from CSP analysis
        pdb_id: PDB identifier for the structure
        out_path: Path to save the modified PDB file
        sasa_results: Optional results dictionary from compute_sasa_occlusion (for chain info and alignment)
        binding_results: Optional results dictionary from interaction analysis (for chain info and alignment)
        
    Returns:
        Path to the created PDB file
    """
    import os
    import tempfile
    import shutil
    
    # Import the medoid selection function from sasa_analysis
    try:
        from .sasa_analysis import find_medoid_model_from_pdb
    except Exception:
        import os as _os, sys as _sys
        _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
        from scripts.sasa_analysis import find_medoid_model_from_pdb
    
    # Input PDB file path
    input_pdb_path = f"./PDB_FILES/{pdb_id}.pdb"
    
    if not os.path.exists(input_pdb_path):
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb_path}")
    
    # Find medoid model if it's a multi-model PDB
    medoid_pdb_path = find_medoid_model_from_pdb(input_pdb_path)
    print(f"[DEBUG] Using PDB file: {medoid_pdb_path}")
    
    # Get chain information from sasa_results or binding_results
    alignment_results = binding_results or sasa_results or {}
    receptor_chain = alignment_results.get('receptor_chain', 'A')
    ligand_chain = alignment_results.get('ligand_chain', 'B')
    
    print(f"[DEBUG] Creating PDB with CSP B-factors:")
    print(f"[DEBUG] Input PDB: {input_pdb_path}")
    print(f"[DEBUG] Medoid PDB: {medoid_pdb_path}")
    print(f"[DEBUG] Output PDB: {out_path}")
    print(f"[DEBUG] Receptor chain: {receptor_chain}")
    print(f"[DEBUG] Ligand chain: {ligand_chain}")
    
    # Create dictionary for quick lookup of CSP values by holo residue number
    csp_lookup = {}
    for result in results:
        if result.csp_A is not None:
            csp_lookup[result.holo_index] = result.csp_A
    
    print(f"[DEBUG] Processing {len(csp_lookup)} residues with CSP data")
    
    # Use sequence alignment to create mapping from CSP sequential positions to PDB residue numbers
    # This is critical because sequential positions (1, 2, 3...) may not match PDB residue numbers
    sequential_to_pdb_map = {}
    
    if alignment_results:
        # Use the existing alignment function that properly aligns sequences
        try:
            sequential_to_pdb_map = create_sequence_alignment_map_from_results(results, alignment_results)
            print(f"[DEBUG] Created alignment mapping with {len(sequential_to_pdb_map)} entries")
            if sequential_to_pdb_map:
                print(f"[DEBUG] First 10 mappings: {list(sequential_to_pdb_map.items())[:10]}")
        except Exception as e:
            print(f"[DEBUG] Warning: Sequence alignment failed: {e}")
            print(f"[DEBUG] Falling back to simple sequential mapping (may be incorrect)")
            # Fallback: create simple sequential mapping (not ideal but better than nothing)
            pdb_residue_numbers = []
            with open(medoid_pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        chain_id = line[21:22].strip()
                        if chain_id == receptor_chain:
                            res_num = int(line[22:26].strip())
                            if res_num not in pdb_residue_numbers:
                                pdb_residue_numbers.append(res_num)
            pdb_residue_numbers.sort()
            for i, pdb_res_num in enumerate(pdb_residue_numbers, 1):
                sequential_to_pdb_map[i] = pdb_res_num
    else:
        print(f"[DEBUG] Warning: No alignment results provided, using simple sequential mapping")
        # Fallback: create simple sequential mapping
        pdb_residue_numbers = []
        with open(medoid_pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21:22].strip()
                    if chain_id == receptor_chain:
                        res_num = int(line[22:26].strip())
                        if res_num not in pdb_residue_numbers:
                            pdb_residue_numbers.append(res_num)
        pdb_residue_numbers.sort()
        for i, pdb_res_num in enumerate(pdb_residue_numbers, 1):
            sequential_to_pdb_map[i] = pdb_res_num
    
    # Create reverse mapping: pdb_residue_number -> csp_value
    pdb_to_csp = {}
    for sequential_pos, pdb_res_num in sequential_to_pdb_map.items():
        if sequential_pos in csp_lookup:
            pdb_to_csp[pdb_res_num] = csp_lookup[sequential_pos]
    
    print(f"[DEBUG] Created PDB-to-CSP mapping with {len(pdb_to_csp)} entries")
    
    # Read medoid PDB and modify B-factors
    modified_lines = []
    atoms_modified = 0
    peptide_atoms_zeroed = 0
    protein_atoms_zeroed = 0
    
    with open(medoid_pdb_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract chain ID and residue number
                chain_id = line[21:22].strip()
                pdb_res_num = int(line[22:26].strip())
                
                # Set B-factor to zero for peptide chain (smaller chain)
                if chain_id == ligand_chain:
                    bfactor_str = "  0.00"
                    modified_line = line[:60] + bfactor_str + line[66:]
                    modified_lines.append(modified_line)
                    peptide_atoms_zeroed += 1
                # Apply CSP values to protein chain residues WITH CSP data
                elif chain_id == receptor_chain and pdb_res_num in pdb_to_csp:
                    # Update B-factor column (columns 60-66) with CSP value
                    csp_value = pdb_to_csp[pdb_res_num]
                    # Format B-factor to 6.2f (6 characters, 2 decimal places)
                    bfactor_str = f"{csp_value:6.2f}"
                    
                    # Replace B-factor in the line
                    modified_line = line[:60] + bfactor_str + line[66:]
                    modified_lines.append(modified_line)
                    atoms_modified += 1
                    
                    if atoms_modified <= 5:  # Debug first 5 atoms
                        atom_name = line[12:16].strip()
                        print(f"[DEBUG] Modified atom {atom_name} in PDB residue {pdb_res_num} (chain {chain_id}): B-factor = {csp_value:.2f}")
                # Set B-factor to zero for protein chain residues WITHOUT CSP data
                elif chain_id == receptor_chain:
                    bfactor_str = "  0.00"
                    modified_line = line[:60] + bfactor_str + line[66:]
                    modified_lines.append(modified_line)
                    protein_atoms_zeroed += 1
                else:
                    # Set B-factor to zero for other chains as well
                    bfactor_str = "  0.00"
                    modified_line = line[:60] + bfactor_str + line[66:]
                    modified_lines.append(modified_line)
            else:
                # Keep non-ATOM lines as-is
                modified_lines.append(line)
    
    # Write modified PDB file
    with open(out_path, 'w') as f:
        f.writelines(modified_lines)
    
    print(f"[DEBUG] PDB file created: {out_path}")
    print(f"[DEBUG] Total protein chain atoms modified with CSP: {atoms_modified}")
    print(f"[DEBUG] Total protein chain atoms set to zero (no CSP): {protein_atoms_zeroed}")
    print(f"[DEBUG] Total peptide chain atoms set to zero: {peptide_atoms_zeroed}")
    if csp_lookup:
        min_csp = min(csp_lookup.values())
        max_csp = max(csp_lookup.values())
        print(f"[DEBUG] CSP range: {min_csp:.2f} to {max_csp:.2f}")
    
    # Clean up temporary files if we created them
    if medoid_pdb_path != input_pdb_path:
        try:
            temp_dir = os.path.dirname(medoid_pdb_path)
            shutil.rmtree(temp_dir)
            print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
        except Exception as e:
            print(f"[DEBUG] Warning: Could not clean up temporary files: {e}")
    
    return out_path


def write_pymol_csp_heatmap_script(
    results: List[CSPResult],
    sasa_results: dict,
    pdb_id: str,
    out_path: str,
    binary_mode: bool = False,
    significance_field: str = 'significant',
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate PyMOL script to visualize CSP values as backbone coloring.
    Creates a modified PDB file with B-factors updated and loads that instead.
    
    Args:
        results: List of CSPResult objects from CSP analysis
        sasa_results: Results dictionary from compute_sasa_occlusion (for chain info)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        binary_mode: If True, color residues red/blue based on CSP significance.
        significance_field: CSPResult attribute used to determine binary significance.
    """
    import os
    
    # Get chain information from sasa_results
    receptor_chain = sasa_results.get('receptor_chain', 'A')
    ligand_chain = sasa_results.get('ligand_chain', 'B')
    
    print(f"[DEBUG] PyMOL CSP Heatmap Generation:")
    print(f"[DEBUG] Receptor chain: {receptor_chain}")
    print(f"[DEBUG] Ligand chain: {ligand_chain}")
    print(f"[DEBUG] Number of CSP results: {len(results)}")
    
    # Create output directory if it doesn't exist
    script_output_dir = os.path.dirname(out_path)
    if script_output_dir:
        os.makedirs(script_output_dir, exist_ok=True)
    
    # Use provided output_dir or derive from script location
    # For file operations, convert relative path to absolute if needed
    if output_dir and output_dir.startswith('./'):
        # Convert relative path to absolute using a stable project root.
        # out_path can be a temporary file path during session generation.
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        pdb_output_dir_abs = os.path.join(project_root, output_dir[2:].rstrip('/'))
    else:
        pdb_output_dir_abs = output_dir or script_output_dir
    
    # Create modified PDB file with CSP B-factors
    modified_pdb_path = os.path.join(pdb_output_dir_abs, f"{pdb_id}_csp.pdb")
    create_pdb_with_csp_bfactors(results, pdb_id, modified_pdb_path, sasa_results=sasa_results)
    
    # Generate PyMOL script - use path relative to project root
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_relative_path = normalized_dir + f"{pdb_id}_csp.pdb"
    else:
        pdb_relative_path = modified_pdb_path
    
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {pdb_relative_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    # Color peptide chain in green
    lines.append(f"color green, structure and chain {ligand_chain}")
    
    # Calculate CSP range for spectrum
    csp_values = [r.csp_A for r in results if r.csp_A is not None]
    significant_residues: Set[int] = set()
    
    if binary_mode:
        significant_residues = _get_significant_residue_numbers(results, sasa_results, significance_field)
        print(f"[DEBUG] Applying binary CSP coloring using field '{significance_field}' ({len(significant_residues)} residues)")
        lines.append(f"# Binary CSP coloring ({significance_field}): red=significant, blue=non-significant")
        lines.append(f"color blue, structure and chain {receptor_chain}")
        for res_num in sorted(significant_residues):
            lines.append(f"color red, structure and chain {receptor_chain} and resi {res_num}")
    else:
        if csp_values:
            min_csp = min(csp_values)
            max_csp = max(csp_values)
            print(f"[DEBUG] CSP range: {min_csp:.2f} to {max_csp:.2f}")
            print(f"[DEBUG] Applying spectrum coloring with min={min_csp:.2f}, max={max_csp:.2f}")
            
            lines.append(f"# CSP range: {min_csp:.2f} to {max_csp:.2f}")
            
            # Apply heatmap coloring: Blue (low) -> White -> Red (high CSP)
            lines.append(f"spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum={min_csp:.2f}, maximum={max_csp:.2f}")
        else:
            print("[DEBUG] No CSP values found, using auto-scaling")
            lines.append(f"spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum=0, maximum=auto")
    
    # Set transparency and visual enhancements
    lines.append("set cartoon_transparency, 0.2, structure")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_ring_mode, 1")
    
    # Add colorbar/legend information
    if binary_mode:
        lines.append("# Coloring legend:")
        lines.append("# Blue: Non-significant CSP residue")
        lines.append("# Red: Significant CSP residue")
    else:
        lines.append("# Colorbar information:")
        lines.append("# Blue regions: Low CSP (minimal perturbation)")
        lines.append("# White regions: Medium CSP (moderate perturbation)")
        lines.append("# Red regions: High CSP (strong perturbation)")
    
    # Add title and summary
    significant_count = sum(1 for r in results if getattr(r, significance_field, False))
    total_csp_count = len([r for r in results if r.csp_A is not None])
    
    lines.append(f"# CSP Heatmap Analysis Summary:")
    lines.append(f"# Receptor chain: {receptor_chain}")
    lines.append(f"# Ligand chain: {ligand_chain}")
    lines.append(f"# Total CSPs: {total_csp_count}")
    lines.append(f"# Significant CSPs: {significant_count}")
    if binary_mode:
        lines.append(f"# Binary coloring ({significance_field}): red=significant, blue=non-significant")
        lines.append(f"# Backbone colored by CSP significance (chemical shift perturbations)")
    else:
        lines.append(f"# Heatmap: Blue (low CSP) -> White -> Red (high CSP)")
        lines.append(f"# Backbone colored by CSP values (chemical shift perturbations)")
    lines.append(f"# Modified PDB file: {modified_pdb_path}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"[DEBUG] PyMOL script saved to: {out_path}")
    print(f"[DEBUG] Modified PDB file: {modified_pdb_path}")
    print(f"[DEBUG] Total PyMOL commands generated: {len(lines)}")
    if csp_values:
        min_csp = min(csp_values)
        max_csp = max(csp_values)
        print(f"[DEBUG] Key commands:")
        print(f"[DEBUG]   - Load: load {modified_pdb_path}, structure")
        print(f"[DEBUG]   - Spectrum: spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum={min_csp:.2f}, maximum={max_csp:.2f}")
        print(f"[DEBUG]   - B-factors pre-loaded from modified PDB file")


def write_pymol_csp_session_file(
    results: List[CSPResult],
    sasa_results: dict,
    pdb_id: str,
    out_path: str,
    binary_mode: bool = False,
    significance_field: str = 'significant',
    output_dir: Optional[str] = None,
) -> bool:
    """
    Generate a PyMOL session file (.pse) with CSP coloring.
    
    Args:
        results: List of CSPResult objects from CSP analysis
        sasa_results: Results dictionary from compute_sasa_occlusion
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL session file
        binary_mode: If True, color CSP residues using binary significance coloring (red/blue).
        significance_field: CSPResult attribute to use for significance classification.
        
    Returns:
        True if successful, False otherwise
    """
    import tempfile
    import subprocess
    
    # Create temporary PyMOL script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as temp_script:
        write_pymol_csp_heatmap_script(
            results,
            sasa_results,
            pdb_id,
            temp_script.name,
            binary_mode=binary_mode,
            significance_field=significance_field,
            output_dir=output_dir,
        )
        temp_script_path = temp_script.name
    
    try:
        # Generate PyMOL session file
        pymol_cmd = [
            'pymol', '-c', '-q',  # Command line mode, quiet
            temp_script_path,      # Run our script
            '-d', f'save {out_path}',  # Save session file
            '-d', 'quit'          # Quit PyMOL
        ]
        result = subprocess.run(pymol_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: PyMOL CSP session generation failed: {result.stderr}")
            print("PyMOL CSP script saved instead. You can run it manually in PyMOL.")
            return False
        else:
            print(f"✓ PyMOL CSP session file saved: {out_path}")
            return True
            
    except FileNotFoundError:
        print("Warning: PyMOL not found. PyMOL CSP script saved instead.")
        return False
    except Exception as e:
        print(f"Warning: Error generating PyMOL CSP session: {e}")
        return False
    finally:
        # Clean up temporary file
        import os
        try:
            os.unlink(temp_script_path)
        except:
            pass


def write_pymol_session_file(
    sasa_results: dict,
    pdb_id: str,
    out_path: str,
    binary_mode: bool = False,
    output_dir: Optional[str] = None,
) -> bool:
    """
    Generate a PyMOL session file (.pse) with delta SASA coloring.
    
    Args:
        sasa_results: Results dictionary from compute_sasa_occlusion
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL session file
        binary_mode: If True, color using binary occlusion categories (red/blue) instead of a spectrum.
    
    Returns:
        True if successful, False otherwise
    """
    import tempfile
    import subprocess
    
    # Create temporary PyMOL script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as temp_script:
        write_pymol_delta_sasa_script(sasa_results, pdb_id, temp_script.name, binary_mode=binary_mode, output_dir=output_dir)
        temp_script_path = temp_script.name
    
    try:
        # Generate PyMOL session file
        pymol_cmd = [
            'pymol', '-c', '-q',  # Command line mode, quiet
            temp_script_path,      # Run our script
            '-d', f'save {out_path}',  # Save session file
            '-d', 'quit'          # Quit PyMOL
        ]
        
        result = subprocess.run(pymol_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: PyMOL session generation failed: {result.stderr}")
            print("PyMOL script saved instead. You can run it manually in PyMOL.")
            return False
        else:
            print(f"✓ PyMOL session file saved: {out_path}")
            return True
            
    except FileNotFoundError:
        print("Warning: PyMOL not found. PyMOL script saved instead.")
        return False
    except Exception as e:
        print(f"Warning: Error generating PyMOL session: {e}")
        return False
    finally:
        # Clean up temporary file
        import os
        try:
            os.unlink(temp_script_path)
        except:
            pass

def write_pymol_occlusion_script(
    sasa_results: dict, 
    pdb_id: str, 
    out_path: str,
    interaction_results: Optional[dict] = None,
    ca_distance_results: Optional[dict] = None,
    any_atom_results: Optional[dict] = None,
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate PyMOL script to visualize binding site residues.
    
    Binding site residues are identified as: is_occluded OR passes_ca_filter OR is_interacting
    OR min_any_atom_distance < 2A (matching the 'occluded_or_ca_or_interaction' positive strategy)
    
    Args:
        sasa_results: Results dictionary from compute_sasa_occlusion
        output_dir: Output directory path (used to construct PDB file paths)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        interaction_results: Optional results dictionary from compute_interaction_filter
        ca_distance_results: Optional results dictionary from compute_ca_distance_filter
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from sasa_results)
    """
    lines = []
    lines.append("reinitialize")
    # Construct PDB path relative to project root
    import os
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        # Ensure it ends with / and use string concatenation (os.path.join doesn't handle ./ well)
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_delta_sasa.pdb"
        # Convert relative path to absolute for file operations
        if pdb_path.startswith('./'):
            # Get project root: go up from out_path until we find the directory containing 'outputs'
            current = os.path.dirname(os.path.abspath(out_path))
            while current != os.path.dirname(current):  # Stop at filesystem root
                if os.path.basename(current) == 'outputs':
                    project_root = os.path.dirname(current)
                    break
                current = os.path.dirname(current)
            else:
                # Fallback: assume project root is 2 levels up from out_path
                project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
            pdb_path_abs = os.path.join(project_root, pdb_path[2:])
        else:
            pdb_path_abs = pdb_path
        if not pdb_path.startswith('.'):
            pdb_path = './' + pdb_path
    else:
        pdb_path = f"./outputs/{pdb_id}/{pdb_id}_delta_sasa.pdb"
        # Get project root: go up from out_path until we find the directory containing 'outputs'
        current = os.path.dirname(os.path.abspath(out_path))
        while current != os.path.dirname(current):  # Stop at filesystem root
            if os.path.basename(current) == 'outputs':
                project_root = os.path.dirname(current)
                break
            current = os.path.dirname(current)
        else:
            # Fallback: assume project root is 2 levels up from out_path
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_path_abs = os.path.join(project_root, pdb_path[2:])
    
    # Create the delta_sasa.pdb file if it doesn't exist
    if sasa_results and not os.path.exists(pdb_path_abs):
        create_pdb_with_delta_sasa_bfactors(sasa_results, pdb_id, pdb_path_abs)
    
    lines.append(f"load {pdb_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    # Get chain information from parameters, fallback to sasa_results if not provided
    if receptor_chain is None:
        receptor_chain = sasa_results.get('receptor_chain', 'A')
    if ligand_chain is None:
        ligand_chain = sasa_results.get('ligand_chain', 'B')
    
    # Identify binding site residues: is_occluded OR passes_ca_filter OR is_interacting
    binding_site_residues: Set[int] = set()
    
    # Add occluded residues
    for info in sasa_results.get('residue_info', []):
        if info.get('is_occluded') and info.get('residue_number') is not None:
            binding_site_residues.add(info['residue_number'])
    
    # Add residues that pass CA distance filter
    if ca_distance_results:
        for info in ca_distance_results.get('residue_info', []):
            if info.get('passes_filter') and info.get('residue_number') is not None:
                binding_site_residues.add(info['residue_number'])
    
    # Add interacting residues (has_hbond OR has_charge_complement OR has_pi_contact)
    if interaction_results:
        for info in interaction_results.get('residue_info', []):
            res_num = info.get('residue_number')
            if res_num is not None:
                is_interacting = (
                    info.get('has_hbond') or 
                    info.get('has_charge_complement') or 
                    info.get('has_pi_contact')
                )
                if is_interacting:
                    binding_site_residues.add(res_num)
    
    # Add residues with min inter-chain atom-atom distance < 2A
    if any_atom_results:
        for info in any_atom_results.get('residue_info', []):
            if info.get('passes_sub_2A_filter') and info.get('residue_number') is not None:
                binding_site_residues.add(info['residue_number'])
    
    # Color all receptor chain residues gray30 (default)
    lines.append(f"color gray30, structure and chain {receptor_chain}")
    
    # Color binding site residues red
    if binding_site_residues:
        for res_num in sorted(binding_site_residues):
            lines.append(f"color red, structure and chain {receptor_chain} and resi {res_num}")
    
    # Color peptide chain cyan
    lines.append(f"color cyan, structure and chain {ligand_chain}")
    
    # Set transparency
    lines.append("set cartoon_transparency, 0.2, structure")
    
    # Add summary comments
    lines.append(f"# Binding Site Analysis Summary:")
    lines.append(f"# Receptor chain: {receptor_chain}")
    lines.append(f"# Ligand chain: {ligand_chain}")
    lines.append(f"# Non-binding receptor residues: gray30")
    lines.append(f"# Binding site residues (red): {len(binding_site_residues)}")
    lines.append(f"# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A")
    if binding_site_residues:
        lines.append(f"# Binding site residue numbers: {sorted(binding_site_residues)}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_pymol_combined_script(
    results: List[CSPResult],
    binding_results: dict,
    pdb_id: str,
    out_path: str,
    binary_mode: bool = False,
    significance_field: str = 'significant',
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    sasa_results: Optional[dict] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate PyMOL script combining CSP and occlusion visualizations.
    
    Args:
        results: CSP analysis results
        binding_results: Residue annotation results (SASA occlusion or interactions)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        binary_mode: If True, color residues red/blue using CSP significance instead of a heatmap.
        significance_field: CSPResult attribute used to determine binary significance.
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from binding_results)
    """
    lines = []
    lines.append("reinitialize")
    # Construct PDB path relative to project root
    import os
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        # Ensure it ends with / and use string concatenation (os.path.join doesn't handle ./ well)
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_delta_sasa.pdb"
        # Convert relative path to absolute for file operations
        if pdb_path.startswith('./'):
            # Get project root: go up from out_path until we find the directory containing 'outputs'
            current = os.path.dirname(os.path.abspath(out_path))
            while current != os.path.dirname(current):  # Stop at filesystem root
                if os.path.basename(current) == 'outputs':
                    project_root = os.path.dirname(current)
                    break
                current = os.path.dirname(current)
            else:
                # Fallback: assume project root is 2 levels up from out_path
                project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
            pdb_path_abs = os.path.join(project_root, pdb_path[2:])
        else:
            pdb_path_abs = pdb_path
        if not pdb_path.startswith('.'):
            pdb_path = './' + pdb_path
    else:
        pdb_path = f"./outputs/{pdb_id}/{pdb_id}_delta_sasa.pdb"
        # Get project root: go up from out_path until we find the directory containing 'outputs'
        current = os.path.dirname(os.path.abspath(out_path))
        while current != os.path.dirname(current):  # Stop at filesystem root
            if os.path.basename(current) == 'outputs':
                project_root = os.path.dirname(current)
                break
            current = os.path.dirname(current)
        else:
            # Fallback: assume project root is 2 levels up from out_path
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_path_abs = os.path.join(project_root, pdb_path[2:])
    
    # Create the delta_sasa.pdb file if it doesn't exist
    # Use sasa_results if provided, otherwise try binding_results (which might be interaction_results)
    if not os.path.exists(pdb_path_abs):
        if sasa_results:
            create_pdb_with_delta_sasa_bfactors(sasa_results, pdb_id, pdb_path_abs)
        elif binding_results and 'residue_info' in binding_results:
            # Check if binding_results has delta_sasa data (it might be sasa_results passed as binding_results)
            try:
                # Try to create it - will fail if binding_results doesn't have the right structure
                create_pdb_with_delta_sasa_bfactors(binding_results, pdb_id, pdb_path_abs)
            except (KeyError, AttributeError, TypeError):
                # If binding_results doesn't have the right structure, skip creation
                # The file should be created by write_pymol_delta_sasa_script later
                pass
    
    lines.append(f"load {pdb_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    # Get chain information from parameters, fallback to binding_results if not provided
    if receptor_chain is None:
        receptor_chain = binding_results.get('receptor_chain', 'A')
    if ligand_chain is None:
        ligand_chain = binding_results.get('ligand_chain', 'B')

    binding_lookup, dataset_type = _build_binding_lookup(binding_results)
    binding_residues = sorted(res for res, flag in binding_lookup.items() if flag)
    
    # Color peptide chain in green
    lines.append(f"color green, structure and chain {ligand_chain}")
    
    # Color protein chain by CSP values or binary significance
    if binary_mode:
        significant_residues = _get_significant_residue_numbers(results, binding_results, significance_field)
        print(f"[DEBUG] Combined view binary coloring ({significance_field}) with {len(significant_residues)} significant residues")
        lines.append(f"# Binary CSP coloring ({significance_field}): red=significant, blue=non-significant")
        lines.append(f"color blue, structure and chain {receptor_chain}")
        for res_num in sorted(significant_residues):
            lines.append(f"color red, structure and chain {receptor_chain} and resi {res_num}")
    else:
            lines.append(f"spectrum b, blue_white_red, structure and chain {receptor_chain}, minimum=0, maximum=auto")
    
    # Highlight occluded residues with spheres
    if binding_residues:
        for binding_residue in binding_residues:
            lines.append(f"color yellow, structure and chain {receptor_chain} and resi {binding_residue}")
            lines.append(f"show spheres, structure and chain {receptor_chain} and resi {binding_residue}")
            lines.append(f"set sphere_scale, 0.4, structure and chain {receptor_chain} and resi {binding_residue}")
    
    # Set transparency
    lines.append("set cartoon_transparency, 0.2, structure")
    
    # Add summary information
    n_binding = len(binding_residues)
    lines.append("# Combined CSP and Binding Analysis:")
    if binary_mode:
        lines.append(f"# Receptor chain: {receptor_chain} (red=significant via {significance_field}, blue=non-significant)")
    else:
        lines.append(f"# Receptor chain: {receptor_chain} (colored by CSP)")
    lines.append(f"# Ligand chain: {ligand_chain} (green)")
    if dataset_type == 'interaction':
        interaction_counts = binding_results.get('n_interacting_residues', n_binding)
        lines.append(f"# Interacting residues (yellow spheres): {interaction_counts}")
        lines.append(f"#   Hydrogen bond residues: {binding_results.get('n_hbond_residues', 0)}")
        lines.append(f"#   Charge complement residues: {binding_results.get('n_charge_residues', 0)}")
        lines.append(f"#   Pi-contact residues: {binding_results.get('n_pi_residues', 0)}")
    elif dataset_type == 'union':
        union_count = binding_results.get('n_union_residues', n_binding)
        lines.append(f"# Union-positive residues (yellow spheres): {union_count}")
        lines.append(f"#   Hydrogen bond residues: {binding_results.get('n_hbond_residues', 0)}")
        lines.append(f"#   Charge complement residues: {binding_results.get('n_charge_residues', 0)}")
        lines.append(f"#   Pi-contact residues: {binding_results.get('n_pi_residues', 0)}")
        lines.append(f"#   SASA occluded residues: {binding_results.get('n_sasa_residues', 0)}")
        lines.append(f"#   CA-distance residues: {binding_results.get('n_ca_residues', 0)}")
    elif dataset_type == 'sasa':
        fraction_occluded = binding_results.get('fraction_occluded', 0.0)
        threshold = binding_results.get('sasa_threshold', 5.0)
        lines.append(f"# Occluded residues (yellow spheres): {n_binding}")
        lines.append(f"# Fraction occluded: {fraction_occluded:.2%}")
        lines.append(f"# SASA threshold: {threshold} Å²")
    else:
        lines.append(f"# Binding residues (yellow spheres): {n_binding}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _get_significant_residue_numbers(
    results: List[CSPResult],
    binding_results: dict,
    significance_field: str = 'significant',
) -> Set[int]:
    """
    Map CSP significance flags to PDB residue numbers for coloring.
    
    Args:
        results: CSP analysis results.
        binding_results: Residue annotation results (provides sequence alignment context).
        significance_field: Attribute name on CSPResult indicating binary significance.
    
    Returns:
        Set of PDB residue numbers corresponding to significant CSP residues.
    """
    significant_residues: Set[int] = set()
    try:
        position_map = create_sequence_alignment_map_from_results(results, binding_results)
        for r in results:
            if getattr(r, significance_field, False):
                pdb_residue_number = position_map.get(r.holo_index)
                if pdb_residue_number is not None:
                    significant_residues.add(pdb_residue_number)
    except Exception as exc:
        print(f"[VISUALIZE WARNING] Failed to map CSP significance with field '{significance_field}': {exc}")
    return significant_residues


def write_pymol_interaction_union_script(
    results: List[CSPResult],
    sasa_results: dict,
    interaction_results: dict,
    ca_distance_results: dict,
    pdb_id: str,
    out_path: str,
    significance_field: str = 'significant',
    receptor_chain: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate a PyMOL script highlighting residues predicted positive by any strategy.
    
    A residue is colored red if any of the following are true:
    - CSP result flagged as significant (per `significance_field`)
    - SASA occlusion marks residue as occluded
    - Interaction analysis detects hydrogen bond, charge complementarity, or pi contact
    - CA distance filter marks residue within proximity threshold
    
    Args:
        results: CSP analysis results.
        sasa_results: SASA occlusion results.
        interaction_results: Interaction filter results.
        ca_distance_results: CA distance filter results.
        pdb_id: PDB identifier.
        out_path: Destination for the PyMOL script.
        significance_field: CSPResult attribute to use for binary classification.
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from interaction_results or sasa_results)
    """
    import os
    
    # Get chains from parameters, fallback to interaction_results or sasa_results if not provided
    if receptor_chain is None:
        receptor_chain = (
            (interaction_results or {}).get('receptor_chain')
            or sasa_results.get('receptor_chain', 'A')
        )
    if ligand_chain is None:
        ligand_chain = (
            (interaction_results or {}).get('ligand_chain')
            or sasa_results.get('ligand_chain', 'B')
        )
    
    # Create output directory if it doesn't exist (don't overwrite output_dir parameter)
    script_output_dir = os.path.dirname(out_path)
    if script_output_dir:
        os.makedirs(script_output_dir, exist_ok=True)
    
    csp_positive = _get_significant_residue_numbers(results, sasa_results, significance_field)
    sasa_positive = {
        info['residue_number']
        for info in sasa_results.get('residue_info', [])
        if info.get('is_occluded') and info.get('residue_number') is not None
    }
    interaction_info = (interaction_results or {}).get('residue_info', [])
    hbond_positive = {
        info['residue_number']
        for info in interaction_info
        if info.get('has_hbond') and info.get('residue_number') is not None
    }
    charge_positive = {
        info['residue_number']
        for info in interaction_info
        if info.get('has_charge_complement') and info.get('residue_number') is not None
    }
    pi_positive = {
        info['residue_number']
        for info in interaction_info
        if info.get('has_pi_contact') and info.get('residue_number') is not None
    }
    interaction_positive = hbond_positive | charge_positive | pi_positive
    ca_positive = {
        info['residue_number']
        for info in (ca_distance_results or {}).get('residue_info', [])
        if info.get('passes_filter') and info.get('residue_number') is not None
    }
    
    union_positive = sasa_positive | interaction_positive | ca_positive
    
    print(f"[DEBUG] Interaction union visualization: {len(union_positive)} residues marked positive")
    
    lines: List[str] = []
    lines.append("reinitialize")
    # Construct PDB path relative to project root
    import os
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        # Ensure it ends with / and use string concatenation (os.path.join doesn't handle ./ well)
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_path = normalized_dir + f"{pdb_id}_delta_sasa.pdb"
        if not pdb_path.startswith('.'):
            pdb_path = './' + pdb_path
    else:
        pdb_path = f"./outputs/{pdb_id}/{pdb_id}_delta_sasa.pdb"
    lines.append(f"load {pdb_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    lines.append(f"color blue, structure and chain {receptor_chain}")
    lines.append(f"color green, structure and chain {ligand_chain}")
    
    for res_num in sorted(union_positive):
        lines.append(f"color red, structure and chain {receptor_chain} and resi {res_num}")
    
    lines.append("set cartoon_transparency, 0.2, structure")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_ring_mode, 1")
    
    lines.append("# Interaction union summary:")
    lines.append(f"# Positive residues (union): {len(union_positive)}")
    lines.append(f"# - CSP significant ({significance_field}): {len(csp_positive)}")
    lines.append(f"# - SASA occluded: {len(sasa_positive)}")
    lines.append(f"# - Interaction (H-bond/charge/pi): {len(interaction_positive)}")
    lines.append(f"#   * H-bond: {len(hbond_positive)}")
    lines.append(f"#   * Charge complementarity: {len(charge_positive)}")
    lines.append(f"#   * Pi contact: {len(pi_positive)}")
    lines.append(f"# - CA proximity filter: {len(ca_positive)}")
    lines.append(f"# Residues colored red if positive by ANY strategy listed above.")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"[DEBUG] PyMOL interaction union script saved to: {out_path}")
    print(f"[DEBUG] Total PyMOL commands generated: {len(lines)}")


def plot_csp_distribution(results: List[CSPResult], out_png: str, title: Optional[str] = None) -> None:
    if not _HAS_PLT:
        return
    
    # Filter results with valid CSP values
    valid_results = [r for r in results if r.csp_A is not None]
    if not valid_results:
        return
    
    # Ensure we're using non-interactive backend
    import matplotlib
    matplotlib.use('Agg', force=True)
    
    # Extract data for plotting
    positions = []
    csp_values = []
    significant_values = []
    significant_positions = []
    
    for r in valid_results:
        # Use apo_index as the sequence position (1-based)
        positions.append(r.apo_index)
        csp_values.append(r.csp_A)
        
        # Separate significant CSPs for highlighting
        if r.significant:
            significant_positions.append(r.apo_index)
            significant_values.append(r.csp_A)
    
    # Create the plot
    plt.figure(figsize=(12, 6))
    
    # Plot all CSPs as a line
    plt.plot(positions, csp_values, 'o-', color='#4C78A8', markersize=4, linewidth=1, alpha=0.7, label='All CSPs')
    
    # Highlight significant CSPs
    if significant_values:
        plt.scatter(significant_positions, significant_values, color='red', s=20, alpha=0.8, label='Significant CSPs', zorder=5)
    
    # Add horizontal line at significance threshold (if we can determine it)
    if significant_values:
        # Estimate threshold as the minimum significant CSP value
        threshold = min(significant_values)
        plt.axhline(y=threshold, color='red', linestyle='--', alpha=0.5, label=f'Threshold ≈ {threshold:.3f}')
    
    # Set x-axis ticks with amino acid sequence
    plt.xticks(positions, [r.apo_aa for r in valid_results], fontsize=8)
    plt.xlabel("Sequence Position (Apo)")
    plt.ylabel("CSP Magnitude")
    plt.title(title or "Chemical Shift Perturbations Along Sequence")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Create three-row residue index system
    residue_indices = [r.apo_index for r in valid_results]
    ones_ticks, tens_ticks, hundreds_ticks = create_three_row_residue_ticks(residue_indices)
    
    # Add three secondary x-axes for different digit places
    ax2 = plt.gca().secondary_xaxis('bottom')
    ax3 = plt.gca().secondary_xaxis('bottom')
    ax4 = plt.gca().secondary_xaxis('bottom')
    
    # Row 1: Ones place (closest to main axis)
    ax2.set_xticks([pos for pos, _ in ones_ticks])
    ax2.set_xticklabels([label for _, label in ones_ticks], fontsize=5)
    ax2.spines['bottom'].set_position(('outward', 15))
    
    # Row 2: Tens place (middle)
    if tens_ticks:
        ax3.set_xticks([pos for pos, _ in tens_ticks])
        ax3.set_xticklabels([label for _, label in tens_ticks], fontsize=5)
        ax3.spines['bottom'].set_position(('outward', 30))
    
    # Row 3: Hundreds place (furthest from main axis)
    if hundreds_ticks:
        ax4.set_xticks([pos for pos, _ in hundreds_ticks])
        ax4.set_xticklabels([label for _, label in hundreds_ticks], fontsize=5)
        ax4.spines['bottom'].set_position(('outward', 45))
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def write_pymol_csp_classification_script(results: List[CSPResult], binding_results: dict, pdb_id: str, out_path: str, significance_field: str = 'significant', receptor_chain: Optional[str] = None, ligand_chain: Optional[str] = None, output_dir: Optional[str] = None) -> None:
    """
    Generate PyMOL script to visualize CSP classification (TP/FP/TN/FN) on protein structure.
    Creates a modified PDB file with B-factors updated and loads that instead.
    
    Args:
        results: List of CSPResult objects from CSP analysis
        binding_results: Results dictionary providing residue annotations (SASA occlusion or interactions)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL script
        significance_field: Which significance field to use ('significant', 'significant_1sd', 'significant_2sd')
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from binding_results)
    """
    import os
    
    # Get chain information from parameters, fallback to binding_results if not provided
    if receptor_chain is None:
        receptor_chain = binding_results.get('receptor_chain', 'A')
    if ligand_chain is None:
        ligand_chain = binding_results.get('ligand_chain', 'B')
    
    print(f"[DEBUG] PyMOL CSP Classification Visualization Generation:")
    print(f"[DEBUG] Receptor chain: {receptor_chain}")
    print(f"[DEBUG] Ligand chain: {ligand_chain}")
    print(f"[DEBUG] Number of CSP results: {len(results)}")
    
    # Create output directory if it doesn't exist
    script_output_dir = os.path.dirname(out_path)
    if script_output_dir:
        os.makedirs(script_output_dir, exist_ok=True)
    
    # Use provided output_dir or derive from script location
    # For file operations, convert relative path to absolute if needed
    if output_dir and output_dir.startswith('./'):
        # Get project root: go up from out_path until we find the directory containing 'outputs'
        current = os.path.dirname(os.path.abspath(out_path))
        while current != os.path.dirname(current):  # Stop at filesystem root
            if os.path.basename(current) == 'outputs':
                project_root = os.path.dirname(current)
                break
            current = os.path.dirname(current)
        else:
            # Fallback: assume project root is 2 levels up from out_path
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(out_path)))
        pdb_output_dir_abs = os.path.join(project_root, output_dir[2:].rstrip('/'))
    else:
        pdb_output_dir_abs = output_dir or script_output_dir
    
    # Create modified PDB file with CSP B-factors
    modified_pdb_path = os.path.join(pdb_output_dir_abs, f"{pdb_id}_csp.pdb")
    create_pdb_with_csp_bfactors(results, pdb_id, modified_pdb_path, binding_results=binding_results)
    
    # Create sequence alignment mapping using CSP, binding, and PDB sequences
    pdb_source_path = f"./PDB_FILES/{pdb_id}.pdb"
    if not os.path.exists(pdb_source_path):
        pdb_source_path = modified_pdb_path
    position_map = _create_csp_binding_pdb_alignment(results, binding_results, pdb_source_path, receptor_chain)

    # Create binding lookup by PDB residue number
    binding_lookup, dataset_type = _build_binding_lookup(binding_results)
    
    # Get all PDB residue numbers for the protein chain from the PDB file
    all_protein_residues: Set[int] = set()
    if os.path.exists(pdb_source_path):
        pdb_seq, pdb_positions = _extract_pdb_sequence(pdb_source_path, receptor_chain)
        all_protein_residues = set(pdb_positions)
    
    # Create classification mapping using PDB residue numbers (not sequential positions)
    classifications: Dict[int, str] = {}  # pdb_residue_number -> classification
    residues_with_csp: Set[int] = set()
    
    for r in results:
        if r.csp_A is not None:
            # Get corresponding PDB residue number from alignment
            pdb_residue_number = position_map.get(r.holo_index)
            if pdb_residue_number is None:
                print(f"[PYMOL WARNING] No alignment found for sequential position {r.holo_index}")
                continue  # Skip residues without alignment
            
            residues_with_csp.add(pdb_residue_number)
            
            # Get binding status for this PDB residue
            is_binding = binding_lookup.get(pdb_residue_number, False)
            
            # Get significance value based on the specified field
            is_significant = getattr(r, significance_field, False)
            
            if is_significant and is_binding:
                classifications[pdb_residue_number] = 'TP'  # True Positive
            elif is_significant and not is_binding:
                classifications[pdb_residue_number] = 'FP'  # False Positive
            elif not is_significant and not is_binding:
                classifications[pdb_residue_number] = 'TN'  # True Negative
            else:  # not is_significant and is_binding
                classifications[pdb_residue_number] = 'FN'  # False Negative
    
    # Identify residues without CSP data
    residues_without_csp = all_protein_residues - residues_with_csp
    
    # Count classifications
    tp_count = list(classifications.values()).count('TP')
    fp_count = list(classifications.values()).count('FP')
    tn_count = list(classifications.values()).count('TN')
    fn_count = list(classifications.values()).count('FN')
    
    print(f"[DEBUG] Classification counts - TP: {tp_count}, FP: {fp_count}, TN: {tn_count}, FN: {fn_count}")
    print(f"[DEBUG] Residues without CSP: {len(residues_without_csp)}")
    
    # Generate PyMOL script - use path relative to project root
    if output_dir:
        # output_dir should be the full path from project root (e.g., "./outputs/2mur_1/")
        normalized_dir = output_dir.replace('\\', '/')
        if not normalized_dir.endswith('/'):
            normalized_dir += '/'
        pdb_relative_path = normalized_dir + f"{pdb_id}_csp.pdb"
    else:
        pdb_relative_path = modified_pdb_path
    
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {pdb_relative_path}, structure")
    lines.append("hide everything, structure")
    lines.append("show cartoon, structure")
    
    # Color peptide chain cyan
    lines.append(f"color cyan, structure and chain {ligand_chain}")
    # Make ligand loops easier to see in cartoon view.
    lines.append(f"set cartoon_tube_radius, 0.45, structure and chain {ligand_chain}")
    
    # Color all protein chain residues gray first (for residues without CSP)
    lines.append(f"color gray, structure and chain {receptor_chain}")
    color_map = _get_classification_colors()
    tp_rgb = _config_hex_to_rgb01(color_map["TP"])
    tn_rgb = _config_hex_to_rgb01(color_map["TN"])
    fp_rgb = _config_hex_to_rgb01(color_map["FP"])
    fn_rgb = _config_hex_to_rgb01(color_map["FN"])
    lines.append(f"set_color tp_color, [{tp_rgb[0]:.4f}, {tp_rgb[1]:.4f}, {tp_rgb[2]:.4f}]")
    lines.append(f"set_color tn_color, [{tn_rgb[0]:.4f}, {tn_rgb[1]:.4f}, {tn_rgb[2]:.4f}]")
    lines.append(f"set_color fp_color, [{fp_rgb[0]:.4f}, {fp_rgb[1]:.4f}, {fp_rgb[2]:.4f}]")
    lines.append(f"set_color fn_color, [{fn_rgb[0]:.4f}, {fn_rgb[1]:.4f}, {fn_rgb[2]:.4f}]")
    
    # Color protein chain residues by classification using PDB residue numbers
    # TP: Bluish green
    tp_residues = [res_num for res_num, cls in classifications.items() if cls == 'TP']
    if tp_residues:
        lines.append(f"# TP (Sig. CSP in Binding Site): {len(tp_residues)} residues")
        for res_num in sorted(tp_residues):
            lines.append(f"color tp_color, structure and chain {receptor_chain} and resi {res_num}")
    
    # FP: Vermillion
    fp_residues = [res_num for res_num, cls in classifications.items() if cls == 'FP']
    if fp_residues:
        lines.append(f"# FP (Sig. CSP -- Allosteric): {len(fp_residues)} residues")
        for res_num in sorted(fp_residues):
            lines.append(f"color fp_color, structure and chain {receptor_chain} and resi {res_num}")
    
    # TN: Blue
    tn_residues = [res_num for res_num, cls in classifications.items() if cls == 'TN']
    if tn_residues:
        lines.append(f"# TN (low CSP -- Allosteric): {len(tn_residues)} residues")
        for res_num in sorted(tn_residues):
            lines.append(f"color tn_color, structure and chain {receptor_chain} and resi {res_num}")
    
    # FN: Reddish purple
    fn_residues = [res_num for res_num, cls in classifications.items() if cls == 'FN']
    if fn_residues:
        lines.append(f"# FN (low CSP in Binding Site): {len(fn_residues)} residues")
        for res_num in sorted(fn_residues):
            lines.append(f"color fn_color, structure and chain {receptor_chain} and resi {res_num}")
    
    # Set transparency and visual enhancements
    lines.append("set cartoon_transparency, 0.2, structure")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_ring_mode, 1")
    
    # Add colorbar/legend information
    lines.append("# CSP Classification Color Scheme:")
    classification_label = "Binding Site"
    if dataset_type == 'interaction':
        classification_label = "Interaction Site"
    elif dataset_type == 'union':
        classification_label = "Union Site"
    elif dataset_type == 'sasa':
        classification_label = "Occluded Site"

    lines.append(f"# TP ({color_map['TP']}): Sig. CSP in {classification_label}")
    lines.append(f"# FP ({color_map['FP']}): Sig. CSP -- Allosteric")
    lines.append(f"# TN ({color_map['TN']}): low CSP -- Allosteric")
    lines.append(f"# FN ({color_map['FN']}): low CSP in {classification_label}")
    lines.append("# Cyan: Peptide chain")
    lines.append("# Gray: Protein residues without CSP data")
    
    # Add title and summary
    lines.append(f"# CSP Classification Analysis Summary:")
    lines.append(f"# Receptor chain: {receptor_chain}")
    lines.append(f"# Ligand chain: {ligand_chain}")
    lines.append(f"# TP (Sig. CSP in {classification_label}): {tp_count}")
    lines.append(f"# FP (Sig. CSP -- Allosteric): {fp_count}")
    lines.append(f"# TN (low CSP -- Allosteric): {tn_count}")
    lines.append(f"# FN (low CSP in {classification_label}): {fn_count}")
    lines.append(f"# Residues without CSP data: {len(residues_without_csp)}")
    lines.append(f"# Structure colored by CSP classification (TP/FP/TN/FN)")
    lines.append(f"# Modified PDB file: {modified_pdb_path}")
    
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    
    print(f"[DEBUG] PyMOL CSP classification script saved to: {out_path}")
    print(f"[DEBUG] Modified PDB file: {modified_pdb_path}")
    print(f"[DEBUG] Total PyMOL commands generated: {len(lines)}")


def write_pymol_csp_classification_session_file(results: List[CSPResult], binding_results: dict, pdb_id: str, out_path: str, significance_field: str = 'significant', receptor_chain: Optional[str] = None, ligand_chain: Optional[str] = None, output_dir: Optional[str] = None) -> bool:
    """
    Generate a PyMOL session file (.pse) with CSP classification coloring.
    
    Args:
        results: List of CSPResult objects from CSP analysis
        binding_results: Results dictionary providing residue annotations (SASA occlusion or interactions)
        pdb_id: PDB identifier for the structure
        out_path: Path to save the PyMOL session file
        significance_field: Which significance field to use ('significant', 'significant_1sd', 'significant_2sd')
        receptor_chain: Receptor chain ID (chain with sequence present in both apo and holo BMRB)
        ligand_chain: Ligand chain ID (if not provided, will try to get from binding_results)
        
    Returns:
        True if successful, False otherwise
    """
    import tempfile
    import subprocess
    
    # Create temporary PyMOL script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as temp_script:
        write_pymol_csp_classification_script(
            results,
            binding_results,
            pdb_id,
            temp_script.name,
            significance_field,
            receptor_chain=receptor_chain,
            ligand_chain=ligand_chain,
            output_dir=output_dir,
        )
        temp_script_path = temp_script.name
    
    try:
        # Generate PyMOL session file
        pymol_cmd = [
            'pymol', '-c', '-q',  # Command line mode, quiet
            temp_script_path,      # Run our script
            '-d', f'save {out_path}',  # Save session file
            '-d', 'quit'          # Quit PyMOL
        ]
        
        result = subprocess.run(pymol_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: PyMOL CSP classification session generation failed: {result.stderr}")
            print("PyMOL CSP classification script saved instead. You can run it manually in PyMOL.")
            return False
        else:
            print(f"✓ PyMOL CSP classification session file saved: {out_path}")
            return True
            
    except FileNotFoundError:
        print("Warning: PyMOL not found. PyMOL CSP classification script saved instead.")
        return False
    except Exception as e:
        print(f"Warning: Error generating PyMOL CSP classification session: {e}")
        return False
    finally:
        # Clean up temporary file
        import os
        try:
            os.unlink(temp_script_path)
        except:
            pass


def create_sequence_alignment_map_from_results(results: List[CSPResult], binding_results: dict) -> Dict[int, int]:
    """
    Create a mapping between sequential CSP positions and PDB residue numbers
    by aligning the amino acid sequences from CSP results and SASA results.
    
    Args:
        results: CSP analysis results with sequential positions
        binding_results: Residue annotation results with PDB residue numbers
        
    Returns:
        Dictionary mapping sequential_position -> pdb_residue_number
    """
    # Extract sequences from CSP data (sequential positions)
    csp_sequence = []
    csp_positions = []
    for r in results:
        if r.holo_aa and r.holo_aa != 'P':  # Skip prolines
            csp_sequence.append(r.holo_aa)
            csp_positions.append(r.holo_index)
    
    # Extract sequences from SASA data (PDB residue numbers)
    occlusion_sequence = []
    occlusion_positions = []
    for info in binding_results.get('residue_info', []):
        residue_name = info.get('residue_name', '').strip()
        residue_number = info.get('residue_number')
        if residue_name and residue_name != 'PRO' and residue_number is not None:  # Skip prolines
            # Convert 3-letter to 1-letter amino acid code
            aa_map = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
            }
            aa_letter = aa_map.get(residue_name, 'X')
            occlusion_sequence.append(aa_letter)
            occlusion_positions.append(residue_number)
    
    # Perform sequence alignment
    try:
        from .align import align_global
    except Exception:
        import os as _os, sys as _sys
        _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
        from scripts.align import align_global
    
    csp_seq_str = ''.join(csp_sequence)
    occlusion_seq_str = ''.join(occlusion_sequence)
    
    print(f"[VISUALIZE ALIGN] CSP sequence length: {len(csp_seq_str)}")
    print(f"[VISUALIZE ALIGN] Occlusion sequence length: {len(occlusion_seq_str)}")
    print(f"[VISUALIZE ALIGN] CSP sequence: {csp_seq_str[:50]}...")
    print(f"[VISUALIZE ALIGN] Occlusion sequence: {occlusion_seq_str[:50]}...")
    
    # Perform global alignment
    aligned_csp, aligned_occlusion, mapping, alignment_score = align_global(csp_seq_str, occlusion_seq_str)
    
    print(f"[VISUALIZE ALIGN] Alignment score: {alignment_score}")
    print(f"[VISUALIZE ALIGN] Mapped pairs: {len(mapping)}")
    
    # Create mapping from sequential position to PDB residue number
    position_map = {}
    for csp_pos, occlusion_pos in mapping:
        if csp_pos <= len(csp_positions) and occlusion_pos <= len(occlusion_positions):
            sequential_position = csp_positions[csp_pos - 1]  # Convert to 1-based
            pdb_residue_number = occlusion_positions[occlusion_pos - 1]  # Convert to 1-based
            position_map[sequential_position] = pdb_residue_number
    
    print(f"[VISUALIZE ALIGN] Created position mapping with {len(position_map)} entries")
    print(f"[VISUALIZE ALIGN] First 10 mappings: {list(position_map.items())[:10]}")
    
    return position_map


def plot_csp_classification_bars(
    results: List[CSPResult],
    binding_results: dict,
    out_png: str,
    title: Optional[str] = None,
    significance_field: str = 'significant',
    include_numeric_residue_ticks: bool = False,
) -> None:
    """
    Create a bar plot showing CSP values per residue, color-coded by classification (TP/FP/TN/FN).
    
    Args:
        results: List of CSPResult objects from CSP analysis
        binding_results: Results dictionary providing residue annotations (SASA occlusion or interactions)
        out_png: Path to save the PNG file
        title: Optional title for the plot
        significance_field: Which significance field to use ('significant', 'significant_1sd', 'significant_2sd')
    """
    if not _HAS_PLT:
        return
    
    # Create sequence alignment mapping
    position_map = create_sequence_alignment_map_from_results(results, binding_results)
    
    # Create occlusion lookup by PDB residue number
    binding_lookup, dataset_type = _build_binding_lookup(binding_results)

    # Filter results with valid CSP values
    valid_results = [r for r in results if r.csp_A is not None]
    if not valid_results:
        return
    
    # Ensure we're using non-interactive backend
    import matplotlib
    matplotlib.use('Agg', force=True)
    
    # Extract data for plotting
    residue_numbers = []
    csp_values = []
    classifications = []
    
    # Color scheme
    colors = _get_classification_colors()
    
    for r in valid_results:
        residue_numbers.append(r.holo_index)
        csp_values.append(r.csp_A)
        
        # Get corresponding PDB residue number from alignment
        pdb_residue_number = position_map.get(r.holo_index)
        if pdb_residue_number is None:
            print(f"[VISUALIZE WARNING] No alignment found for sequential position {r.holo_index}")
            # Default to non-binding if no alignment found
            is_binding = False
        else:
            is_binding = binding_lookup.get(pdb_residue_number, False)

        # Get significance value based on the specified field
        is_significant = getattr(r, significance_field, False)

        # Classify based on significance and binding/site annotation
        if is_significant and is_binding:
            classifications.append('TP')  # True Positive
        elif is_significant and not is_binding:
            classifications.append('FP')  # False Positive
        elif not is_significant and not is_binding:
            classifications.append('TN')  # True Negative
        else:  # not is_significant and is_binding
            classifications.append('FN')  # False Negative
    
    # Calculate significance threshold (minimum significant CSP value)
    significant_csps = [r.csp_A for r in valid_results if getattr(r, significance_field, False)]
    threshold = min(significant_csps) if significant_csps else 0.0
    
    # Create the plot
    plt.figure(figsize=(14, 8))
    
    # Create color array for bars
    bar_colors = [colors[cls] for cls in classifications]
    
    # Create bar plot
    bars = plt.bar(residue_numbers, csp_values, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Add threshold line
    plt.axhline(y=threshold, color='black', linestyle='--', linewidth=2, alpha=0.8, label=f'Threshold = {threshold:.3f}')
    
    # Count classifications for legend
    tp_count = classifications.count('TP')
    fp_count = classifications.count('FP')
    tn_count = classifications.count('TN')
    fn_count = classifications.count('FN')
    
    # Create custom legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=colors['TP'], alpha=0.8, label=f'(TP) Sig. CSP in Binding Site ({tp_count})'),
        plt.Rectangle((0,0),1,1, facecolor=colors['FP'], alpha=0.8, label=f'(FP) Sig. CSP -- Allosteric ({fp_count})'),
        plt.Rectangle((0,0),1,1, facecolor=colors['TN'], alpha=0.8, label=f'(TN) low CSP -- Allosteric ({tn_count})'),
        plt.Rectangle((0,0),1,1, facecolor=colors['FN'], alpha=0.8, label=f'(FN) low CSP in Binding Site ({fn_count})')
    ]
    
    plt.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=28,
    )
    
    # Set x-axis ticks with amino acid sequence (letters only on first row).
    ax_main = plt.gca()
    aa_labels = [str(r.holo_aa)[0] if r.holo_aa else '' for r in valid_results]
    ax_main.set_xticks(residue_numbers)
    ax_main.set_xticklabels(aa_labels, fontsize=14)
    ax_main.tick_params(axis='x', which='major', pad=2, labelsize=14)
    ax_main.tick_params(axis='y', which='major', labelsize=22)
    
    if include_numeric_residue_ticks:
        # Create three-row residue index system
        residue_indices = [r.holo_index for r in valid_results]
        ones_ticks, tens_ticks, hundreds_ticks = create_three_row_residue_ticks(residue_indices)

        # Add secondary x-axes for digit-place rows.
        # Only create rows that will actually be populated to avoid default numeric ticks.
        ax2 = ax_main.secondary_xaxis('bottom')

        # Row 1: Ones place (closest to main axis)
        ax2.set_xticks([pos for pos, _ in ones_ticks])
        ax2.set_xticklabels([label for _, label in ones_ticks], fontsize=9)
        ax2.spines['bottom'].set_position(('outward', 18))
        ax2.tick_params(axis='x', which='major', pad=1, labelsize=9)

        # Row 2: Tens place (middle)
        if tens_ticks:
            ax3 = ax_main.secondary_xaxis('bottom')
            ax3.set_xticks([pos for pos, _ in tens_ticks])
            ax3.set_xticklabels([label for _, label in tens_ticks], fontsize=9)
            ax3.spines['bottom'].set_position(('outward', 36))
            ax3.tick_params(axis='x', which='major', pad=1, labelsize=9)

        # Row 3: Hundreds place (furthest from main axis)
        if hundreds_ticks:
            ax4 = ax_main.secondary_xaxis('bottom')
            ax4.set_xticks([pos for pos, _ in hundreds_ticks])
            ax4.set_xticklabels([label for _, label in hundreds_ticks], fontsize=9)
            ax4.spines['bottom'].set_position(('outward', 54))
            ax4.tick_params(axis='x', which='major', pad=1, labelsize=9)
    
    # Add grid for better readability
    plt.grid(True, alpha=0.3, axis='y')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
    
    print(f"[DEBUG] CSP classification bar plot saved: {out_png}")
    print(f"[DEBUG] Classification counts - TP: {tp_count}, FP: {fp_count}, TN: {tn_count}, FN: {fn_count}")


def plot_csp_histogram(results: List[CSPResult], out_png: str, title: Optional[str] = None) -> None:
    """Create a histogram of CSP values (original functionality)."""
    if not _HAS_PLT:
        return
    vals = [r.csp_A for r in results if r.csp_A is not None]
    if not vals:
        return
    
    # Ensure we're using non-interactive backend
    import matplotlib
    matplotlib.use('Agg', force=True)
    
    plt.figure(figsize=(6, 3))
    plt.hist(vals, bins=20, color="#4C78A8", edgecolor="white")
    plt.xlabel("CSP (A)")
    plt.ylabel("Count")
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _build_per_atom_classification_data(
    results: List[CSPResult],
    binding_results: dict,
    atom: str,
    threshold_config: Optional[Dict] = None,
) -> Optional[Dict[str, object]]:
    """
    Helper to compute per-atom deltas, thresholds, and TP/FP/TN/FN labels.
    """
    atom_upper = atom.upper()
    atom_stats = compute_atom_deltas_with_offset(results, atom_upper, threshold_config)
    if not atom_stats.per_residue:
        return None

    # Sequence alignment and binding lookup
    position_map = create_sequence_alignment_map_from_results(results, binding_results)
    binding_lookup, dataset_type = _build_binding_lookup(binding_results)

    residue_numbers: List[int] = []
    values: List[float] = []
    classifications: List[str] = []

    for entry in atom_stats.per_residue:
        holo_index = int(entry["holo_index"])
        delta_abs = float(entry["delta_abs"])
        significant = bool(entry.get("significant", False))

        pdb_residue_number = position_map.get(holo_index)
        if pdb_residue_number is None:
            is_binding = False
        else:
            is_binding = binding_lookup.get(pdb_residue_number, False)

        residue_numbers.append(holo_index)
        values.append(delta_abs)

        if significant and is_binding:
            classifications.append("TP")
        elif significant and not is_binding:
            classifications.append("FP")
        elif not significant and not is_binding:
            classifications.append("TN")
        else:
            classifications.append("FN")

    if not residue_numbers:
        return None

    return {
        "atom": atom_upper,
        "residue_numbers": residue_numbers,
        "values": values,
        "classifications": classifications,
        "offset": atom_stats.offset,
        "cutoff": atom_stats.cutoff,
        "dataset_type": dataset_type,
    }


def plot_per_atom_classification_panels(
    results_hn: List[CSPResult],
    results_ca: Optional[List[CSPResult]],
    binding_results: dict,
    out_png: str,
    title: Optional[str] = None,
    threshold_config: Optional[Dict] = None,
    pdb_id: Optional[str] = None,
    structure_pdb_path: Optional[str] = None,
    receptor_chain: Optional[str] = None,
) -> None:
    """
    Create a multi-panel bar plot showing per-atom perturbation classifications
    (H-only, N-only, CA-only) using TP/FP/TN/FN coloring.
    """
    if not _HAS_PLT:
        return

    # Ensure non-interactive backend
    import matplotlib
    matplotlib.use('Agg', force=True)

    atom_data_list: List[Dict[str, object]] = []

    # H and N from H/N CSP results
    for atom in ("H", "N"):
        data = _build_per_atom_classification_data(results_hn, binding_results, atom, threshold_config)
        if data is not None and data.get("values"):
            atom_data_list.append(data)

    # CA from CA-inclusive results, if available
    if results_ca:
        ca_data = _build_per_atom_classification_data(results_ca, binding_results, "CA", threshold_config)
        if ca_data is not None and ca_data.get("values"):
            atom_data_list.append(ca_data)

    if not atom_data_list:
        return

    n_panels = len(atom_data_list)
    fig, axes = plt.subplots(n_panels, 1, figsize=(14, 4 * n_panels), sharex=True)
    if n_panels == 1:
        axes_list = [axes]
    else:
        axes_list = list(axes)

    colors = _get_classification_colors()

    # Use residue indices from the first panel for tick helpers
    all_residue_indices: List[int] = list(atom_data_list[0]["residue_numbers"])  # type: ignore
    ss_labels = _resolve_secondary_structure_by_holo_index(
        results_hn,
        binding_results,
        all_residue_indices,
        pdb_id=pdb_id,
        structure_pdb_path=structure_pdb_path,
        receptor_chain=receptor_chain,
    )
    ss_lookup = {idx: ss for idx, ss in zip(all_residue_indices, ss_labels)}

    for ax, atom_data in zip(axes_list, atom_data_list):
        atom_label = str(atom_data["atom"])
        residue_numbers = list(atom_data["residue_numbers"])  # type: ignore
        values = list(atom_data["values"])  # type: ignore
        classifications = list(atom_data["classifications"])  # type: ignore
        cutoff = float(atom_data["cutoff"])  # type: ignore
        bar_colors = [colors.get(cls, '#7f8c8d') for cls in classifications]
        ax.bar(residue_numbers, values, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        panel_ss_labels = [ss_lookup.get(idx, "C") for idx in residue_numbers]
        _draw_secondary_structure_track(ax, residue_numbers, panel_ss_labels, inside_axis=True)
        non_negative_ticks = [tick for tick in ax.get_yticks() if tick >= 0]
        ax.set_yticks(non_negative_ticks)

        if cutoff > 0.0:
            ax.axhline(y=cutoff, color='black', linestyle='--', linewidth=1.5, alpha=0.8,
                       label=f'Threshold = {cutoff:.3f}')

        tp_count = classifications.count('TP')
        fp_count = classifications.count('FP')
        tn_count = classifications.count('TN')
        fn_count = classifications.count('FN')

        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor=colors['TP'], alpha=0.8,
                          label=f'(TP) Sig. Δ{atom_label} in binding site ({tp_count})'),
            plt.Rectangle((0, 0), 1, 1, facecolor=colors['FP'], alpha=0.8,
                          label=f'(FP) Sig. Δ{atom_label} -- Allosteric ({fp_count})'),
            plt.Rectangle((0, 0), 1, 1, facecolor=colors['TN'], alpha=0.8,
                          label=f'(TN) low Δ{atom_label} -- Allosteric ({tn_count})'),
            plt.Rectangle((0, 0), 1, 1, facecolor=colors['FN'], alpha=0.8,
                          label=f'(FN) low Δ{atom_label} in binding site ({fn_count})'),
        ]

        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        ax.set_ylabel(f"|Δ{atom_label}| (ppm)", fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')

        precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
        recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        mcc_denominator = ((tp_count + fp_count) * (tp_count + fn_count) * (tn_count + fp_count) * (tn_count + fn_count)) ** 0.5
        mcc_score = ((tp_count * tn_count) - (fp_count * fn_count)) / mcc_denominator if mcc_denominator > 0 else 0.0

        metrics_text = (
            f'Precision: {precision:.3f}\n'
            f'Recall: {recall:.3f}\n'
            f'F1-Score: {f1_score:.3f}\n'
            f'MCC: {mcc_score:.3f}'
        )
        ax.text(
            0.02,
            0.98,
            metrics_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=9,
        )

        ax.set_title(f"{atom_label}-only perturbation classification", fontsize=13, fontweight='bold')

    # Shared x-axis labelling on the bottom panel
    bottom_ax = axes_list[-1]
    if all_residue_indices:
        bottom_ax.set_xticks(all_residue_indices)
        bottom_ax.set_xticklabels(
            [r.holo_aa for r in results_hn if r.holo_index in all_residue_indices],
            fontsize=8,
        )

        ones_ticks, tens_ticks, hundreds_ticks = create_three_row_residue_ticks(all_residue_indices)
        ax2 = bottom_ax.secondary_xaxis('bottom')
        ax3 = bottom_ax.secondary_xaxis('bottom')
        ax4 = bottom_ax.secondary_xaxis('bottom')

        ax2.set_xticks([pos for pos, _ in ones_ticks])
        ax2.set_xticklabels([label for _, label in ones_ticks], fontsize=5)
        ax2.spines['bottom'].set_position(('outward', 15))

        if tens_ticks:
            ax3.set_xticks([pos for pos, _ in tens_ticks])
            ax3.set_xticklabels([label for _, label in tens_ticks], fontsize=5)
            ax3.spines['bottom'].set_position(('outward', 30))

        if hundreds_ticks:
            ax4.set_xticks([pos for pos, _ in hundreds_ticks])
            ax4.set_xticklabels([label for _, label in hundreds_ticks], fontsize=5)
            ax4.spines['bottom'].set_position(('outward', 45))
            ax4.set_xlabel('Residue', fontsize=11)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
