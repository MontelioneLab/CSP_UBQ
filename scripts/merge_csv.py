"""
Master CSV file generator with sequence alignment.

This module consolidates all CSV files for each target into a single master CSV file.
The CSP sequence (from csp_table.csv) serves as the ground truth reference, and all 
other CSV files are aligned to it.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple, Optional, Any

# Support running as a script or module
try:
    from .align import align_global
    from .config import paths
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.align import align_global
    from scripts.config import paths


def load_csp_reference(csp_table_path: str) -> Tuple[str, List[int], Dict[int, Dict[str, Any]]]:
    """
    Load CSP table and extract the reference sequence.
    
    Args:
        csp_table_path: Path to csp_table.csv
        
    Returns:
        Tuple of (reference_sequence, sequential_positions, csp_data)
        - reference_sequence: String of amino acids from holo_aa column
        - sequential_positions: List of sequential positions from holo_resi
        - csp_data: Dict mapping sequential_position to all row data
    """
    if not os.path.exists(csp_table_path):
        raise FileNotFoundError(f"CSP table not found: {csp_table_path}")
    
    reference_sequence = []
    sequential_positions = []
    csp_data = {}
    
    # Amino acid mapping from 1-letter to 3-letter codes
    aa_1to3 = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    
    with open(csp_table_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                holo_resi = int(row['holo_resi'])
                holo_aa = row['holo_aa'].strip()
                
                # Only include rows with valid amino acid data
                if holo_aa and holo_aa != 'P' and holo_aa != '':  # Skip prolines and empty
                    reference_sequence.append(holo_aa)
                    sequential_positions.append(holo_resi)
                    csp_data[holo_resi] = dict(row)
            except (ValueError, KeyError) as e:
                if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
                    print(f"Warning: Skipping invalid row in CSP table: {e}")
                continue
    
    ref_seq_str = ''.join(reference_sequence)
    
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Loaded CSP reference: {len(ref_seq_str)} residues")
        print(f"[MERGE] Reference sequence: {ref_seq_str[:50]}...")
    
    return ref_seq_str, sequential_positions, csp_data


def load_and_align_csv(csv_path: str, ref_sequence: str, ref_positions: List[int], 
                      csv_type: str = "unknown") -> Dict[int, Dict[str, Any]]:
    """
    Load other CSV files and align to reference sequence.
    
    Args:
        csv_path: Path to CSV file to align
        ref_sequence: Reference sequence from CSP table
        ref_positions: Sequential positions from CSP table
        csv_type: Type of CSV for logging purposes
        
    Returns:
        Dict mapping sequential_position to aligned row data
    """
    if not os.path.exists(csv_path):
        print(f"Warning: CSV file not found: {csv_path}")
        return {}
    
    # Extract sequence from CSV file
    csv_sequence = []
    csv_positions = []
    csv_data = {}
    
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # Try different column name patterns
                residue_number = None
                residue_name = None
                
                if 'residue_number' in row and 'residue_name' in row:
                    residue_number = int(row['residue_number'])
                    residue_name = row['residue_name'].strip()
                elif 'resi' in row and 'aa' in row:
                    residue_number = int(row['resi'])
                    residue_name = row['aa'].strip()
                else:
                    # Skip rows without recognizable sequence columns
                    continue
                
                if residue_name and residue_name != 'PRO':  # Skip prolines
                    # Convert 3-letter to 1-letter if needed
                    if len(residue_name) == 3:
                        aa_3to1 = {
                            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                        }
                        aa_letter = aa_3to1.get(residue_name, 'X')
                    else:
                        aa_letter = residue_name
                    
                    csv_sequence.append(aa_letter)
                    csv_positions.append(residue_number)
                    csv_data[residue_number] = dict(row)
                    
            except (ValueError, KeyError) as e:
                print(f"Warning: Skipping invalid row in {csv_type}: {e}")
                continue
    
    csv_seq_str = ''.join(csv_sequence)
    
    if not csv_seq_str:
        print(f"Warning: No valid sequence found in {csv_type}")
        return {}
    
    # Perform alignment
    aligned_ref, aligned_csv, mapping, score = align_global(ref_sequence, csv_seq_str)
    
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Aligned {csv_type}: {len(csv_seq_str)} residues, score: {score:.1f}")
        print(f"[MERGE] {csv_type} sequence: {csv_seq_str[:50]}...")
    
    # Map aligned positions back to sequential positions
    aligned_data = {}
    for ref_pos, csv_pos in mapping:
        if ref_pos <= len(ref_positions) and csv_pos <= len(csv_positions):
            sequential_position = ref_positions[ref_pos - 1]  # Convert to 1-based
            csv_position = csv_positions[csv_pos - 1]  # Convert to 1-based
            
            if csv_position in csv_data:
                aligned_data[sequential_position] = csv_data[csv_position]
    
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Mapped {len(aligned_data)} positions from {csv_type}")
    
    return aligned_data


def compute_classification(row: Dict[str, Any]) -> str:
    """
    Compute classification (TP, TN, FP, FN) for a residue based on CSP significance
    and ground truth using occluded_or_ca_or_interaction strategy.
    
    Args:
        row: Dictionary containing row data with all necessary fields
        
    Returns:
        Classification string: 'TP', 'TN', 'FP', 'FN', or '' (empty if data insufficient)
    """
    def _to_clean_str(value: Any) -> str:
        """Convert mixed scalar values to a stripped string safely."""
        if value is None:
            return ''
        if isinstance(value, str):
            return value.strip()
        return str(value).strip()

    def _to_bool(value: Any) -> bool:
        """Parse common bool-like values robustly (str/bool/numpy scalar)."""
        if value is None:
            return False
        if isinstance(value, bool):
            return value
        value_str = _to_clean_str(value).lower()
        return value_str in ('true', '1', 'yes', 'y', 't')

    # Extract significant value (prediction)
    significant_str = _to_clean_str(row.get('significant', ''))
    try:
        significant = int(significant_str) == 1
    except (ValueError, TypeError):
        # Missing or invalid significant value - cannot classify
        return ''
    
    # Determine ground truth using occluded_or_ca_or_interaction strategy
    # Check occlusion status
    is_occluded = _to_bool(row.get('is_occluded_occlusion', ''))
    
    # Check CA distance filter status
    passes_ca_filter = _to_bool(row.get('passes_filter_distance', ''))
    
    # Check interaction status
    has_hbond = _to_bool(row.get('has_hbond_interaction', ''))
    has_charge_complement = _to_bool(row.get('has_charge_complement_interaction', ''))
    has_pi_contact = _to_bool(row.get('has_pi_contact_interaction', ''))
    
    is_interacting = has_hbond or has_charge_complement or has_pi_contact
    
    # Check min_any_atom_distance < 2A (direct contact)
    is_sub_2A = False
    passes_sub_2A_raw = row.get('passes_sub_2A_filter_any_atom', '')
    passes_sub_2A_str = _to_clean_str(passes_sub_2A_raw)
    if passes_sub_2A_str:
        is_sub_2A = _to_bool(passes_sub_2A_raw)
    else:
        min_any_str = row.get('min_any_atom_distance', '')
        if min_any_str not in ('', None):
            try:
                is_sub_2A = float(min_any_str) < 2.0
            except (ValueError, TypeError):
                pass
    
    # Ground truth is True if ANY of occlusion, CA distance, interaction, or min_any_atom < 2A is True
    is_positive = is_occluded or passes_ca_filter or is_interacting or is_sub_2A
    
    # Classify based on significance and ground truth
    if significant and is_positive:
        return 'TP'  # True Positive
    elif significant and not is_positive:
        return 'FP'  # False Positive
    elif not significant and not is_positive:
        return 'TN'  # True Negative
    else:  # not significant and is_positive
        return 'FN'  # False Negative


def _residue_info_to_aligned_data(
    residue_info: list,
    ref_sequence: str,
    ref_positions: List[int],
    value_key: str,
    extra_keys: Optional[List[str]] = None,
) -> Dict[int, Dict[str, Any]]:
    """
    Convert residue_info from interaction_analysis to aligned_data format.
    Maps residue_number -> {value_key: value, ...} aligned to ref_positions.
    """
    if not residue_info:
        return {}
    keys_to_copy = [value_key] + (extra_keys or [])
    # Build csv_sequence from residue_info (residue_number, residue_name)
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    csv_sequence = []
    csv_positions = []
    csv_data = {}
    for info in residue_info:
        res_num = info.get('residue_number')
        res_name = info.get('residue_name', '').strip()
        if res_num is None or not res_name or res_name == 'PRO':
            continue
        aa_letter = aa_3to1.get(res_name, 'X')
        csv_sequence.append(aa_letter)
        csv_positions.append(res_num)
        csv_data[res_num] = {k: info.get(k) for k in keys_to_copy}
    csv_seq_str = ''.join(csv_sequence)
    if not csv_seq_str:
        return {}
    aligned_ref, aligned_csv, mapping, _ = align_global(ref_sequence, csv_seq_str)
    aligned_data = {}
    for ref_pos, csv_pos in mapping:
        if ref_pos <= len(ref_positions) and csv_pos <= len(csv_positions):
            seq_pos = ref_positions[ref_pos - 1]
            csv_pos_val = csv_positions[csv_pos - 1]
            if csv_pos_val in csv_data:
                aligned_data[seq_pos] = csv_data[csv_pos_val]
    return aligned_data


def _compute_missing_distance_filters(
    target_dir: str,
    ref_sequence: str,
    ref_positions: List[int],
    csp_data: Dict[int, Dict[str, Any]],
    aligned_data: Dict[str, Dict[int, Dict[str, Any]]],
    distance_threshold: float = 6.0,
) -> None:
    """
    Compute nn_distance and any_atom filters when CSV files don't exist.
    Writes the CSVs and adds results to aligned_data in-place.
    """
    try:
        from scripts.interaction_analysis import (
            compute_nn_distance_filter,
            compute_min_atom_distance_filter,
            write_nn_distance_csv,
            write_any_atom_distance_csv,
        )
    except ImportError:
        try:
            from .interaction_analysis import (
                compute_nn_distance_filter,
                compute_min_atom_distance_filter,
                write_nn_distance_csv,
                write_any_atom_distance_csv,
            )
        except ImportError:
            if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
                print("[MERGE] interaction_analysis not available, skipping N-N/any-atom computation")
            return

    holo_pdb = None
    for pos in ref_positions:
        if pos in csp_data:
            holo_pdb = (csp_data[pos].get('holo_pdb') or '').strip()
            break
    if not holo_pdb:
        return
    pdb_path = os.path.join(paths.pdb_cache_dir, f"{holo_pdb}.pdb")
    if not os.path.exists(pdb_path):
        if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
            print(f"[MERGE] PDB not found for {holo_pdb}, skipping N-N and any-atom distance computation")
        return

    need_nn = 'nn_distance' not in aligned_data or not aligned_data['nn_distance']
    need_any = 'any_atom' not in aligned_data or not aligned_data['any_atom']
    if not need_nn and not need_any:
        return

    if need_nn:
        nn_results = compute_nn_distance_filter(
            pdb_path,
            distance_threshold=distance_threshold,
            receptor_chain_id=None,
            ligand_chain_id=None,
        )
        if 'error' not in nn_results and nn_results.get('residue_info'):
            nn_csv = os.path.join(target_dir, "nn_distance_filter.csv")
            write_nn_distance_csv(nn_results['residue_info'], nn_csv)
            aligned_data['nn_distance'] = _residue_info_to_aligned_data(
                nn_results['residue_info'], ref_sequence, ref_positions, 'min_nn_distance'
            )

    if need_any:
        any_results = compute_min_atom_distance_filter(
            pdb_path,
            distance_threshold=distance_threshold,
            receptor_chain_id=None,
            ligand_chain_id=None,
        )
        if 'error' not in any_results and any_results.get('residue_info'):
            any_csv = os.path.join(target_dir, "any_atom_distance_filter.csv")
            write_any_atom_distance_csv(any_results['residue_info'], any_csv)
            aligned_data['any_atom'] = _residue_info_to_aligned_data(
                any_results['residue_info'], ref_sequence, ref_positions, 'min_any_atom_distance',
                extra_keys=['passes_sub_2A_filter']
            )


def merge_all_csv_files(target_dir: str, output_path: str) -> None:
    """
    Main merge function that consolidates all CSV files.
    
    Args:
        target_dir: Directory containing CSV files (e.g., outputs/1cf4)
        output_path: Path for output master CSV file
    """
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Starting merge for target: {target_dir}")
    
    # Load CSP table as reference
    csp_table_path = os.path.join(target_dir, "csp_table.csv")
    ref_sequence, ref_positions, csp_data = load_csp_reference(csp_table_path)
    
    if not ref_sequence:
        raise ValueError("No valid CSP reference sequence found")
    
    # Auto-detect CSV files to merge
    csv_files = []
    csv_types = {
        'occlusion_analysis.csv': 'occlusion',
        'interaction_filter.csv': 'interaction',
        'ca_distance_filter.csv': 'distance',
        'nn_distance_filter.csv': 'nn_distance',
        'any_atom_distance_filter.csv': 'any_atom',
        '1d_analysis.csv': 'one_d',
    }
    
    for filename in os.listdir(target_dir):
        if filename.endswith('.csv') and filename != 'master_alignment.csv':
            # Skip confusion matrix files, CSP table files, and grid search files
            if (filename.startswith('confusion_matrix') or 
                filename.startswith('csp_table') or 
                filename.startswith('offset_grid')):
                continue
            # Only process files that are in csv_types
            if filename in csv_types:
                csv_files.append(filename)
    
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Found CSV files to merge: {csv_files}")
    
    # Load and align each CSV file
    aligned_data = {}
    
    for csv_file in csv_files:
        csv_path = os.path.join(target_dir, csv_file)
        csv_type = csv_types[csv_file]  # Safe since we only added files from csv_types
        
        aligned_csv_data = load_and_align_csv(csv_path, ref_sequence, ref_positions, csv_type)
        aligned_data[csv_type] = aligned_csv_data

    # Compute N-N and any-atom distance filters when CSV files don't exist
    distance_threshold = 6.0
    if 'distance' in aligned_data and aligned_data['distance']:
        for pos, row in aligned_data['distance'].items():
            dt = row.get('distance_threshold')
            if dt is not None:
                try:
                    distance_threshold = float(dt)
                except (ValueError, TypeError):
                    pass
                break
    _compute_missing_distance_filters(
        target_dir, ref_sequence, ref_positions, csp_data, aligned_data, distance_threshold
    )
    
    # Merge all data into master CSV
    write_master_csv(csp_data, aligned_data, output_path, ref_positions)
    
    if os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes"):
        print(f"[MERGE] Master CSV written to: {output_path}")


def write_master_csv(csp_data: Dict[int, Dict[str, Any]], 
                    aligned_data: Dict[str, Dict[int, Dict[str, Any]]], 
                    output_path: str, ref_positions: List[int]) -> None:
    """
    Write the merged master CSV file.
    
    Args:
        csp_data: CSP table data
        aligned_data: Aligned data from other CSV files
        output_path: Output file path
        ref_positions: List of sequential positions
    """
    # Collect all unique column names
    all_columns = set()
    
    # Add CSP columns
    for row_data in csp_data.values():
        all_columns.update(row_data.keys())
    
    # Add columns from aligned data with prefixes
    for csv_type, data in aligned_data.items():
        for row_data in data.values():
            for col_name in row_data.keys():
                # Add prefix to avoid column name conflicts; for 1D analysis we
                # keep the original column names (they are already specific).
                if col_name in ['residue_number', 'residue_name', 'resi', 'aa']:
                    continue
                if csv_type == 'one_d':
                    all_columns.add(col_name)
                else:
                    prefixed_name = f"{col_name}_{csv_type}"
                    all_columns.add(prefixed_name)
    
    # Define column order
    primary_columns = [
        'sequential_position', 'pdb_residue_number', 'apo_bmrb', 'holo_bmrb', 'holo_pdb', 
        'chain', 'apo_resi', 'apo_aa', 'holo_resi', 'holo_aa'
    ]
    
    csp_columns = [
        'H_apo', 'N_apo', 'H_holo', 'N_holo', 'dH', 'dN', 'csp_A', 'significant', 'significant_1sd', 'significant_2sd'
    ]
    
    # 1D single-atom metrics (if present)
    one_d_columns = [
        'H_offset', 'N_offset', 'CA_offset',
        'CA_apo', 'CA_holo',
        'dH_1d', 'CSP_H_1d', 'z_H_1d', 'csp_H_1d_significant',
        'dN_1d', 'CSP_N_1d', 'z_N_1d', 'csp_N_1d_significant',
        'dCA_1d', 'CSP_CA_1d', 'z_CA_1d', 'csp_CA_1d_significant',
    ]
    
    # Add classification column after CSP columns
    classification_column = ['classification']

    # Explicit distance columns (min N-N and min any-atom) for cleaner downstream use
    distance_columns = ['min_n_distance', 'min_any_atom_distance']
    
    # Sort remaining columns (exclude distance_columns to avoid duplicates)
    remaining_columns = sorted(
        all_columns
        - set(primary_columns)
        - set(csp_columns)
        - set(classification_column)
        - set(one_d_columns)
        - set(distance_columns)
    )
    
    column_order = primary_columns + csp_columns + one_d_columns + classification_column + distance_columns + remaining_columns
    
    # Write CSV
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=column_order)
        writer.writeheader()
        
        for seq_pos in ref_positions:
            row = {}
            
            # Add sequential position
            row['sequential_position'] = seq_pos
            
            # Add CSP data
            if seq_pos in csp_data:
                csp_row = csp_data[seq_pos]
                row.update(csp_row)
                # Extract PDB residue number from CSP data if available
                if 'holo_resi' in csp_row:
                    row['pdb_residue_number'] = csp_row['holo_resi']
            
            # Add aligned data from other CSV files
            for csv_type, data in aligned_data.items():
                if seq_pos in data:
                    csv_row = data[seq_pos]
                    for col_name, value in csv_row.items():
                        if col_name in ['residue_number', 'residue_name', 'resi', 'aa']:
                            continue
                        if csv_type == 'one_d':
                            # For 1D analysis, keep original column names so they
                            # appear as direct fields (no prefix).
                            row[col_name] = value
                        else:
                            prefixed_name = f"{col_name}_{csv_type}"
                            row[prefixed_name] = value
            
            # Compute classification for this row
            row['classification'] = compute_classification(row)

            # Add min_n_distance and min_any_atom_distance (from nn_distance and any_atom aligned data)
            if 'nn_distance' in aligned_data and seq_pos in aligned_data['nn_distance']:
                val = aligned_data['nn_distance'][seq_pos].get('min_nn_distance')
                if val is not None:
                    row['min_n_distance'] = val
            if 'any_atom' in aligned_data and seq_pos in aligned_data['any_atom']:
                val = aligned_data['any_atom'][seq_pos].get('min_any_atom_distance')
                if val is not None:
                    row['min_any_atom_distance'] = val
            
            # Format float values to 4 decimal places
            for key, value in row.items():
                if isinstance(value, (float, int)) and not isinstance(value, bool):
                    row[key] = f"{value:.4f}"
            
            writer.writerow(row)


def main() -> None:
    """CLI interface for standalone usage."""
    parser = argparse.ArgumentParser(description="Merge CSV files with sequence alignment")
    parser.add_argument("--target", required=True, 
                       help="Target directory (e.g., outputs/1cf4)")
    parser.add_argument("--output", 
                       help="Output path (default: master_alignment.csv)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.target):
        print(f"Error: Target directory not found: {args.target}")
        sys.exit(1)
    
    output_path = args.output or os.path.join(args.target, "master_alignment.csv")
    
    try:
        merge_all_csv_files(args.target, output_path)
        print(f"Successfully created master CSV: {output_path}")
    except Exception as e:
        print(f"Error creating master CSV: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
