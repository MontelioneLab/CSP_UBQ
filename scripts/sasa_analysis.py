"""
SASA-based backbone occlusion analysis for protein-peptide complexes.

This module provides functionality to identify protein residues whose backbone NH groups
are occluded by peptide binding using Solvent Accessible Surface Area (SASA) calculations.
"""

import os
import numpy as np
import mdtraj as md
import tempfile
import shutil
import csv
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# Import utility functions
from util import get_longest_chain, get_shortest_chain

# MDTraj shrake_rupley returns SASA in nm²; convert to Å² for threshold comparison and output
NM2_TO_ANG2 = 100.0  # 1 nm² = 100 Å²


def extract_models_from_pdb(pdb_path: str) -> List[str]:
    """
    Extract individual models from a multi-model PDB file.
    
    Args:
        pdb_path: Path to the multi-model PDB file
        
    Returns:
        List of paths to individual model PDB files
    """
    model_files = []
    current_model = None
    model_count = 0
    
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    # Collect header lines (before first MODEL)
    header_lines = []
    for line in lines:
        if line.startswith('MODEL'):
            break
        if not line.startswith('ENDMDL'):
            header_lines.append(line)
    
    # Create temporary directory for model files
    temp_dir = tempfile.mkdtemp()
    
    for line in lines:
        if line.startswith('MODEL'):
            # Start of new model
            if current_model is not None:
                current_model.close()
            
            model_count += 1
            model_file_path = os.path.join(temp_dir, f'model_{model_count}.pdb')
            current_model = open(model_file_path, 'w')
            model_files.append(model_file_path)
            
            # Write header lines first
            current_model.write('\n'.join(header_lines))
            if header_lines:
                current_model.write('\n')
            
        elif line.startswith('ENDMDL'):
            # End of current model
            if current_model is not None:
                current_model.close()
                current_model = None
                
        elif current_model is not None:
            # Write line to current model file
            current_model.write(line)
    
    # Close last model if still open
    if current_model is not None:
        current_model.close()
    
    return model_files


def find_medoid_model_from_pdb(pdb_path: str) -> str:
    """
    Find the medoid model from a multi-model PDB file.
    
    Args:
        pdb_path: Path to the multi-model PDB file
        
    Returns:
        Path to the medoid model PDB file
    """
    try:
        # First check if it's a multi-model PDB
        with open(pdb_path, 'r') as f:
            has_models = any(line.startswith('MODEL') for line in f)
        
        if not has_models:
            print("[DEBUG] Single model PDB file, using as-is")
            return pdb_path
        
        print("[DEBUG] Multi-model PDB detected, extracting models...")
        model_files = extract_models_from_pdb(pdb_path)
        print(f"[DEBUG] Extracted {len(model_files)} models")
        
        if len(model_files) == 1:
            print("[DEBUG] Only one model found, using it")
            return model_files[0]
        
        # Load all models and calculate RMSD matrix
        ensemble = []
        for model_file in model_files:
            try:
                traj = md.load_pdb(model_file)
                ensemble.append(traj)
            except Exception as e:
                print(f"[DEBUG] Warning: Could not load model {model_file}: {e}")
                continue
        
        if len(ensemble) < 2:
            print("[DEBUG] Not enough valid models for medoid calculation, using first model")
            return model_files[0]
        
        print(f"[DEBUG] Calculating RMSD matrix for {len(ensemble)} models...")
        
        # Calculate RMSD matrix
        n_structures = len(ensemble)
        rmsd_matrix = np.zeros((n_structures, n_structures))
        
        for i in range(n_structures):
            for j in range(i+1, n_structures):
                try:
                    rmsd_arr = md.rmsd(ensemble[i], ensemble[j])
                    rmsd_val = float(np.asarray(rmsd_arr).flat[0])
                    rmsd_matrix[i, j] = rmsd_val
                    rmsd_matrix[j, i] = rmsd_val
                except Exception as e:
                    print(f"[DEBUG] Warning: Could not calculate RMSD between models {i} and {j}: {e}")
                    rmsd_matrix[i, j] = 0.0
                    rmsd_matrix[j, i] = 0.0
        
        # Find medoid (model with minimum sum of RMSDs)
        sum_rmsd = np.sum(rmsd_matrix, axis=1)
        medoid_index = np.argmin(sum_rmsd)
        
        print(f"[DEBUG] Selected medoid model {medoid_index + 1} (sum RMSD: {sum_rmsd[medoid_index]:.3f})")
        
        return model_files[medoid_index]
        
    except Exception as e:
        print(f"[DEBUG] Error finding medoid model: {e}")
        print("[DEBUG] Falling back to original PDB file")
        return pdb_path


def compute_sasa_occlusion(holo_pdb_path: str, sasa_threshold: float = 0.0, 
                          receptor_chain_id: str = None, ligand_chain_id: str = None,
                          mode: str = 'standard', probe_radius: float = 0.0, 
                          nh_mode: str = 'sum') -> Dict:
    """
    Compute SASA-based occlusion analysis for a receptor-ligand complex.
    
    Identifies receptor residues whose backbone NH groups are occluded by the ligand
    by comparing SASA values between the full complex and receptor-only structure.
    
    Args:
        holo_pdb_path: Path to the holo PDB file (receptor-ligand complex)
        sasa_threshold: Minimum delta SASA (Å²) to consider NH group as occluded (delta_sasa > threshold)
        receptor_chain_id: Optional receptor chain ID (if None, auto-detect longest chain)
        ligand_chain_id: Optional ligand chain ID (if None, auto-detect shortest chain)
        mode: SASA computation mode ('standard', 'fast', or 'residue' for per-residue)
        probe_radius: Probe radius for SASA calculation (0.0 exposes NH better for CSP)
        nh_mode: NH calculation mode ('sum' for N+H, 'n_only' for N-only)
        
    Returns:
        Dictionary containing:
        - occluded_residues: List of residue numbers with occluded NH groups
        - delta_sasa_values: Array of delta SASA for all residues
        - receptor_chain: ID of receptor chain
        - ligand_chain: ID of ligand chain
        - n_occluded: Count of occluded residues
        - fraction_occluded: Fraction of total residues occluded
        - residue_info: List of dicts with detailed residue information
        - avg_percent_drop: Average percentage drop in SASA
        - avg_percent_burial: Average percentage burial
    """
    try:
        # Find medoid model if it's a multi-model PDB
        medoid_pdb_path = find_medoid_model_from_pdb(holo_pdb_path)
        print(f"[DEBUG] Using PDB file: {medoid_pdb_path}")
        
        # Load the holo structure
        holo_traj = md.load_pdb(medoid_pdb_path)
        print(f"[DEBUG] Loaded holo structure with {holo_traj.n_atoms} atoms, {holo_traj.n_residues} residues")
        
        # Identify chains
        if receptor_chain_id is None or ligand_chain_id is None:
            # Auto-detect chains
            receptor_chain = get_longest_chain(medoid_pdb_path)
            ligand_chain = get_shortest_chain(medoid_pdb_path)
            
            # Validation for ambiguous receptor chains
            large_chains = [c for c in holo_traj.topology.chains if c.n_residues > 50]
            if len(large_chains) > 1:
                chain_info = [(c.chain_id, c.n_residues) for c in large_chains]
                raise ValueError(f"Ambiguous receptor chain; specify manually. Large chains: {chain_info}")
        else:
            receptor_chain = receptor_chain_id
            ligand_chain = ligand_chain_id
        
        print(f"[DEBUG] Identified receptor chain: {receptor_chain}, ligand chain: {ligand_chain}")
        
        # Print sequence lengths for confirmation
        receptor_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == receptor_chain)
        ligand_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == ligand_chain)
        print(f"[DEBUG] Receptor chain length: {receptor_chain_obj.n_residues} residues")
        print(f"[DEBUG] Ligand chain length: {ligand_chain_obj.n_residues} residues")
        
        if receptor_chain == ligand_chain:
            raise ValueError("Only one chain found in structure - cannot perform occlusion analysis")
        
        # Get all chain IDs
        all_chains = set([atom.residue.chain.chain_id for atom in holo_traj.topology.atoms])
        print(f"[DEBUG] All chains in structure: {all_chains}")
        if receptor_chain not in all_chains or ligand_chain not in all_chains:
            raise ValueError(f"Chain identification failed. Found chains: {all_chains}")
        
        # Create receptor-only structure by removing ligand chain
        try:
            receptor_chain_index = next(i for i, chain in enumerate(holo_traj.topology.chains) if chain.chain_id == receptor_chain)
        except StopIteration:
            raise ValueError(f"Receptor chain '{receptor_chain}' not found in topology. Available chains: {[chain.chain_id for chain in holo_traj.topology.chains]}")
        
        receptor_atoms = holo_traj.topology.select(f'chainid {receptor_chain_index}')
        print(f"[DEBUG] Selected {len(receptor_atoms)} atoms for receptor chain {receptor_chain} (index {receptor_chain_index})")
        receptor_traj = holo_traj.atom_slice(receptor_atoms, inplace=False)
        print(f"[DEBUG] Created receptor-only trajectory with {receptor_traj.n_atoms} atoms, {receptor_traj.n_residues} residues")
        
        # Compute SASA for both structures
        # Use probe radius of 1.4 Å (water molecule radius)
        import time
        start_time = time.time()
        
        if mode == 'residue':
            print(f"[DEBUG] Computing per-residue SASA (probe_radius={probe_radius})...")
            holo_sasa = md.shrake_rupley(holo_traj, probe_radius=probe_radius, mode='residue')
            receptor_sasa = md.shrake_rupley(receptor_traj, probe_radius=probe_radius, mode='residue')
            print(f"[DEBUG] Residue-level SASA arrays: holo shape {holo_sasa.shape}, receptor shape {receptor_sasa.shape}")
            print(f"[DEBUG] holo_sasa values (per residue): {holo_sasa}")
            print(f"[DEBUG] receptor_sasa values (per residue): {receptor_sasa}")
        elif mode == 'fast':
            print(f"[DEBUG] Computing SASA for holo structure (fast mode, probe_radius={probe_radius})...")
            holo_sasa = md.shrake_rupley(holo_traj, probe_radius=probe_radius, mode='residue')
            print(f"[DEBUG] holo_sasa values (per residue): {holo_sasa}")
            print(f"[DEBUG] Computing SASA for receptor-only structure (fast mode, probe_radius={probe_radius})...")
            receptor_sasa = md.shrake_rupley(receptor_traj, probe_radius=probe_radius, mode='residue')
            print(f"[DEBUG] receptor_sasa values (per residue): {receptor_sasa}")
        else:
            print(f"[DEBUG] Computing SASA for holo structure (probe_radius={probe_radius})...")
            holo_sasa = md.shrake_rupley(holo_traj, probe_radius=probe_radius)
            print(f"[DEBUG] holo_sasa values: {holo_sasa}")
            print(f"[DEBUG] Computing SASA for receptor-only structure (probe_radius={probe_radius})...")
            receptor_sasa = md.shrake_rupley(receptor_traj, probe_radius=probe_radius)
            print(f"[DEBUG] receptor_sasa values: {receptor_sasa}")
        
        computation_time = time.time() - start_time
        print(f"[DEBUG] SASA computation took {computation_time:.2f} seconds")
        print(f"[DEBUG] SASA arrays: holo shape {holo_sasa.shape}, receptor shape {receptor_sasa.shape}")
        # Warning for high probe radius
        if probe_radius > 1.0:
            print(f"[WARNING] High probe radius ({probe_radius} Å) may zero H SASA; consider probe=0.0 or n_only mode.")
        
        # Handle residue mode differently
        if mode == 'residue':
            print("[DEBUG] Processing residue-level SASA analysis...")
            
            # Get receptor residues and their indices
            receptor_residues = []
            receptor_res_indices = []
            for i, residue in enumerate(holo_traj.topology.residues):
                if residue.chain.chain_id == receptor_chain:
                    receptor_residues.append(residue.resSeq)
                    receptor_res_indices.append(i)
            
            print(f"[DEBUG] Found {len(receptor_residues)} receptor residues")
            
            # Calculate delta SASA for each residue
            delta_sasa_values = []
            occluded_residues = []
            residue_info = []
            
            for res_num in receptor_residues:
                # Find residue index in holo structure
                holo_res_idx = None
                for i, residue in enumerate(holo_traj.topology.residues):
                    if residue.resSeq == res_num and residue.chain.chain_id == receptor_chain:
                        holo_res_idx = i
                        break
                
                # Find residue index in receptor-only structure
                receptor_res_idx = None
                for i, residue in enumerate(receptor_traj.topology.residues):
                    if residue.resSeq == res_num:
                        receptor_res_idx = i
                        break
                
                if holo_res_idx is not None and receptor_res_idx is not None:
                    # Get residue SASA values (MDTraj returns nm²)
                    holo_res_sasa = holo_sasa[0, holo_res_idx]
                    receptor_res_sasa = receptor_sasa[0, receptor_res_idx]
                    
                    # Calculate delta SASA and convert to Å² for threshold comparison
                    delta_sasa_nm2 = receptor_res_sasa - holo_res_sasa
                    delta_sasa = delta_sasa_nm2 * NM2_TO_ANG2  # Å²
                    delta_sasa_values.append(delta_sasa)
                    
                    # Calculate percent burial (exposed_nh_sasa in Å²)
                    exposed_nh_sasa = 25.0  # Typical exposed N+H ~25 Å²
                    percent_burial = (delta_sasa / exposed_nh_sasa) * 100 if exposed_nh_sasa > 0 else 0
                    
                    # Determine occlusion based on sasa_threshold (Å²)
                    is_occluded = delta_sasa > sasa_threshold
                    if is_occluded:
                        occluded_residues.append(res_num)
                    
                    # Get residue name
                    residue_name = holo_traj.topology.residue(holo_res_idx).name
                    
                    residue_info.append({
                        'residue_number': res_num,
                        'residue_name': residue_name,
                        'chain_id': receptor_chain,
                        'delta_sasa': delta_sasa,
                        'percent_drop': 0.0,  # Not applicable for residue mode
                        'percent_burial': percent_burial,
                        'is_occluded': is_occluded,
                        'is_proline': residue_name == 'PRO',
                        'nh_mode': nh_mode,
                        'probe_radius': probe_radius,
                        'holo_res_sasa': holo_res_sasa,
                        'receptor_res_sasa': receptor_res_sasa,
                        'contact_confirmed': None  # Will be set later
                    })
                else:
                    print(f"[DEBUG] Warning: Could not find residue {res_num} in both structures")
            
            # Skip to contact validation
            print(f"[DEBUG] Processed {len(residue_info)} residues total")
            print(f"[DEBUG] Found {len(occluded_residues)} occluded residues")
            
        else:
            # Original atomic-level processing
            # Get receptor chain atoms in the original structure
            receptor_chain_atoms = holo_traj.topology.select(f'chainid {receptor_chain_index}')
            print(f"[DEBUG] Found {len(receptor_chain_atoms)} receptor chain atoms in original structure")
            print(f"[DEBUG] Total atoms in entire holo_traj: {holo_traj.n_atoms}")
            
            # Extract backbone N and H atoms for receptor chain
            backbone_n_atoms = []
            backbone_h_atoms = []
            residue_numbers = []
            residue_names = []
            
            for atom_idx in receptor_chain_atoms:
                atom = holo_traj.topology.atom(atom_idx)
                if atom.residue.chain.chain_id == receptor_chain:
                    if atom.name == 'N':
                        backbone_n_atoms.append(atom_idx)
                        residue_numbers.append(atom.residue.resSeq)
                        residue_names.append(atom.residue.name)
                    elif atom.name == 'H':
                        backbone_h_atoms.append(atom_idx)
            
            print(f"[DEBUG] Found {len(backbone_n_atoms)} backbone N atoms and {len(backbone_h_atoms)} backbone H atoms")
            print(f"[DEBUG] Residue numbers: {residue_numbers[:10]}... (showing first 10)")
            print(f"[DEBUG] Residue names: {residue_names[:10]}... (showing first 10)")
            
            # Warning about missing H atoms
            if len(backbone_h_atoms) < len(backbone_n_atoms) * 0.8:
                print(f"[WARNING] Missing H atoms in {len(backbone_n_atoms) - len(backbone_h_atoms)} residues; ensure PDB has hydrogens (e.g., via pdbfixer).")
            
            # Create dictionary for fast H atom lookup
            res_to_h = {holo_traj.topology.atom(h_idx).residue.resSeq: h_idx for h_idx in backbone_h_atoms}
            
            # Calculate delta SASA for each residue
            delta_sasa_values = []
            occluded_residues = []
            residue_info = []
            
            for i, (n_atom_idx, res_num, res_name) in enumerate(zip(backbone_n_atoms, residue_numbers, residue_names)):
                if i < 5:  # Debug first 5 residues
                    print(f"[DEBUG] Processing residue {i}: {res_name}{res_num}, N atom idx: {n_atom_idx}")
                
                # Skip proline residues (no backbone H atom)
                if res_name == 'PRO':
                    delta_sasa_values.append(0.0)
                    residue_info.append({
                        'residue_number': res_num,
                        'residue_name': res_name,
                        'chain_id': receptor_chain,
                        'delta_sasa': 0.0,
                        'percent_drop': 0.0,
                        'percent_burial': 0.0,
                        'is_occluded': False,
                        'is_proline': True,
                        'nh_mode': nh_mode,
                        'probe_radius': probe_radius,
                        'n_sasa_holo': 0.0,
                        'h_sasa_holo': 0.0,
                        'n_sasa_receptor': 0.0,
                        'h_sasa_receptor': 0.0,
                        'contact_confirmed': None
                    })
                    continue
            
                # Find corresponding H atom for this residue
                h_atom_idx = res_to_h.get(res_num)
                if h_atom_idx is not None:
                    h_atom = holo_traj.topology.atom(h_atom_idx)
                    if h_atom.residue.chain.chain_id != receptor_chain:
                        h_atom_idx = None
                
                if h_atom_idx is None:
                    if i < 5:
                        print(f"[DEBUG] No H atom found for residue {res_name}{res_num}")
                    # No H atom found, skip this residue
                    delta_sasa_values.append(0.0)
                    residue_info.append({
                        'residue_number': res_num,
                        'residue_name': res_name,
                        'chain_id': receptor_chain,
                        'delta_sasa': 0.0,
                        'percent_drop': 0.0,
                        'percent_burial': 0.0,
                        'is_occluded': False,
                        'is_proline': False,
                        'note': 'No backbone H atom found',
                        'nh_mode': nh_mode,
                        'probe_radius': probe_radius,
                        'n_sasa_holo': 0.0,
                        'h_sasa_holo': 0.0,
                        'n_sasa_receptor': 0.0,
                        'h_sasa_receptor': 0.0,
                        'contact_confirmed': None
                    })
                    continue
            
                if i < 5:
                    print(f"[DEBUG] Found H atom for residue {res_name}{res_num}: H atom idx: {h_atom_idx}")
                
                # Calculate SASA for N and H atoms
                n_sasa_holo = holo_sasa[0, n_atom_idx]
                h_sasa_holo = holo_sasa[0, h_atom_idx]
            
                if i < 5:
                    print(f"[DEBUG] Holo SASA - N: {n_sasa_holo:.3f}, H: {h_sasa_holo:.3f}")
                
                # Find corresponding atoms in receptor-only structure
                n_atom_receptor_idx = None
                h_atom_receptor_idx = None
                
                # Get all atoms in the sliced trajectory (only one chain at index 0)
                receptor_chain_atoms_sliced = list(range(receptor_traj.n_atoms))
                
                if i < 5:
                    print(f"[DEBUG] Searching in {len(receptor_chain_atoms_sliced)} atoms in sliced trajectory")
                
                for atom_idx in receptor_chain_atoms_sliced:
                    atom = receptor_traj.topology.atom(atom_idx)
                    if atom.name == 'N' and atom.residue.resSeq == res_num:
                        n_atom_receptor_idx = atom_idx
                    elif atom.name == 'H' and atom.residue.resSeq == res_num:
                        h_atom_receptor_idx = atom_idx
                
                if n_atom_receptor_idx is None or h_atom_receptor_idx is None:
                    if i < 5:
                        print(f"[DEBUG] Could not find atoms in receptor-only structure: N={n_atom_receptor_idx}, H={h_atom_receptor_idx}")
                    delta_sasa_values.append(0.0)
                    residue_info.append({
                        'residue_number': res_num,
                        'residue_name': res_name,
                        'chain_id': receptor_chain,
                        'delta_sasa': 0.0,
                        'percent_drop': 0.0,
                        'percent_burial': 0.0,
                        'is_occluded': False,
                        'is_proline': False,
                        'note': 'Could not find atoms in receptor-only structure',
                        'nh_mode': nh_mode,
                        'probe_radius': probe_radius,
                        'n_sasa_holo': 0.0,
                        'h_sasa_holo': 0.0,
                        'n_sasa_receptor': 0.0,
                        'h_sasa_receptor': 0.0,
                        'contact_confirmed': None
                    })
                    continue
            
                n_sasa_receptor = receptor_sasa[0, n_atom_receptor_idx]
                h_sasa_receptor = receptor_sasa[0, h_atom_receptor_idx]
                
                if i < 5:
                    print(f"[DEBUG] Receptor SASA - N: {n_sasa_receptor:.3f}, H: {h_sasa_receptor:.3f}")
                
                # Calculate delta SASA (receptor - holo); MDTraj returns nm²
                if nh_mode == 'n_only':
                    delta_sasa_nm2 = n_sasa_receptor - n_sasa_holo
                else:  # 'sum' mode
                    delta_sasa_nm2 = (n_sasa_receptor + h_sasa_receptor) - (n_sasa_holo + h_sasa_holo)
                delta_sasa = delta_sasa_nm2 * NM2_TO_ANG2  # Convert to Å²
                delta_sasa_values.append(delta_sasa)
                
                # Calculate percent drop (total_receptor_sasa in nm², delta_sasa_nm2 for ratio)
                if nh_mode == 'n_only':
                    total_receptor_sasa = n_sasa_receptor
                else:
                    total_receptor_sasa = n_sasa_receptor + h_sasa_receptor
                percent_drop = (delta_sasa_nm2 / total_receptor_sasa) * 100 if total_receptor_sasa > 0 else 0
                
                # Calculate percent burial (exposed_nh_sasa in Å²)
                exposed_nh_sasa = 25.0  # Typical exposed N+H ~25 Å²
                percent_burial = (delta_sasa / exposed_nh_sasa) * 100 if exposed_nh_sasa > 0 else 0
                
                # Determine occlusion based on sasa_threshold (Å²)
                is_occluded = delta_sasa >= sasa_threshold
                
                if i < 5:
                    print(f"[DEBUG] Delta SASA: {delta_sasa:.3f}, Percent drop: {percent_drop:.1f}%, Percent burial: {percent_burial:.1f}%")
                    print(f"[DEBUG] Occlusion threshold: {sasa_threshold:.3f} Å², Is occluded: {is_occluded}")
                if is_occluded:
                    occluded_residues.append(res_num)
                
                residue_info.append({
                    'residue_number': res_num,
                    'residue_name': res_name,
                    'chain_id': receptor_chain,
                    'delta_sasa': delta_sasa,
                    'percent_drop': percent_drop,
                    'percent_burial': percent_burial,
                    'is_occluded': is_occluded,
                    'is_proline': False,
                    'nh_mode': nh_mode,
                    'probe_radius': probe_radius,
                    'n_sasa_holo': n_sasa_holo,
                    'h_sasa_holo': h_sasa_holo,
                    'n_sasa_receptor': n_sasa_receptor,
                    'h_sasa_receptor': h_sasa_receptor,
                    'contact_confirmed': None  # Will be set later
                })
        
        print(f"[DEBUG] Processed {len(residue_info)} residues total")
        print(f"[DEBUG] Found {len(occluded_residues)} occluded residues")
        
        # Contact validation for occluded residues
        if occluded_residues:
            print("[DEBUG] Performing contact validation for occluded residues...")
            try:
                ligand_chain_index = next(i for i, chain in enumerate(holo_traj.topology.chains) if chain.chain_id == ligand_chain)
                ligand_atoms = holo_traj.topology.select(f'chainid {ligand_chain_index} and (name CA or name CB)')
                
                for res_num in occluded_residues:
                    receptor_nh_atoms = holo_traj.topology.select(f'residue {res_num} and chainid {receptor_chain_index} and (name N or name H)')
                    
                    if len(receptor_nh_atoms) > 0 and len(ligand_atoms) > 0:
                        contacts = md.compute_distances(holo_traj, [ligand_atoms, receptor_nh_atoms], periodic=False)[0, 0]
                        contact_confirmed = any(contact < 4.5 for contact in contacts)
                        
                        # Update residue_info with contact confirmation
                        for info in residue_info:
                            if info['residue_number'] == res_num:
                                info['contact_confirmed'] = contact_confirmed
                                break
                        
                        if contact_confirmed:
                            print(f"[DEBUG] Residue {res_num}: contact confirmed (<4.5 Å)")
                        else:
                            print(f"[DEBUG] Residue {res_num}: no close contacts found")
                            
            except Exception as e:
                print(f"[DEBUG] Warning: Contact validation failed: {e}")
                # Set contact_confirmed to None for all residues
                for info in residue_info:
                    info['contact_confirmed'] = None
        else:
            # Set contact_confirmed to None for all residues
            for info in residue_info:
                info['contact_confirmed'] = None
        
        # Calculate statistics
        n_occluded = len(occluded_residues)
        total_residues = len(residue_info)
        fraction_occluded = n_occluded / total_residues if total_residues > 0 else 0.0
        
        # Calculate average percent drop
        successful_residues = [info for info in residue_info if 'percent_drop' in info]
        avg_percent_drop = np.mean([info['percent_drop'] for info in successful_residues]) if successful_residues else 0.0
        
        # Calculate average percent burial
        burial_residues = [info for info in residue_info if 'percent_burial' in info]
        avg_percent_burial = np.mean([info['percent_burial'] for info in burial_residues]) if burial_residues else 0.0
        
        # Clean up temporary files if we created them
        if medoid_pdb_path != holo_pdb_path:
            try:
                temp_dir = os.path.dirname(medoid_pdb_path)
                shutil.rmtree(temp_dir)
                print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
            except Exception as e:
                print(f"[DEBUG] Warning: Could not clean up temporary files: {e}")
        
        return {
            'occluded_residues': occluded_residues,
            'delta_sasa_values': np.array(delta_sasa_values),
            'receptor_chain': receptor_chain,
            'ligand_chain': ligand_chain,
            'n_occluded': n_occluded,
            'fraction_occluded': fraction_occluded,
            'residue_info': residue_info,
            'sasa_threshold': sasa_threshold,
            'avg_percent_drop': avg_percent_drop,
            'avg_percent_burial': avg_percent_burial,
            'computation_time': computation_time,
            'probe_radius': probe_radius,
            'nh_mode': nh_mode
        }
        
    except Exception as e:
        print(f"Error in compute_sasa_occlusion: {e}")
        
        # Clean up temporary files if we created them
        try:
            if 'medoid_pdb_path' in locals() and medoid_pdb_path != holo_pdb_path:
                temp_dir = os.path.dirname(medoid_pdb_path)
                shutil.rmtree(temp_dir)
                print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
        except Exception as cleanup_e:
            print(f"[DEBUG] Warning: Could not clean up temporary files: {cleanup_e}")
        
        return {
            'occluded_residues': [],
            'delta_sasa_values': np.array([]),
            'receptor_chain': None,
            'ligand_chain': None,
            'n_occluded': 0,
            'fraction_occluded': 0.0,
            'residue_info': [],
            'sasa_threshold': sasa_threshold,
            'error': str(e)
        }


def compute_sasa_occlusion_nh_atom(
    holo_pdb_path: str,
    min_threshold: float = 0.5,
    receptor_chain_id: Optional[str] = None,
    ligand_chain_id: Optional[str] = None,
    probe_radius: float = 0.0,
) -> Dict:
    """
    Compute SASA-based occlusion using backbone N+H atom SASA (atom-level).

    Uses md.shrake_rupley with mode='atom' for both holo and receptor-only.
    For each receptor residue: sum SASA of backbone N and H atoms only (skip PRO for H).
    Per residue: delta_sasa = sum(SASA_receptor[N,H]) - sum(SASA_holo[N,H]).
    MDTraj returns nm²; delta_sasa is converted to Å² for comparison.
    is_occluded = delta_sasa > min_threshold (min_threshold in Å²).

    Returns the same structure as compute_sasa_occlusion.
    """
    try:
        medoid_pdb_path = find_medoid_model_from_pdb(holo_pdb_path)
        holo_traj = md.load_pdb(medoid_pdb_path)

        if receptor_chain_id is None or ligand_chain_id is None:
            receptor_chain = get_longest_chain(medoid_pdb_path)
            ligand_chain = get_shortest_chain(medoid_pdb_path)
            large_chains = [c for c in holo_traj.topology.chains if c.n_residues > 50]
            if len(large_chains) > 1:
                chain_info = [(c.chain_id, c.n_residues) for c in large_chains]
                raise ValueError(f"Ambiguous receptor chain; specify manually. Large chains: {chain_info}")
        else:
            receptor_chain = receptor_chain_id
            ligand_chain = ligand_chain_id

        if receptor_chain == ligand_chain:
            raise ValueError("Only one chain found in structure - cannot perform occlusion analysis")

        receptor_chain_index = next(
            i for i, chain in enumerate(holo_traj.topology.chains) if chain.chain_id == receptor_chain
        )
        receptor_atoms = holo_traj.topology.select(f"chainid {receptor_chain_index}")
        receptor_traj = holo_traj.atom_slice(receptor_atoms, inplace=False)

        holo_sasa = md.shrake_rupley(holo_traj, probe_radius=probe_radius, mode="atom")
        receptor_sasa = md.shrake_rupley(receptor_traj, probe_radius=probe_radius, mode="atom")

        receptor_chain_atoms = holo_traj.topology.select(f"chainid {receptor_chain_index}")
        backbone_n_atoms = []
        backbone_h_atoms = []
        residue_numbers = []
        residue_names = []

        for atom_idx in receptor_chain_atoms:
            atom = holo_traj.topology.atom(atom_idx)
            if atom.residue.chain.chain_id == receptor_chain:
                if atom.name == "N":
                    backbone_n_atoms.append(atom_idx)
                    residue_numbers.append(atom.residue.resSeq)
                    residue_names.append(atom.residue.name)
                elif atom.name == "H":
                    backbone_h_atoms.append(atom_idx)

        res_to_h = {holo_traj.topology.atom(h_idx).residue.resSeq: h_idx for h_idx in backbone_h_atoms}

        res_to_rec_n = {}
        res_to_rec_h = {}
        for i, atom in enumerate(receptor_traj.topology.atoms):
            if atom.name == "N":
                res_to_rec_n[atom.residue.resSeq] = i
            elif atom.name == "H":
                res_to_rec_h[atom.residue.resSeq] = i

        delta_sasa_values = []
        occluded_residues = []
        residue_info = []

        for n_atom_idx, res_num, res_name in zip(backbone_n_atoms, residue_numbers, residue_names):
            n_sasa_holo = holo_sasa[0, n_atom_idx]
            h_atom_idx = res_to_h.get(res_num)

            h_sasa_holo = holo_sasa[0, h_atom_idx] if h_atom_idx is not None else 0.0
            rec_n_idx = res_to_rec_n.get(res_num)
            rec_h_idx = res_to_rec_h.get(res_num)
            n_sasa_receptor = receptor_sasa[0, rec_n_idx] if rec_n_idx is not None else 0.0
            h_sasa_receptor = receptor_sasa[0, rec_h_idx] if rec_h_idx is not None else 0.0

            if res_name == "PRO":
                h_sasa_holo = 0.0
                h_sasa_receptor = 0.0

            nh_holo = n_sasa_holo + h_sasa_holo
            nh_receptor = n_sasa_receptor + h_sasa_receptor
            delta_sasa_nm2 = nh_receptor - nh_holo
            delta_sasa = delta_sasa_nm2 * NM2_TO_ANG2  # Convert to Å² for threshold and output
            is_occluded = delta_sasa > min_threshold

            delta_sasa_values.append(delta_sasa)
            if is_occluded:
                occluded_residues.append(res_num)

            exposed_nh_sasa = 25.0  # Å²
            percent_burial = (delta_sasa / exposed_nh_sasa) * 100 if exposed_nh_sasa > 0 else 0
            total_receptor_sasa = nh_receptor
            percent_drop = (delta_sasa_nm2 / total_receptor_sasa) * 100 if total_receptor_sasa > 0 else 0

            residue_info.append({
                "residue_number": res_num,
                "residue_name": res_name,
                "chain_id": receptor_chain,
                "delta_sasa": delta_sasa,
                "percent_drop": percent_drop,
                "percent_burial": percent_burial,
                "is_occluded": is_occluded,
                "is_proline": res_name == "PRO",
                "nh_mode": "sum",
                "probe_radius": probe_radius,
                "n_sasa_holo": n_sasa_holo,
                "h_sasa_holo": h_sasa_holo if h_atom_idx is not None else 0.0,
                "n_sasa_receptor": n_sasa_receptor,
                "h_sasa_receptor": h_sasa_receptor if h_atom_idx is not None else 0.0,
                "contact_confirmed": None,
            })

        n_occluded = len(occluded_residues)
        total_residues = len(residue_info)
        fraction_occluded = n_occluded / total_residues if total_residues > 0 else 0.0
        avg_percent_drop = np.mean([info["percent_drop"] for info in residue_info]) if residue_info else 0.0
        avg_percent_burial = np.mean([info["percent_burial"] for info in residue_info]) if residue_info else 0.0

        if medoid_pdb_path != holo_pdb_path:
            try:
                temp_dir = os.path.dirname(medoid_pdb_path)
                shutil.rmtree(temp_dir)
            except Exception:
                pass

        return {
            "occluded_residues": occluded_residues,
            "delta_sasa_values": np.array(delta_sasa_values),
            "receptor_chain": receptor_chain,
            "ligand_chain": ligand_chain,
            "n_occluded": n_occluded,
            "fraction_occluded": fraction_occluded,
            "residue_info": residue_info,
            "sasa_threshold": min_threshold,
            "avg_percent_drop": avg_percent_drop,
            "avg_percent_burial": avg_percent_burial,
            "probe_radius": probe_radius,
            "nh_mode": "sum",
        }
    except Exception as e:
        return {
            "occluded_residues": [],
            "delta_sasa_values": np.array([]),
            "receptor_chain": None,
            "ligand_chain": None,
            "n_occluded": 0,
            "fraction_occluded": 0.0,
            "residue_info": [],
            "sasa_threshold": min_threshold,
            "error": str(e),
        }


def write_occlusion_analysis_csv(residue_info: List[Dict], output_path: str) -> None:
    """
    Write occlusion analysis results to CSV file.
    
    Args:
        residue_info: List of residue information dictionaries from compute_sasa_occlusion
        output_path: Path to output CSV file
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'residue_number', 'residue_name', 'chain_id', 'delta_sasa', 'percent_drop', 'percent_burial',
            'is_occluded', 'is_proline', 'nh_mode', 'probe_radius',             'n_sasa_holo', 'h_sasa_holo',
            'n_sasa_receptor', 'h_sasa_receptor', 'contact_confirmed', 'note'
        ])
        
        for info in residue_info:
            writer.writerow([
                info['residue_number'],
                info['residue_name'],
                info['chain_id'],
                f"{info['delta_sasa']:.4f}",
                f"{info.get('percent_drop', 0.0):.4f}",
                f"{info.get('percent_burial', 0.0):.4f}",
                info['is_occluded'],
                info.get('is_proline', False),
                info.get('nh_mode', ''),
                f"{info.get('probe_radius', 0.0):.4f}",
                f"{info.get('n_sasa_holo', 0.0):.4f}",
                f"{info.get('h_sasa_holo', 0.0):.4f}",
                f"{info.get('n_sasa_receptor', 0.0):.4f}",
                f"{info.get('h_sasa_receptor', 0.0):.4f}",
                info.get('contact_confirmed', ''),
                info.get('note', '')
            ])


def get_occlusion_summary(sasa_results: Dict) -> str:
    """
    Generate a summary string of occlusion analysis results.
    
    Args:
        sasa_results: Results dictionary from compute_sasa_occlusion
        
    Returns:
        Formatted summary string
    """
    if 'error' in sasa_results:
        return f"Occlusion analysis failed: {sasa_results['error']}"
    
    # Calculate average delta SASA
    avg_delta = np.mean(sasa_results['delta_sasa_values']) if len(sasa_results['delta_sasa_values']) > 0 else 0.0
    
    summary = f"""
SASA Occlusion Analysis Summary:
- Receptor chain: {sasa_results['receptor_chain']}
- Ligand chain: {sasa_results['ligand_chain']}
- Total residues analyzed: {len(sasa_results['residue_info'])}
- Occluded residues: {sasa_results['n_occluded']}
- Fraction occluded: {sasa_results['fraction_occluded']:.2%}
- SASA threshold: {sasa_results['sasa_threshold']} Å²
- Probe radius: {sasa_results.get('probe_radius', 0.0):.1f} Å
- NH mode: {sasa_results.get('nh_mode', 'sum')}
- Average delta SASA: {avg_delta:.3f} Å²
- Average percent drop: {sasa_results.get('avg_percent_drop', 0.0):.1f}%
- Average percent burial: {sasa_results.get('avg_percent_burial', 0.0):.1f}%
- Computation time: {sasa_results.get('computation_time', 0.0):.2f} seconds
- Occluded residue numbers: {sasa_results['occluded_residues']}
"""
    return summary.strip()
