"""
Interaction analysis for protein-peptide complexes.

This module provides functionality to identify protein residues that interact with peptide residues
through hydrogen bonds or charge complementarity using structural analysis.
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
from sasa_analysis import find_medoid_model_from_pdb

try:
    from config import ca_distance_analysis
except Exception:
    try:
        from scripts.config import ca_distance_analysis
    except Exception:
        ca_distance_analysis = None  # type: ignore


# Define charged residue mappings
POSITIVE_RESIDUES = ['LYS', 'ARG', 'HIS']
NEGATIVE_RESIDUES = ['ASP', 'GLU']

# Aromatic ring definitions for pi-contact detection (each sub-list represents one ring)
AROMATIC_RING_DEFINITIONS = {
    'PHE': [
        ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    ],
    'TYR': [
        ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    ],
    'TRP': [
        ['CG', 'CD1', 'NE1', 'CE2', 'CD2'],
        ['CE2', 'CZ2', 'CH2', 'CZ3', 'CE3'],
    ],
    'HIS': [
        ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    ],
}

# Hydrogen bond criteria
HBOND_DISTANCE_THRESHOLD = 2.5  # Å (H-heavy atom distance)
HBOND_ANGLE_THRESHOLD = 120.0   # degrees

# Direct contact threshold for binding site definition (min inter-chain atom-atom distance)
# Uses config value; fallback if config unavailable
def _get_direct_contact_threshold() -> float:
    if ca_distance_analysis is not None:
        return getattr(ca_distance_analysis, 'direct_contact_threshold', 2.0)
    return 2.0


def _build_bond_map(topology: md.Topology) -> Dict[int, List[int]]:
    """
    Build adjacency mapping of atom index -> bonded atom indices.
    """
    bond_map: Dict[int, List[int]] = defaultdict(list)
    for bond in topology.bonds:
        atom1, atom2 = bond
        idx1 = atom1.index
        idx2 = atom2.index
        bond_map[idx1].append(idx2)
        bond_map[idx2].append(idx1)
    return bond_map


def compute_interaction_filter(holo_pdb_path: str, distance_threshold: float = 2.0,
                             receptor_chain_id: str = None, ligand_chain_id: str = None,
                             pi_distance_threshold: float = 6.0) -> Dict:
    """
    Compute interaction analysis for a receptor-ligand complex.
    
    Identifies receptor residues that interact with ligand residues through:
    1. Hydrogen bonds (using MDTraj's baker_hubbard function and custom validation)
    2. Charge complementarity (opposite charged residues in proximity)
    
    Args:
        holo_pdb_path: Path to the holo PDB file (receptor-ligand complex)
        distance_threshold: Distance threshold for charge complementarity (Å)
        receptor_chain_id: Optional receptor chain ID (if None, auto-detect longest chain)
        ligand_chain_id: Optional ligand chain ID (if None, auto-detect shortest chain)
        
    Returns:
        Dictionary containing:
        - residue_info: List of dicts with detailed residue interaction information
        - receptor_chain: ID of receptor chain
        - ligand_chain: ID of ligand chain
        - n_hbond_residues: Count of residues with hydrogen bonds
        - n_charge_residues: Count of residues with charge complementarity
        - n_pi_residues: Count of residues with pi contacts
        - n_multiple_residues: Count of residues engaging in multiple interaction types
        - n_interacting_residues: Count of residues with at least one interaction type
        - fraction_interacting: Fraction of total residues with interactions
        - hbond_pairs: List of hydrogen bond pairs
        - charge_pairs: List of charge complementary pairs
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
        
        print(f"[DEBUG] === CHAIN IDENTIFICATION ===")
        print(f"[DEBUG] Identified receptor chain: {receptor_chain}, ligand chain: {ligand_chain}")
        
        # Print sequence lengths for confirmation
        receptor_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == receptor_chain)
        ligand_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == ligand_chain)
        print(f"[DEBUG] Receptor chain length: {receptor_chain_obj.n_residues} residues")
        print(f"[DEBUG] Ligand chain length: {ligand_chain_obj.n_residues} residues")
        
        # Debug: Show all chains and their lengths
        print(f"[DEBUG] All chains in structure:")
        for chain in holo_traj.topology.chains:
            print(f"[DEBUG]   Chain {chain.chain_id}: {chain.n_residues} residues")
        
        if receptor_chain == ligand_chain:
            raise ValueError("Only one chain found in structure - cannot perform interaction analysis")
        
        # Get all chain IDs
        all_chains = set([atom.residue.chain.chain_id for atom in holo_traj.topology.atoms])
        print(f"[DEBUG] All chains in structure: {all_chains}")
        if receptor_chain not in all_chains or ligand_chain not in all_chains:
            raise ValueError(f"Chain identification failed. Found chains: {all_chains}")
        
        # Detect hydrogen bonds
        print("[DEBUG] Detecting hydrogen bonds...")
        hbond_pairs = detect_hydrogen_bonds(holo_traj, receptor_chain, ligand_chain)
        # for hbond_pair in hbond_pairs:
        #     print(f"[DEBUG] hbond_pair donor_chain: {hbond_pair['donor_chain']}")
        #     print(f"[DEBUG] hbond_pair acceptor_chain: {hbond_pair['acceptor_chain']}")
        #     print(f"[DEBUG] hbond_pair donor_atom: {hbond_pair['donor_atom']}")
        #     print(f"[DEBUG] hbond_pair acceptor_atom: {hbond_pair['acceptor_atom']}")
        #     print(f"[DEBUG] hbond_pair hydrogen_atom: {hbond_pair['hydrogen_atom']}")
        #     print(f"[DEBUG] hbond_pair distance: {hbond_pair['distance']}")
        # raise Exception("hbond_pairs")
        print(f"[DEBUG] Found {len(hbond_pairs)} hydrogen bonds")
        
        # Detect charge complementarity
        print("[DEBUG] Detecting charge complementarity...")
        charge_pairs = detect_charge_complementarity(holo_traj, receptor_chain, ligand_chain, distance_threshold)
        print(f"[DEBUG] Found {len(charge_pairs)} charge complementary pairs")
        
        # Detect aromatic pi-contact proximity
        print("[DEBUG] Detecting aromatic pi contacts...")
        pi_contact_pairs = detect_pi_contacts(
            holo_traj,
            receptor_chain,
            ligand_chain,
            pi_distance_threshold
        )
        print(f"[DEBUG] Found {len(pi_contact_pairs)} NH-aromatic pi contacts")
        
        # Analyze receptor residues
        print("[DEBUG] Analyzing receptor residues...")
        residue_info = analyze_protein_residues(
            holo_traj,
            receptor_chain,
            hbond_pairs,
            charge_pairs,
            pi_contact_pairs
        )
        # raise Exception("residue_info")
        # Calculate statistics
        n_hbond_residues = sum(1 for info in residue_info if info['has_hbond'])
        n_charge_residues = sum(1 for info in residue_info if info['has_charge_complement'])
        n_pi_residues = sum(1 for info in residue_info if info['has_pi_contact'])
        n_multiple_residues = sum(
            1 for info in residue_info
            if sum([info['has_hbond'], info['has_charge_complement'], info['has_pi_contact']]) >= 2
        )
        total_residues = len(residue_info)
        total_interacting = sum(
            1 for info in residue_info
            if info['has_hbond'] or info['has_charge_complement'] or info['has_pi_contact']
        )
        fraction_interacting = total_interacting / total_residues if total_residues > 0 else 0.0
        
        # Clean up temporary files if we created them
        if medoid_pdb_path != holo_pdb_path:
            try:
                temp_dir = os.path.dirname(medoid_pdb_path)
                shutil.rmtree(temp_dir)
                print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
            except Exception as e:
                print(f"[DEBUG] Warning: Could not clean up temporary files: {e}")
        
        return {
            'residue_info': residue_info,
            'receptor_chain': receptor_chain,
            'ligand_chain': ligand_chain,
            'n_hbond_residues': n_hbond_residues,
            'n_charge_residues': n_charge_residues,
            'n_multiple_residues': n_multiple_residues,
            'fraction_interacting': fraction_interacting,
            'hbond_pairs': hbond_pairs,
            'charge_pairs': charge_pairs,
            'pi_contact_pairs': pi_contact_pairs,
            'distance_threshold': distance_threshold,
            'pi_distance_threshold': pi_distance_threshold,
            'n_pi_residues': n_pi_residues,
            'n_interacting_residues': total_interacting
        }
        
    except Exception as e:
        print(f"Error in compute_interaction_filter: {e}")
        
        # Clean up temporary files if we created them
        try:
            if 'medoid_pdb_path' in locals() and medoid_pdb_path != holo_pdb_path:
                temp_dir = os.path.dirname(medoid_pdb_path)
                shutil.rmtree(temp_dir)
                print(f"[DEBUG] Cleaned up temporary files in {temp_dir}")
        except Exception as cleanup_e:
            print(f"[DEBUG] Warning: Could not clean up temporary files: {cleanup_e}")
        
        return {
            'residue_info': [],
            'receptor_chain': None,
            'ligand_chain': None,
            'n_hbond_residues': 0,
            'n_charge_residues': 0,
            'n_multiple_residues': 0,
            'fraction_interacting': 0.0,
            'hbond_pairs': [],
            'charge_pairs': [],
            'pi_contact_pairs': [],
            'distance_threshold': distance_threshold,
            'pi_distance_threshold': pi_distance_threshold,
            'n_pi_residues': 0,
            'n_interacting_residues': 0,
            'error': str(e)
        }


def detect_hydrogen_bonds(traj: md.Trajectory, receptor_chain: str, ligand_chain: str) -> List[Dict]:
    """
    Detect hydrogen bonds between receptor and ligand chains using explicit hydrogen atoms.
    
    Args:
        traj: MDTraj trajectory object
        receptor_chain: Receptor chain ID
        ligand_chain: Ligand chain ID
        
    Returns:
        List of hydrogen bond dictionaries with donor/acceptor information
    """
    hbond_pairs = []
    
    try:
        # Use MDTraj's baker_hubbard function
        hbonds = md.baker_hubbard(traj, freq=0.0, exclude_water=True)
        
        print(f"[DEBUG] MDTraj baker_hubbard returned: {type(hbonds)}")
        print(f"[DEBUG] MDTraj baker_hubbard length: {len(hbonds) if hasattr(hbonds, '__len__') else 'no length'}")
        
        if len(hbonds) > 0:
            print(f"[DEBUG] First element type: {type(hbonds[0])}")
            print(f"[DEBUG] First element length: {len(hbonds[0]) if hasattr(hbonds[0], '__len__') else 'no length'}")
            if len(hbonds[0]) > 0:
                print(f"[DEBUG] First hbond type: {type(hbonds[0][0])}")
                print(f"[DEBUG] First hbond shape: {hbonds[0][0].shape if hasattr(hbonds[0][0], 'shape') else 'no shape'}")
                print(f"[DEBUG] First hbond content: {hbonds[0][0]}")
        
        print(f"[DEBUG] MDTraj found {len(hbonds)} total hydrogen bonds")
        
        # Check if we have any hydrogen bonds to process
        if len(hbonds[0]) == 0:
            print(f"[DEBUG] No hydrogen bonds found by MDTraj, falling back to custom detection")
            raise Exception("No hydrogen bonds found by MDTraj")
        else:
            # MDTraj baker_hubbard returns a flat array, need to reshape it
            # Each hydrogen bond is 3 consecutive indices: donor, acceptor, hydrogen
            hbond_array = hbonds[0]
            if len(hbond_array) % 3 != 0:
                print(f"[DEBUG] Warning: Unexpected hbond array length {len(hbond_array)}, not divisible by 3")
                raise Exception("Unexpected hbond array length")
            else:
                # Reshape into groups of 3
                for i, triplet in enumerate(hbonds):
                    donor_idx, hydrogen_idx, acceptor_idx = triplet
                    print(f"[DEBUG] Processing hbond {i}: donor={donor_idx}, acceptor={acceptor_idx}, hydrogen={hydrogen_idx}")
                    
                    donor_atom = traj.topology.atom(donor_idx)
                    acceptor_atom = traj.topology.atom(acceptor_idx)
                    hydrogen_atom = traj.topology.atom(hydrogen_idx)
                    
                    donor_chain = donor_atom.residue.chain.chain_id
                    acceptor_chain = acceptor_atom.residue.chain.chain_id
                    hydrogen_chain = hydrogen_atom.residue.chain.chain_id

                    donor_residue = donor_atom.residue.resSeq
                    acceptor_residue = acceptor_atom.residue.resSeq
                    hydrogen_residue = hydrogen_atom.residue.resSeq

                    print(f"[DEBUG] Donor residue: index={donor_idx}, name={donor_atom.name}, element={donor_atom.element.symbol if donor_atom.element else 'Unknown'}, chain={donor_chain}, residue={donor_residue}")
                    print(f"[DEBUG] Acceptor residue: index={acceptor_idx}, name={acceptor_atom.name}, element={acceptor_atom.element.symbol if acceptor_atom.element else 'Unknown'}, chain={acceptor_chain}, residue={acceptor_residue}")
                    print(f"[DEBUG] Hydrogen residue: index={hydrogen_idx}, name={hydrogen_atom.name}, element={hydrogen_atom.element.symbol if hydrogen_atom.element else 'Unknown'}, chain={hydrogen_chain}, residue={hydrogen_residue}")
                    # raise Exception("atom types")
                    # Check if this is a cross-chain hydrogen bond
                    if (donor_chain == receptor_chain and acceptor_chain == ligand_chain) or \
                       (donor_chain == ligand_chain and acceptor_chain == receptor_chain):
                        
                        # Validate with explicit hydrogen atom
                        if hydrogen_idx < traj.n_atoms:
                            hydrogen_atom = traj.topology.atom(hydrogen_idx)
                            hydrogen_name = hydrogen_atom.name

                            # Calculate H-acceptor distance (should be ~1.8 Å for good H-bonds)
                            h_acceptor_dist = np.linalg.norm(traj.xyz[0, hydrogen_idx] - traj.xyz[0, acceptor_idx]) * 10.0 # Å to nm
                            
                            # Additional validation with custom criteria
                            if validate_hydrogen_bond(traj, donor_idx, acceptor_idx, hydrogen_idx):
                                hbond_info = {
                                    'donor_residue': donor_atom.residue.resSeq,
                                    'donor_chain': donor_chain,
                                    'donor_atom': donor_atom.name,
                                    'acceptor_residue': acceptor_atom.residue.resSeq,
                                    'acceptor_chain': acceptor_chain,
                                    'acceptor_atom': acceptor_atom.name,
                                    'hydrogen_atom': hydrogen_name,
                                    'distance': h_acceptor_dist  # H-acceptor distance
                                }
                                hbond_pairs.append(hbond_info)
                                print(f"[DEBUG]   Added H-bond: {donor_chain}{donor_atom.residue.resSeq}({donor_atom.name})-{hydrogen_name} -> {acceptor_chain}{acceptor_atom.residue.resSeq}({acceptor_atom.name}) - H-A dist: {h_acceptor_dist:.2f}Å")
                            else:
                                print(f"[DEBUG]   Rejected H-bond: {donor_chain}{donor_atom.residue.resSeq}({donor_atom.name})-{hydrogen_name} -> {acceptor_chain}{acceptor_atom.residue.resSeq}({acceptor_atom.name}) - H-A dist: {h_acceptor_dist:.2f}Å (failed validation)")
                        else:
                            print(f"[DEBUG]   Skipped H-bond: invalid hydrogen index {hydrogen_idx}")
                    
    except Exception as e:
        print(f"[DEBUG] Warning: Hydrogen bond detection failed: {e}")
        print(f"[DEBUG] Falling back to custom detection using explicit hydrogens")
        # Fallback to custom detection using explicit hydrogens
        raise Exception("Hydrogen bond detection failed")
    
    print(f"[DEBUG] Found {len(hbond_pairs)} cross-chain hydrogen bonds")
    return hbond_pairs


def validate_hydrogen_bond(traj: md.Trajectory, donor_idx: int, acceptor_idx: int, hydrogen_idx: int) -> bool:
    """
    Validate hydrogen bond using H-acceptor distance and angle criteria.
    
    Args:
        traj: MDTraj trajectory object
        donor_idx: Donor atom index
        acceptor_idx: Acceptor atom index
        hydrogen_idx: Hydrogen atom index
        
    Returns:
        True if hydrogen bond meets criteria
    """
    try:
        # Calculate H-acceptor distance (should be ~1.8 Å for good H-bonds)
        hydrogen_pos = traj.xyz[0, hydrogen_idx]
        acceptor_pos = traj.xyz[0, acceptor_idx]
        h_acceptor_dist = np.linalg.norm(hydrogen_pos - acceptor_pos)
        
        if h_acceptor_dist > HBOND_DISTANCE_THRESHOLD:
            return False
        
        # Calculate D-H···A angle
        donor_pos = traj.xyz[0, donor_idx]
        
        # Vector from hydrogen to donor
        vec_HD = donor_pos - hydrogen_pos
        # Vector from hydrogen to acceptor
        vec_HA = acceptor_pos - hydrogen_pos
        
        # Calculate angle at hydrogen
        norm_HD = np.linalg.norm(vec_HD)
        norm_HA = np.linalg.norm(vec_HA)
        
        if norm_HD == 0 or norm_HA == 0:
            return False
        
        cos_angle = np.dot(vec_HD, vec_HA) / (norm_HD * norm_HA)
        angle = np.arccos(np.clip(cos_angle, -1.0, 1.0)) * 180.0 / np.pi
        
        return angle >= HBOND_ANGLE_THRESHOLD
        
    except Exception:
        return False

def is_reasonable_hbond_residue_pair(protein_residue: str, peptide_residue: str) -> bool:
    """
    Check if two residues could reasonably form hydrogen bonds.
    
    Args:
        protein_residue: Protein residue name
        peptide_residue: Peptide residue name
        
    Returns:
        True if the residues could form hydrogen bonds
    """
    # Define residues that can participate in hydrogen bonds
    hbond_capable = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                     'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    
    return protein_residue in hbond_capable and peptide_residue in hbond_capable


def is_reasonable_hbond_pair(donor_atom: str, acceptor_atom: str) -> bool:
    """
    Check if a donor-acceptor pair is reasonable for hydrogen bonding.
    
    Args:
        donor_atom: Donor atom name
        acceptor_atom: Acceptor atom name
        
    Returns:
        True if the pair is reasonable for hydrogen bonding
    """
    # Define reasonable hydrogen bond patterns
    reasonable_pairs = [
        # Backbone-backbone
        ('N', 'O'),  # Peptide bond
        # Sidechain-backbone
        ('NE', 'O'), ('NH1', 'O'), ('NH2', 'O'), ('ND1', 'O'), ('ND2', 'O'), ('NE2', 'O'), ('NZ', 'O'),
        ('OG', 'O'), ('OH', 'O'),
        # Backbone-sidechain
        ('N', 'OE1'), ('N', 'OE2'), ('N', 'OD1'), ('N', 'OD2'), ('N', 'OG'), ('N', 'OG1'), ('N', 'OH'),
        # Sidechain-sidechain
        ('NE', 'OE1'), ('NE', 'OE2'), ('NE', 'OD1'), ('NE', 'OD2'), ('NE', 'OG'), ('NE', 'OG1'), ('NE', 'OH'),
        ('NH1', 'OE1'), ('NH1', 'OE2'), ('NH1', 'OD1'), ('NH1', 'OD2'), ('NH1', 'OG'), ('NH1', 'OG1'), ('NH1', 'OH'),
        ('NH2', 'OE1'), ('NH2', 'OE2'), ('NH2', 'OD1'), ('NH2', 'OD2'), ('NH2', 'OG'), ('NH2', 'OG1'), ('NH2', 'OH'),
        ('ND1', 'OE1'), ('ND1', 'OE2'), ('ND1', 'OD1'), ('ND1', 'OD2'), ('ND1', 'OG'), ('ND1', 'OG1'), ('ND1', 'OH'),
        ('ND2', 'OE1'), ('ND2', 'OE2'), ('ND2', 'OD1'), ('ND2', 'OD2'), ('ND2', 'OG'), ('ND2', 'OG1'), ('ND2', 'OH'),
        ('NE2', 'OE1'), ('NE2', 'OE2'), ('NE2', 'OD1'), ('NE2', 'OD2'), ('NE2', 'OG'), ('NE2', 'OG1'), ('NE2', 'OH'),
        ('NZ', 'OE1'), ('NZ', 'OE2'), ('NZ', 'OD1'), ('NZ', 'OD2'), ('NZ', 'OG'), ('NZ', 'OG1'), ('NZ', 'OH'),
        ('OG', 'OE1'), ('OG', 'OE2'), ('OG', 'OD1'), ('OG', 'OD2'), ('OG', 'OG'), ('OG', 'OG1'), ('OG', 'OH'),
        ('OH', 'OE1'), ('OH', 'OE2'), ('OH', 'OD1'), ('OH', 'OD2'), ('OH', 'OG'), ('OH', 'OG1'), ('OH', 'OH'),
    ]
    
    return (donor_atom, acceptor_atom) in reasonable_pairs


def detect_charge_complementarity(traj: md.Trajectory, receptor_chain: str, ligand_chain: str, 
                                distance_threshold: float) -> List[Dict]:
    """
    Detect charge complementary pairs between receptor and ligand chains.
    
    Args:
        traj: MDTraj trajectory object
        receptor_chain: Receptor chain ID
        ligand_chain: Ligand chain ID
        distance_threshold: Distance threshold for charge complementarity
        
    Returns:
        List of charge complementary pair dictionaries
    """
    charge_pairs = []
    
    # Get charged residues from each chain
    receptor_charged = []
    ligand_charged = []
    
    print(f"[DEBUG] === DETECTING CHARGE COMPLEMENTARITY ===")
    print(f"[DEBUG] Distance threshold: {distance_threshold} Å")
    
    for residue in traj.topology.residues:
        if residue.chain.chain_id == receptor_chain:
            if residue.name in POSITIVE_RESIDUES + NEGATIVE_RESIDUES:
                receptor_charged.append(residue)
                charge_type = 'positive' if residue.name in POSITIVE_RESIDUES else 'negative'
                print(f"[DEBUG]   Receptor charged residue: {residue.name}{residue.resSeq} ({charge_type})")
        elif residue.chain.chain_id == ligand_chain:
            if residue.name in POSITIVE_RESIDUES + NEGATIVE_RESIDUES:
                ligand_charged.append(residue)
                charge_type = 'positive' if residue.name in POSITIVE_RESIDUES else 'negative'
                print(f"[DEBUG]   Ligand charged residue: {residue.name}{residue.resSeq} ({charge_type})")
    
    print(f"[DEBUG] Found {len(receptor_charged)} charged receptor residues, {len(ligand_charged)} charged ligand residues")
    
    # Check for complementary charge pairs
    for receptor_res in receptor_charged:
        for ligand_res in ligand_charged:
            receptor_charge = 'positive' if receptor_res.name in POSITIVE_RESIDUES else 'negative'
            ligand_charge = 'positive' if ligand_res.name in POSITIVE_RESIDUES else 'negative'
            
            # Check if charges are complementary
            if receptor_charge != ligand_charge:
                # Calculate minimum distance between sidechain atoms
                min_distance = calculate_residue_distance(traj, receptor_res, ligand_res) * 10.0 # Å to nm
                
                print(f"[DEBUG]   Checking {receptor_res.name}{receptor_res.resSeq} ({receptor_charge}) <-> {ligand_res.name}{ligand_res.resSeq} ({ligand_charge}): {min_distance:.2f} Å")
                
                if min_distance <= distance_threshold:
                    charge_info = {
                        'receptor_residue': receptor_res.resSeq,
                        'receptor_chain': receptor_res.chain.chain_id,
                        'receptor_residue_name': receptor_res.name,
                        'receptor_charge': receptor_charge,
                        'ligand_residue': ligand_res.resSeq,
                        'ligand_chain': ligand_res.chain.chain_id,
                        'ligand_residue_name': ligand_res.name,
                        'ligand_charge': ligand_charge,
                        'distance': min_distance
                    }
                    charge_pairs.append(charge_info)
                    print(f"[DEBUG]     ✓ CHARGE PAIR ADDED: {receptor_res.name}{receptor_res.resSeq} <-> {ligand_res.name}{ligand_res.resSeq}")
                else:
                    print(f"[DEBUG]     ✗ Too far: {min_distance:.2f} Å > {distance_threshold} Å")
            else:
                print(f"[DEBUG]   Skipping same charge: {receptor_res.name}{receptor_res.resSeq} ({receptor_charge}) <-> {ligand_res.name}{ligand_res.resSeq} ({ligand_charge})")
    
    print(f"[DEBUG] Total charge complementary pairs found: {len(charge_pairs)}")
    # raise Exception("charge_pairs")
    return charge_pairs


def _compute_ring_centers(traj: md.Trajectory, residue: md.core.topology.Residue) -> List[Dict]:
    """
    Compute aromatic ring centers for a residue.
    
    Args:
        traj: MDTraj trajectory object
        residue: Residue to evaluate
        
    Returns:
        List of dictionaries containing ring center information
    """
    ring_definitions = AROMATIC_RING_DEFINITIONS.get(residue.name, [])
    if not ring_definitions:
        return []
    
    atom_lookup = {atom.name: atom for atom in residue.atoms}
    ring_centers = []
    
    for ring_index, atom_names in enumerate(ring_definitions):
        atom_indices = []
        missing_atom = False
        for atom_name in atom_names:
            atom = atom_lookup.get(atom_name)
            if atom is None:
                missing_atom = True
                break
            atom_indices.append(atom.index)
        
        if missing_atom or not atom_indices:
            continue
        
        # MDTraj stores coordinates in nanometers; convert to angstroms when needed later
        coords = traj.xyz[0, atom_indices, :]
        center = coords.mean(axis=0)
        
        ring_centers.append({
            'center': center,
            'residue': residue.resSeq,
            'residue_name': residue.name,
            'chain_id': residue.chain.chain_id,
            'atom_names': atom_names,
            'ring_index': ring_index
        })
    
    return ring_centers


def detect_pi_contacts(traj: md.Trajectory, receptor_chain: str, ligand_chain: str,
                       distance_threshold: float) -> List[Dict]:
    """
    Detect NH-to-aromatic pi contacts between receptor and ligand.
    
    Args:
        traj: MDTraj trajectory object
        receptor_chain: Receptor chain ID
        ligand_chain: Ligand chain ID
        distance_threshold: Distance threshold in angstroms
        
    Returns:
        List of pi-contact dictionaries with distance information
    """
    pi_contacts: Dict[Tuple[int, int, int], Dict] = {}
    
    # Pre-compute bond adjacency once
    bond_map = _build_bond_map(traj.topology)

    # Pre-compute aromatic ring centers for ligand residues
    ring_centers = []
    for residue in traj.topology.residues:
        if residue.chain.chain_id != ligand_chain:
            continue
        ring_centers.extend(_compute_ring_centers(traj, residue))
    
    if not ring_centers:
        print("[DEBUG] No aromatic rings detected in ligand chain")
        return []
    
    # Iterate over receptor residues to find NH groups
    for residue in traj.topology.residues:
        if residue.chain.chain_id != receptor_chain:
            continue
        
        n_atom = next((atom for atom in residue.atoms if atom.name == 'N'), None)
        if n_atom is None:
            continue
        
        # Identify hydrogens bonded to the backbone nitrogen
        hydrogen_atoms: List[md.core.topology.Atom] = []
        for bonded_idx in bond_map.get(n_atom.index, []):
            bonded_atom = traj.topology.atom(bonded_idx)
            if bonded_atom.element is not None and bonded_atom.element.symbol == 'H':
                hydrogen_atoms.append(bonded_atom)
        
        if not hydrogen_atoms:
            continue
        
        n_pos = traj.xyz[0, n_atom.index]
        
        for h_atom in hydrogen_atoms:
            h_pos = traj.xyz[0, h_atom.index]
            midpoint = (n_pos + h_pos) / 2.0
            
            for ring in ring_centers:
                distance = np.linalg.norm(midpoint - ring['center']) * 10.0  # nm -> Å
                if distance <= distance_threshold:
                    key = (residue.resSeq, ring['residue'], ring['ring_index'])
                    existing = pi_contacts.get(key)
                    
                    if existing is None or distance < existing['distance']:
                        pi_contacts[key] = {
                            'receptor_residue': residue.resSeq,
                            'receptor_chain': receptor_chain,
                            'receptor_residue_name': residue.name,
                            'ligand_residue': ring['residue'],
                            'ligand_chain': ring['chain_id'],
                            'ligand_residue_name': ring['residue_name'],
                            'distance': distance,
                            'ring_atoms': ','.join(ring['atom_names'])
                        }
    
    return list(pi_contacts.values())


def calculate_residue_distance(traj: md.Trajectory, res1: md.core.topology.Residue, 
                             res2: md.core.topology.Residue) -> float:
    """
    Calculate minimum distance between two residues using sidechain atoms.
    
    Args:
        traj: MDTraj trajectory object
        res1: First residue
        res2: Second residue
        
    Returns:
        Minimum distance between sidechain atoms
    """
    min_distance = float('inf')
    
    # Get sidechain atoms for each residue
    res1_atoms = [atom for atom in res1.atoms if atom.name not in ['N', 'CA', 'C', 'O', 'H']]
    res2_atoms = [atom for atom in res2.atoms if atom.name not in ['N', 'CA', 'C', 'O', 'H']]
    
    # If no sidechain atoms, use CA atoms
    if not res1_atoms:
        res1_atoms = [atom for atom in res1.atoms if atom.name == 'CA']
    if not res2_atoms:
        res2_atoms = [atom for atom in res2.atoms if atom.name == 'CA']
    
    # Calculate minimum distance
    for atom1 in res1_atoms:
        for atom2 in res2_atoms:
            distance = np.linalg.norm(traj.xyz[0, atom1.index] - traj.xyz[0, atom2.index])
            min_distance = min(min_distance, distance)
    
    return min_distance


def analyze_protein_residues(traj: md.Trajectory, receptor_chain: str, 
                           hbond_pairs: List[Dict], charge_pairs: List[Dict],
                           pi_contact_pairs: List[Dict]) -> List[Dict]:
    """
    Analyze receptor residues for interaction patterns.
    
    Args:
        traj: MDTraj trajectory object
        receptor_chain: Receptor chain ID
        hbond_pairs: List of hydrogen bond pairs
        charge_pairs: List of charge complementary pairs
        
    Returns:
        List of residue information dictionaries
    """
    residue_info = []
    
    # Get all receptor residues
    receptor_residues = [res for res in traj.topology.residues if res.chain.chain_id == receptor_chain]
    
    print(f"[DEBUG] Analyzing {len(receptor_residues)} receptor residues")
    print(f"[DEBUG] Input: {len(hbond_pairs)} hydrogen bonds, {len(charge_pairs)} charge pairs, {len(pi_contact_pairs)} pi contacts")
    
    # Debug: Print all input interactions
    print(f"[DEBUG] === INPUT INTERACTIONS ===")
    print(f"[DEBUG] Hydrogen bonds:")
    for i, hbond in enumerate(hbond_pairs):
        print(f"[DEBUG]   HB{i+1}: {hbond['donor_chain']}{hbond['donor_residue']}({hbond['donor_atom']}) -> {hbond['acceptor_chain']}{hbond['acceptor_residue']}({hbond['acceptor_atom']}) - {hbond['distance']:.2f}Å")
    
    print(f"[DEBUG] Charge pairs:")
    for i, charge in enumerate(charge_pairs):
        print(f"[DEBUG]   CP{i+1}: {charge['receptor_chain']}{charge['receptor_residue']}({charge['receptor_residue_name']}) <-> {charge['ligand_chain']}{charge['ligand_residue']}({charge['ligand_residue_name']}) - {charge['distance']:.2f}Å")
    
    # Create mappings for quick lookup - FIXED LOGIC
    hbond_residues = {}
    charge_residues = {}
    pi_residues = {}
    
    # Process hydrogen bonds - only count receptor residues as having hbonds
    print(f"[DEBUG] === PROCESSING HYDROGEN BONDS ===")
    for hbond in hbond_pairs:
        if hbond['donor_chain'] == receptor_chain:
            # Receptor residue is donor
            receptor_res_num = hbond['donor_residue']
            ligand_res_num = hbond['acceptor_residue']
            if receptor_res_num not in hbond_residues:
                hbond_residues[receptor_res_num] = []
            hbond_residues[receptor_res_num].append(hbond)
            print(f"[DEBUG]   Receptor {receptor_res_num} is DONOR -> Ligand {ligand_res_num}")
        elif hbond['acceptor_chain'] == receptor_chain:
            # Receptor residue is acceptor
            receptor_res_num = hbond['acceptor_residue']
            ligand_res_num = hbond['donor_residue']
            if receptor_res_num not in hbond_residues:
                hbond_residues[receptor_res_num] = []
            hbond_residues[receptor_res_num].append(hbond)
            print(f"[DEBUG]   Receptor {receptor_res_num} is ACCEPTOR <- Ligand {ligand_res_num}")
        else:
            print(f"[DEBUG]   WARNING: H-bond not involving receptor chain {receptor_chain}: {hbond}")
    
    # Process charge pairs - only count receptor residues
    print(f"[DEBUG] === PROCESSING CHARGE PAIRS ===")
    for charge_pair in charge_pairs:
        receptor_res_num = charge_pair['receptor_residue']
        ligand_res_num = charge_pair['ligand_residue']
        if receptor_res_num not in charge_residues:
            charge_residues[receptor_res_num] = []
        charge_residues[receptor_res_num].append(charge_pair)
        print(f"[DEBUG]   Receptor {receptor_res_num} <-> Ligand {ligand_res_num}")
    # raise Exception("charge_residues")
    
    print(f"[DEBUG] Receptor residues with hydrogen bonds: {len(hbond_residues)}")
    print(f"[DEBUG] Receptor residues with charge complementarity: {len(charge_residues)}")
    
    # Debug: Show which receptor residues have interactions
    print(f"[DEBUG] === PROCESSING PI CONTACTS ===")
    for pi_contact in pi_contact_pairs:
        receptor_res_num = pi_contact['receptor_residue']
        ligand_res_num = pi_contact['ligand_residue']
        if receptor_res_num not in pi_residues:
            pi_residues[receptor_res_num] = []
        pi_residues[receptor_res_num].append(pi_contact)
        print(f"[DEBUG]   Receptor {receptor_res_num} <-> Ligand {ligand_res_num} (pi)")
    
    print(f"[DEBUG] Receptor residues with pi contacts: {len(pi_residues)}")
    
    print(f"[DEBUG] === RECEPTOR RESIDUES WITH INTERACTIONS ===")
    interacting_residues = set(list(hbond_residues.keys()) + list(charge_residues.keys()) + list(pi_residues.keys()))
    for res_num in sorted(interacting_residues):
        hb_count = len(hbond_residues[res_num]) if res_num in hbond_residues else 0
        ch_count = len(charge_residues[res_num]) if res_num in charge_residues else 0
        pi_count = len(pi_residues[res_num]) if res_num in pi_residues else 0
        print(f"[DEBUG]   Residue {res_num}: {hb_count} H-bonds, {ch_count} charge pairs, {pi_count} pi contacts")
    
    # Debug: Show the actual dictionary contents
    print(f"[DEBUG] === HYDROGEN BOND RESIDUES DICTIONARY ===")
    for res_num, hbonds in hbond_residues.items():
        print(f"[DEBUG]   Residue {res_num}: {len(hbonds)} H-bonds")
        for hbond in hbonds:
            print(f"[DEBUG]     {hbond}")
    
    print(f"[DEBUG] === CHARGE RESIDUES DICTIONARY ===")
    for res_num, charges in charge_residues.items():
        print(f"[DEBUG]   Residue {res_num}: {len(charges)} charge pairs")
        for charge in charges:
            print(f"[DEBUG]     {charge}")
    
    print(f"[DEBUG] === PI RESIDUES DICTIONARY ===")
    for res_num, pi_contacts in pi_residues.items():
        print(f"[DEBUG]   Residue {res_num}: {len(pi_contacts)} pi contacts")
        for pi_contact in pi_contacts:
            print(f"[DEBUG]     {pi_contact}")
    
    # Debug: Check for residues that appear in both dictionaries
    print(f"[DEBUG] === CHECKING FOR OVERLAP ===")
    hbond_keys = set(hbond_residues.keys())
    charge_keys = set(charge_residues.keys())
    pi_keys = set(pi_residues.keys())
    overlap_all = hbond_keys.union(charge_keys).intersection(pi_keys)
    print(f"[DEBUG] Residues with H-bond & charge: {sorted(hbond_keys.intersection(charge_keys))}")
    print(f"[DEBUG] Residues with H-bond & pi: {sorted(hbond_keys.intersection(pi_keys))}")
    print(f"[DEBUG] Residues with charge & pi: {sorted(charge_keys.intersection(pi_keys))}")
    print(f"[DEBUG] Residues with all three: {sorted(overlap_all)}")
    print(f"[DEBUG] Residues only in H-bond dict: {sorted(hbond_keys - charge_keys - pi_keys)}")
    print(f"[DEBUG] Residues only in charge dict: {sorted(charge_keys - hbond_keys - pi_keys)}")
    print(f"[DEBUG] Residues only in pi dict: {sorted(pi_keys - hbond_keys - charge_keys)}")
    
    # Analyze each receptor residue
    print(f"[DEBUG] === ANALYZING EACH RECEPTOR RESIDUE ===")
    for residue in receptor_residues:
        res_num = residue.resSeq
        res_name = residue.name
        
        has_hbond = res_num in hbond_residues
        has_charge_complement = res_num in charge_residues
        has_pi_contact = res_num in pi_residues
        
        # Determine interaction category
        interaction_flags = {
            'hbond': has_hbond,
            'charge': has_charge_complement,
            'pi': has_pi_contact
        }
        flag_count = sum(interaction_flags.values())
        
        if flag_count == 0:
            interaction_category = 'none'
        elif flag_count == 1:
            if has_hbond:
                interaction_category = 'hbond_only'
            elif has_charge_complement:
                interaction_category = 'charge_only'
            else:
                interaction_category = 'pi_only'
        elif flag_count == 2:
            if not has_pi_contact:
                interaction_category = 'hbond_and_charge'
            elif not has_charge_complement:
                interaction_category = 'hbond_and_pi'
            else:
                interaction_category = 'charge_and_pi'
        else:
            interaction_category = 'hbond_charge_pi'
        
        # Get partner residues - FIXED LOGIC with proper deduplication
        partner_residues = set()
        hbond_partner_residues = set()
        charge_partner_residues = set()
        pi_partner_residues = set()
        
        if has_hbond:
            print(f"[DEBUG]   Residue {res_num} ({res_name}) has {len(hbond_residues.get(res_num, []))} H-bonds:")
            for hbond in hbond_residues.get(res_num, []):
                if hbond['donor_residue'] == res_num:
                    partner_res = hbond['acceptor_residue']
                    partner_chain = hbond['acceptor_chain']
                    print(f"[DEBUG]     As DONOR: -> {partner_chain}{partner_res}")
                    partner_residues.add(partner_res)
                    hbond_partner_residues.add(partner_res)
                else:
                    partner_res = hbond['donor_residue']
                    partner_chain = hbond['donor_chain']
                    print(f"[DEBUG]     As ACCEPTOR: <- {partner_chain}{partner_res}")
                    partner_residues.add(partner_res)
                    hbond_partner_residues.add(partner_res)
        
        if has_charge_complement:
            print(f"[DEBUG]   Residue {res_num} ({res_name}) has {len(charge_residues.get(res_num, []))} charge pairs:")
            for charge_pair in charge_residues.get(res_num, []):
                partner_res = charge_pair['ligand_residue']
                partner_chain = charge_pair['ligand_chain']
                print(f"[DEBUG]     Charge pair: <-> {partner_chain}{partner_res}")
                partner_residues.add(partner_res)
                charge_partner_residues.add(partner_res)
    
        if has_pi_contact:
            print(f"[DEBUG]   Residue {res_num} ({res_name}) has {len(pi_residues.get(res_num, []))} pi contacts:")
            for pi_contact in pi_residues.get(res_num, []):
                partner_res = pi_contact['ligand_residue']
                partner_chain = pi_contact['ligand_chain']
                print(f"[DEBUG]     Pi contact: <-> {partner_chain}{partner_res}")
                partner_residues.add(partner_res)
                pi_partner_residues.add(partner_res)
        
        partner_residues = sorted(partner_residues)  # Convert set to sorted list
        
        # Debug: Show final partner residues for interacting residues
        if has_hbond or has_charge_complement or has_pi_contact:
            print(f"[DEBUG]   Residue {res_num} ({res_name}) final partners: {sorted(partner_residues)}")
        
        # Debug: Show the actual values being set
        print(f"[DEBUG]   Residue {res_num} ({res_name}): has_hbond={has_hbond}, has_charge_complement={has_charge_complement}, has_pi_contact={has_pi_contact}, category={interaction_category}")
        
        residue_info.append({
            'residue_number': res_num,
            'residue_name': res_name,
            'chain_id': receptor_chain,
            'has_hbond': has_hbond,
            'has_charge_complement': has_charge_complement,
            'has_pi_contact': has_pi_contact,
            'interaction_category': interaction_category,
            'partner_residues': partner_residues,
            'hbond_count': len(hbond_residues.get(res_num, [])),
            'charge_pair_count': len(charge_residues.get(res_num, [])),
            'pi_contact_count': len(pi_residues.get(res_num, [])),
            'hbond_partner_residues': sorted(hbond_partner_residues),
            'charge_partner_residues': sorted(charge_partner_residues),
            'pi_partner_residues': sorted(pi_partner_residues),
            'is_charged': res_name in POSITIVE_RESIDUES + NEGATIVE_RESIDUES,
            'charge_type': 'positive' if res_name in POSITIVE_RESIDUES else \
                          'negative' if res_name in NEGATIVE_RESIDUES else 'neutral'
        })
    
    print(f"[DEBUG] === FINAL SUMMARY ===")
    total_interacting = sum(1 for info in residue_info if info['has_hbond'] or info['has_charge_complement'] or info['has_pi_contact'])
    print(f"[DEBUG] Total interacting residues: {total_interacting}")
    # raise Exception("residue_info")
    return residue_info


def write_interaction_analysis_csv(residue_info: List[Dict], output_path: str) -> None:
    """
    Write interaction analysis results to CSV file.
    
    Args:
        residue_info: List of residue information dictionaries from compute_interaction_filter
        output_path: Path to output CSV file
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'residue_number', 'residue_name', 'chain_id',
            'has_hbond', 'has_charge_complement', 'has_pi_contact',
            'interaction_category', 'partner_residues',
            'hbond_count', 'charge_pair_count', 'pi_contact_count',
            'hbond_partner_residues', 'charge_partner_residues', 'pi_partner_residues',
            'is_charged', 'charge_type'
        ])
        
        for info in residue_info:
            writer.writerow([
                info['residue_number'],
                info['residue_name'],
                info['chain_id'],
                info['has_hbond'],
                info['has_charge_complement'],
                info['has_pi_contact'],
                info['interaction_category'],
                ','.join(map(str, info['partner_residues'])),
                info['hbond_count'],
                info['charge_pair_count'],
                info['pi_contact_count'],
                ','.join(map(str, info['hbond_partner_residues'])),
                ','.join(map(str, info['charge_partner_residues'])),
                ','.join(map(str, info['pi_partner_residues'])),
                info['is_charged'],
                info['charge_type']
            ])


def get_interaction_summary(interaction_results: Dict) -> str:
    """
    Generate a summary string of interaction analysis results.
    
    Args:
        interaction_results: Results dictionary from compute_interaction_filter
        
    Returns:
        Formatted summary string
    """
    if 'error' in interaction_results:
        return f"Interaction analysis failed: {interaction_results['error']}"
    
    summary = f"""
Interaction Analysis Summary:
- Receptor chain: {interaction_results['receptor_chain']}
- Ligand chain: {interaction_results['ligand_chain']}
- Total residues analyzed: {len(interaction_results['residue_info'])}
- Residues with hydrogen bonds: {interaction_results['n_hbond_residues']}
- Residues with charge complementarity: {interaction_results['n_charge_residues']}
- Residues with pi contacts: {interaction_results['n_pi_residues']}
- Residues with multiple interaction types: {interaction_results['n_multiple_residues']}
- Total interacting residues: {interaction_results['n_interacting_residues']}
- Fraction interacting: {interaction_results['fraction_interacting']:.2%}
- Distance threshold (charge): {interaction_results['distance_threshold']} Å
- Distance threshold (pi): {interaction_results['pi_distance_threshold']} Å
- Total hydrogen bonds: {len(interaction_results['hbond_pairs'])}
- Total charge pairs: {len(interaction_results['charge_pairs'])}
- Total pi contacts: {len(interaction_results['pi_contact_pairs'])}
"""
    return summary.strip()


def compute_ca_distance_filter(holo_pdb_path: str, distance_threshold: float = 6.0,
                              receptor_chain_id: str = None, ligand_chain_id: str = None) -> Dict:
    """
    Compute CA-CA distance filter for a receptor-ligand complex.
    
    Identifies receptor residues where the CA atom is within the specified distance
    threshold of any ligand CA atom.
    
    Args:
        holo_pdb_path: Path to the holo PDB file (receptor-ligand complex)
        distance_threshold: Distance threshold for CA-CA proximity (Å)
        receptor_chain_id: Optional receptor chain ID (if None, auto-detect longest chain)
        ligand_chain_id: Optional ligand chain ID (if None, auto-detect shortest chain)
        
    Returns:
        Dictionary containing:
        - residue_info: List of dicts with detailed residue distance information
        - receptor_chain: ID of receptor chain
        - ligand_chain: ID of ligand chain
        - n_passing: Count of residues passing the distance filter
        - fraction_passing: Fraction of total residues passing the filter
        - distance_threshold: Distance threshold used
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
        
        print(f"[DEBUG] === CA DISTANCE FILTER ===")
        print(f"[DEBUG] Identified receptor chain: {receptor_chain}, ligand chain: {ligand_chain}")
        
        # Print sequence lengths for confirmation
        receptor_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == receptor_chain)
        ligand_chain_obj = next(c for c in holo_traj.topology.chains if c.chain_id == ligand_chain)
        print(f"[DEBUG] Receptor chain length: {receptor_chain_obj.n_residues} residues")
        print(f"[DEBUG] Ligand chain length: {ligand_chain_obj.n_residues} residues")
        
        if receptor_chain == ligand_chain:
            raise ValueError("Only one chain found in structure - cannot perform CA distance analysis")
        
        # Get all chain IDs
        all_chains = set([atom.residue.chain.chain_id for atom in holo_traj.topology.atoms])
        print(f"[DEBUG] All chains in structure: {all_chains}")
        if receptor_chain not in all_chains or ligand_chain not in all_chains:
            raise ValueError(f"Chain identification failed. Found chains: {all_chains}")
        
        # Get receptor and ligand residues
        receptor_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == receptor_chain]
        ligand_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == ligand_chain]
        
        print(f"[DEBUG] Analyzing {len(receptor_residues)} receptor residues against {len(ligand_residues)} ligand residues")
        
        # Analyze each receptor residue
        residue_info = []
        n_passing = 0
        
        for receptor_res in receptor_residues:
            # Skip non-standard residues
            if receptor_res.name not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                continue
            
            # Find CA atom for this receptor residue
            receptor_ca = None
            for atom in receptor_res.atoms:
                if atom.name == 'CA':
                    receptor_ca = atom
                    break
            
            if receptor_ca is None:
                print(f"[DEBUG] Warning: No CA atom found for receptor residue {receptor_res.name}{receptor_res.resSeq}")
                continue
            
            # Calculate minimum distance to any ligand CA atom
            min_distance = float('inf')
            nearest_ligand_residue = None
            
            for ligand_res in ligand_residues:
                # Skip non-standard residues
                if ligand_res.name not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                          'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    continue
                
                # Find CA atom for this ligand residue
                ligand_ca = None
                for atom in ligand_res.atoms:
                    if atom.name == 'CA':
                        ligand_ca = atom
                        break
                
                if ligand_ca is None:
                    continue
                
                # Calculate distance (MDTraj uses nm internally, convert to Å)
                distance = np.linalg.norm(holo_traj.xyz[0, receptor_ca.index] - holo_traj.xyz[0, ligand_ca.index]) * 10.0
                
                if distance < min_distance:
                    min_distance = distance
                    nearest_ligand_residue = ligand_res.resSeq
            
            # Determine if residue passes the filter
            passes_filter = min_distance < distance_threshold
            if passes_filter:
                n_passing += 1
            
            residue_info.append({
                'residue_number': receptor_res.resSeq,
                'residue_name': receptor_res.name,
                'chain_id': receptor_chain,
                'min_ca_distance': min_distance,
                'passes_filter': passes_filter,
                'nearest_ligand_residue': nearest_ligand_residue,
                'distance_threshold': distance_threshold
            })
        
        fraction_passing = n_passing / len(residue_info) if residue_info else 0.0
        
        print(f"[DEBUG] CA distance filter results:")
        print(f"[DEBUG] - Total receptor residues analyzed: {len(residue_info)}")
        print(f"[DEBUG] - Residues passing filter: {n_passing}")
        print(f"[DEBUG] - Fraction passing: {fraction_passing:.2%}")
        print(f"[DEBUG] - Distance threshold: {distance_threshold} Å")
        
        return {
            'residue_info': residue_info,
            'receptor_chain': receptor_chain,
            'ligand_chain': ligand_chain,
            'n_passing': n_passing,
            'fraction_passing': fraction_passing,
            'distance_threshold': distance_threshold
        }
        
    except Exception as e:
        print(f"Error in compute_ca_distance_filter: {e}")
        return {
            'error': str(e),
            'residue_info': [],
            'receptor_chain': '',
            'ligand_chain': '',
            'n_passing': 0,
            'fraction_passing': 0.0,
            'distance_threshold': distance_threshold
        }


def write_ca_distance_csv(residue_info: List[Dict], output_path: str) -> None:
    """
    Write CA distance filter results to CSV file.
    
    Args:
        residue_info: List of residue information dictionaries from compute_ca_distance_filter
        output_path: Path to output CSV file
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'residue_number', 'residue_name', 'chain_id', 'min_ca_distance', 
            'passes_filter', 'nearest_ligand_residue', 'distance_threshold'
        ])
        
        for info in residue_info:
            writer.writerow([
                info['residue_number'],
                info['residue_name'],
                info['chain_id'],
                f"{info['min_ca_distance']:.4f}",
                info['passes_filter'],
                info['nearest_ligand_residue'],
                f"{info['distance_threshold']:.4f}"
            ])


def get_ca_distance_summary(ca_distance_results: Dict) -> str:
    """
    Generate a summary string of CA distance filter results.
    
    Args:
        ca_distance_results: Results dictionary from compute_ca_distance_filter
        
    Returns:
        Formatted summary string
    """
    if 'error' in ca_distance_results:
        return f"CA distance filter failed: {ca_distance_results['error']}"
    
    summary = f"""
CA Distance Filter Summary:
- Receptor chain: {ca_distance_results['receptor_chain']}
- Ligand chain: {ca_distance_results['ligand_chain']}
- Total residues analyzed: {len(ca_distance_results['residue_info'])}
- Residues passing filter: {ca_distance_results['n_passing']}
- Fraction passing: {ca_distance_results['fraction_passing']:.2%}
- Distance threshold: {ca_distance_results['distance_threshold']} Å
"""
    return summary.strip()


STANDARD_RESIDUES = frozenset([
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
])


def compute_nn_distance_filter(holo_pdb_path: str, distance_threshold: float = 6.0,
                               receptor_chain_id: str = None, ligand_chain_id: str = None) -> Dict:
    """
    Compute N-N (backbone nitrogen) distance filter for a receptor-ligand complex.
    
    Identifies receptor residues where the backbone N atom is within the specified
    distance threshold of any ligand backbone N atom.
    
    Args:
        holo_pdb_path: Path to the holo PDB file (receptor-ligand complex)
        distance_threshold: Distance threshold for N-N proximity (Å)
        receptor_chain_id: Optional receptor chain ID (if None, auto-detect longest chain)
        ligand_chain_id: Optional ligand chain ID (if None, auto-detect shortest chain)
        
    Returns:
        Dictionary containing:
        - residue_info: List of dicts with detailed residue distance information
        - receptor_chain: ID of receptor chain
        - ligand_chain: ID of ligand chain
        - n_passing: Count of residues passing the distance filter
        - fraction_passing: Fraction of total residues passing the filter
        - distance_threshold: Distance threshold used
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
            raise ValueError("Only one chain found in structure - cannot perform N-N distance analysis")
        
        receptor_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == receptor_chain]
        ligand_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == ligand_chain]
        
        residue_info = []
        n_passing = 0
        
        for receptor_res in receptor_residues:
            if receptor_res.name not in STANDARD_RESIDUES:
                continue
            
            receptor_n = None
            for atom in receptor_res.atoms:
                if atom.name == 'N':
                    receptor_n = atom
                    break
            
            if receptor_n is None:
                continue
            
            min_distance = float('inf')
            nearest_ligand_residue = None
            
            for ligand_res in ligand_residues:
                if ligand_res.name not in STANDARD_RESIDUES:
                    continue
                
                ligand_n = None
                for atom in ligand_res.atoms:
                    if atom.name == 'N':
                        ligand_n = atom
                        break
                
                if ligand_n is None:
                    continue
                
                distance = np.linalg.norm(holo_traj.xyz[0, receptor_n.index] - holo_traj.xyz[0, ligand_n.index]) * 10.0
                if distance < min_distance:
                    min_distance = distance
                    nearest_ligand_residue = ligand_res.resSeq
            
            passes_filter = min_distance < distance_threshold
            if passes_filter:
                n_passing += 1
            
            residue_info.append({
                'residue_number': receptor_res.resSeq,
                'residue_name': receptor_res.name,
                'chain_id': receptor_chain,
                'min_nn_distance': min_distance,
                'passes_filter': passes_filter,
                'nearest_ligand_residue': nearest_ligand_residue,
                'distance_threshold': distance_threshold
            })
        
        fraction_passing = n_passing / len(residue_info) if residue_info else 0.0
        
        return {
            'residue_info': residue_info,
            'receptor_chain': receptor_chain,
            'ligand_chain': ligand_chain,
            'n_passing': n_passing,
            'fraction_passing': fraction_passing,
            'distance_threshold': distance_threshold
        }
        
    except Exception as e:
        print(f"Error in compute_nn_distance_filter: {e}")
        return {
            'error': str(e),
            'residue_info': [],
            'receptor_chain': '',
            'ligand_chain': '',
            'n_passing': 0,
            'fraction_passing': 0.0,
            'distance_threshold': distance_threshold
        }


def compute_min_atom_distance_filter(holo_pdb_path: str, distance_threshold: float = 6.0,
                                    receptor_chain_id: str = None, ligand_chain_id: str = None) -> Dict:
    """
    Compute minimum inter-chain atomic distance filter (any atom type).
    
    For each receptor residue, computes the minimum distance between any atom in that
    residue and any atom in any ligand residue.
    
    Args:
        holo_pdb_path: Path to the holo PDB file (receptor-ligand complex)
        distance_threshold: Distance threshold for proximity (Å)
        receptor_chain_id: Optional receptor chain ID (if None, auto-detect longest chain)
        ligand_chain_id: Optional ligand chain ID (if None, auto-detect shortest chain)
        
    Returns:
        Dictionary containing:
        - residue_info: List of dicts with detailed residue distance information
        - receptor_chain: ID of receptor chain
        - ligand_chain: ID of ligand chain
        - n_passing: Count of residues passing the distance filter
        - fraction_passing: Fraction of total residues passing the filter
        - distance_threshold: Distance threshold used
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
            raise ValueError("Only one chain found in structure - cannot perform any-atom distance analysis")
        
        receptor_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == receptor_chain]
        ligand_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == ligand_chain]
        
        residue_info = []
        n_passing = 0
        
        for receptor_res in receptor_residues:
            if receptor_res.name not in STANDARD_RESIDUES:
                continue
            
            min_distance = float('inf')
            nearest_ligand_residue = None
            
            for receptor_atom in receptor_res.atoms:
                for ligand_res in ligand_residues:
                    if ligand_res.name not in STANDARD_RESIDUES:
                        continue
                    for ligand_atom in ligand_res.atoms:
                        distance = np.linalg.norm(
                            holo_traj.xyz[0, receptor_atom.index] - holo_traj.xyz[0, ligand_atom.index]
                        ) * 10.0
                        if distance < min_distance:
                            min_distance = distance
                            nearest_ligand_residue = ligand_res.resSeq
            
            passes_filter = min_distance < distance_threshold
            passes_sub_2A_filter = min_distance < _get_direct_contact_threshold()
            if passes_filter:
                n_passing += 1
            
            residue_info.append({
                'residue_number': receptor_res.resSeq,
                'residue_name': receptor_res.name,
                'chain_id': receptor_chain,
                'min_any_atom_distance': min_distance,
                'passes_filter': passes_filter,
                'passes_sub_2A_filter': passes_sub_2A_filter,
                'nearest_ligand_residue': nearest_ligand_residue,
                'distance_threshold': distance_threshold
            })
        
        fraction_passing = n_passing / len(residue_info) if residue_info else 0.0
        
        return {
            'residue_info': residue_info,
            'receptor_chain': receptor_chain,
            'ligand_chain': ligand_chain,
            'n_passing': n_passing,
            'fraction_passing': fraction_passing,
            'distance_threshold': distance_threshold
        }
        
    except Exception as e:
        print(f"Error in compute_min_atom_distance_filter: {e}")
        return {
            'error': str(e),
            'residue_info': [],
            'receptor_chain': '',
            'ligand_chain': '',
            'n_passing': 0,
            'fraction_passing': 0.0,
            'distance_threshold': distance_threshold
        }


def write_nn_distance_csv(residue_info: List[Dict], output_path: str) -> None:
    """Write N-N distance filter results to CSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'residue_number', 'residue_name', 'chain_id', 'min_nn_distance',
            'passes_filter', 'nearest_ligand_residue', 'distance_threshold'
        ])
        for info in residue_info:
            writer.writerow([
                info['residue_number'],
                info['residue_name'],
                info['chain_id'],
                f"{info['min_nn_distance']:.4f}",
                info['passes_filter'],
                info['nearest_ligand_residue'],
                f"{info['distance_threshold']:.4f}"
            ])


def write_any_atom_distance_csv(residue_info: List[Dict], output_path: str) -> None:
    """Write any-atom distance filter results to CSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'residue_number', 'residue_name', 'chain_id', 'min_any_atom_distance',
            'passes_filter', 'passes_sub_2A_filter', 'nearest_ligand_residue', 'distance_threshold'
        ])
        for info in residue_info:
            writer.writerow([
                info['residue_number'],
                info['residue_name'],
                info['chain_id'],
                f"{info['min_any_atom_distance']:.4f}",
                info['passes_filter'],
                info.get('passes_sub_2A_filter', False),
                info['nearest_ligand_residue'],
                f"{info['distance_threshold']:.4f}"
            ])


def get_nn_distance_summary(nn_distance_results: Dict) -> str:
    """Generate a summary string of N-N distance filter results."""
    if 'error' in nn_distance_results:
        return f"N-N distance filter failed: {nn_distance_results['error']}"
    summary = f"""
N-N Distance Filter Summary:
- Receptor chain: {nn_distance_results['receptor_chain']}
- Ligand chain: {nn_distance_results['ligand_chain']}
- Total residues analyzed: {len(nn_distance_results['residue_info'])}
- Residues passing filter: {nn_distance_results['n_passing']}
- Fraction passing: {nn_distance_results['fraction_passing']:.2%}
- Distance threshold: {nn_distance_results['distance_threshold']} Å
"""
    return summary.strip()


def get_any_atom_distance_summary(any_atom_results: Dict) -> str:
    """Generate a summary string of any-atom distance filter results."""
    if 'error' in any_atom_results:
        return f"Any-atom distance filter failed: {any_atom_results['error']}"
    summary = f"""
Any-Atom Distance Filter Summary:
- Receptor chain: {any_atom_results['receptor_chain']}
- Ligand chain: {any_atom_results['ligand_chain']}
- Total residues analyzed: {len(any_atom_results['residue_info'])}
- Residues passing filter: {any_atom_results['n_passing']}
- Fraction passing: {any_atom_results['fraction_passing']:.2%}
- Distance threshold: {any_atom_results['distance_threshold']} Å
"""
    return summary.strip()
