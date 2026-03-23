"""
Identify Peptide Binding Site Residues

This script analyzes a PDB file containing a protein-peptide complex and identifies
which peptide residues are in the binding site according to the criteria established
in pipeline.py and confusion_matrix_analysis.py:169-170.

Criteria: Peptide residues that interact with protein residues meeting:
is_occluded OR passes_ca_filter OR is_interacting

Usage:
    python scripts/identify_peptide_binding_site.py <pdb_file> [--protein-chain CHAIN] [--peptide-chain CHAIN] [--output OUTPUT.csv]
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Set, Dict, List, Optional

# Support running as a script (python scripts/identify_peptide_binding_site.py) or module (python -m scripts.identify_peptide_binding_site)
try:
    from .config import sasa_analysis, ca_distance_analysis
    from .sasa_analysis import compute_sasa_occlusion, find_medoid_model_from_pdb
    from .interaction_analysis import compute_interaction_filter, compute_ca_distance_filter, compute_min_atom_distance_filter
    from .util import get_longest_chain, get_shortest_chain
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import sasa_analysis, ca_distance_analysis
    from scripts.sasa_analysis import compute_sasa_occlusion, find_medoid_model_from_pdb
    from scripts.interaction_analysis import compute_interaction_filter, compute_ca_distance_filter, compute_min_atom_distance_filter
    from scripts.util import get_longest_chain, get_shortest_chain


def identify_peptide_binding_site_residues(
    pdb_path: str,
    receptor_chain_id: Optional[str] = None,
    ligand_chain_id: Optional[str] = None,
) -> Dict:
    """
    Identify ligand residues in the binding site.
    
    Args:
        pdb_path: Path to PDB file
        receptor_chain_id: Optional receptor chain ID (auto-detect if None)
        ligand_chain_id: Optional ligand chain ID (auto-detect if None)
        
    Returns:
        Dictionary containing:
        - peptide_binding_site_residues: Set of ligand residue numbers
        - receptor_chain: Receptor chain ID
        - ligand_chain: Ligand chain ID
        - protein_binding_site_residues: Set of receptor residue numbers meeting criteria
        - peptide_residue_info: List of dicts with ligand residue details
    """
    # Find medoid model if it's a multi-model PDB
    medoid_pdb_path = find_medoid_model_from_pdb(pdb_path)
    
    # Auto-detect chains if not provided
    if receptor_chain_id is None or ligand_chain_id is None:
        detected_receptor_chain = get_longest_chain(medoid_pdb_path)
        detected_ligand_chain = get_shortest_chain(medoid_pdb_path)
        
        if receptor_chain_id is None:
            receptor_chain_id = detected_receptor_chain
        if ligand_chain_id is None:
            ligand_chain_id = detected_ligand_chain
    
    print(f"[INFO] Using receptor chain: {receptor_chain_id}, ligand chain: {ligand_chain_id}")
    
    # Run SASA occlusion analysis
    print("[INFO] Running SASA occlusion analysis...")
    sasa_results = compute_sasa_occlusion(
        medoid_pdb_path,
        sasa_threshold=sasa_analysis.sasa_threshold,
        probe_radius=0.0,
        nh_mode='sum',
        mode='residue',
        receptor_chain_id=receptor_chain_id,
        ligand_chain_id=ligand_chain_id
    )
    
    # Run CA distance filter analysis
    print("[INFO] Running CA distance filter analysis...")
    ca_distance_results = compute_ca_distance_filter(
        medoid_pdb_path,
        distance_threshold=ca_distance_analysis.ca_distance_threshold,
        receptor_chain_id=receptor_chain_id,
        ligand_chain_id=ligand_chain_id
    )
    
    # Run interaction analysis
    print("[INFO] Running interaction analysis...")
    interaction_results = compute_interaction_filter(
        medoid_pdb_path,
        distance_threshold=4.5,
        pi_distance_threshold=6.0,
        receptor_chain_id=receptor_chain_id,
        ligand_chain_id=ligand_chain_id
    )
    
    # Run any-atom distance filter analysis (for min_any_atom_distance < 2A criterion)
    print("[INFO] Running any-atom distance filter analysis...")
    any_atom_results = compute_min_atom_distance_filter(
        medoid_pdb_path,
        distance_threshold=ca_distance_analysis.ca_distance_threshold,
        receptor_chain_id=receptor_chain_id,
        ligand_chain_id=ligand_chain_id
    )
    
    # Create lookup maps for protein residues
    sasa_map = {}
    for info in sasa_results.get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is not None:
            sasa_map[res_num] = {
                'is_occluded': info.get('is_occluded', False)
            }
    
    ca_map = {}
    for info in ca_distance_results.get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is not None:
            ca_map[res_num] = {
                'passes_filter': info.get('passes_filter', False),
                'nearest_ligand_residue': info.get('nearest_ligand_residue')
            }
    
    interaction_map = {}
    for info in interaction_results.get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is not None:
            has_hbond = info.get('has_hbond', False)
            has_charge_complement = info.get('has_charge_complement', False)
            has_pi_contact = info.get('has_pi_contact', False)
            is_interacting = has_hbond or has_charge_complement or has_pi_contact
            
            interaction_map[res_num] = {
                'is_interacting': is_interacting,
                'has_hbond': has_hbond,
                'has_charge_complement': has_charge_complement,
                'has_pi_contact': has_pi_contact,
                'partner_residues': info.get('partner_residues', [])
            }
    
    any_atom_map = {}
    for info in (any_atom_results or {}).get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is not None:
            any_atom_map[res_num] = {
                'passes_sub_2A_filter': info.get('passes_sub_2A_filter', False),
                'nearest_ligand_residue': info.get('nearest_ligand_residue')
            }
    
    # Identify protein residues meeting criteria: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance < 2A
    protein_binding_site_residues: Set[int] = set()
    
    all_protein_residues = set(sasa_map.keys()) | set(ca_map.keys()) | set(interaction_map.keys()) | set(any_atom_map.keys())
    
    for res_num in all_protein_residues:
        is_occluded = sasa_map.get(res_num, {}).get('is_occluded', False)
        passes_ca_filter = ca_map.get(res_num, {}).get('passes_filter', False)
        is_interacting = interaction_map.get(res_num, {}).get('is_interacting', False)
        passes_sub_2A_filter = any_atom_map.get(res_num, {}).get('passes_sub_2A_filter', False)
        
        if is_occluded or passes_ca_filter or is_interacting or passes_sub_2A_filter:
            protein_binding_site_residues.add(res_num)
    
    # Collect peptide residues from binding site protein residues
    peptide_binding_site_residues: Set[int] = set()
    peptide_residue_info_map: Dict[int, Dict] = {}
    
    for protein_res_num in protein_binding_site_residues:
        # From interaction results: get partner_residues
        interaction_data = interaction_map.get(protein_res_num, {})
        partner_residues = interaction_data.get('partner_residues', [])
        
        for peptide_res_num in partner_residues:
            peptide_binding_site_residues.add(peptide_res_num)
            
            # Track interaction types for this peptide residue
            if peptide_res_num not in peptide_residue_info_map:
                peptide_residue_info_map[peptide_res_num] = {
                    'interaction_types': set(),
                    'protein_partners': []
                }
            
            # Add interaction types
            if interaction_data.get('has_hbond', False):
                peptide_residue_info_map[peptide_res_num]['interaction_types'].add('hbond')
            if interaction_data.get('has_charge_complement', False):
                peptide_residue_info_map[peptide_res_num]['interaction_types'].add('charge')
            if interaction_data.get('has_pi_contact', False):
                peptide_residue_info_map[peptide_res_num]['interaction_types'].add('pi')
            
            peptide_residue_info_map[peptide_res_num]['protein_partners'].append(protein_res_num)
        
        # From CA distance results: get nearest_ligand_residue
        ca_data = ca_map.get(protein_res_num, {})
        if ca_data.get('passes_filter', False):
            nearest_ligand_residue = ca_data.get('nearest_ligand_residue')
            if nearest_ligand_residue is not None:
                peptide_binding_site_residues.add(nearest_ligand_residue)
                
                if nearest_ligand_residue not in peptide_residue_info_map:
                    peptide_residue_info_map[nearest_ligand_residue] = {
                        'interaction_types': set(),
                        'protein_partners': []
                    }
                
                peptide_residue_info_map[nearest_ligand_residue]['interaction_types'].add('ca_distance')
                peptide_residue_info_map[nearest_ligand_residue]['protein_partners'].append(protein_res_num)
        
        # From any-atom distance results: get nearest_ligand_residue when passes_sub_2A_filter
        any_atom_data = any_atom_map.get(protein_res_num, {})
        if any_atom_data.get('passes_sub_2A_filter', False):
            nearest_ligand_residue = any_atom_data.get('nearest_ligand_residue')
            if nearest_ligand_residue is not None:
                peptide_binding_site_residues.add(nearest_ligand_residue)
                
                if nearest_ligand_residue not in peptide_residue_info_map:
                    peptide_residue_info_map[nearest_ligand_residue] = {
                        'interaction_types': set(),
                        'protein_partners': []
                    }
                
                peptide_residue_info_map[nearest_ligand_residue]['interaction_types'].add('any_atom_sub_2A')
                peptide_residue_info_map[nearest_ligand_residue]['protein_partners'].append(protein_res_num)
    
    # Build peptide residue info list (get residue names from peptide chain if available)
    import mdtraj as md
    holo_traj = md.load_pdb(medoid_pdb_path)
    ligand_residues = [res for res in holo_traj.topology.residues if res.chain.chain_id == ligand_chain_id]
    ligand_residue_name_map = {res.resSeq: res.name for res in ligand_residues}
    
    peptide_residue_info = []
    for peptide_res_num in sorted(peptide_binding_site_residues):
        info = peptide_residue_info_map.get(peptide_res_num, {})
        interaction_types = sorted(info.get('interaction_types', set()))
        protein_partners = sorted(info.get('protein_partners', []))
        residue_name = ligand_residue_name_map.get(peptide_res_num, 'UNK')
        
        peptide_residue_info.append({
            'peptide_residue_number': peptide_res_num,
            'residue_name': residue_name,
            'chain_id': ligand_chain_id,
            'interaction_types': interaction_types,
            'interaction_type_str': '+'.join(interaction_types) if interaction_types else 'none',
            'protein_partners': protein_partners
        })
    
    return {
        'peptide_binding_site_residues': peptide_binding_site_residues,
        'receptor_chain': receptor_chain_id,
        'ligand_chain': ligand_chain_id,
        'protein_binding_site_residues': protein_binding_site_residues,
        'peptide_residue_info': peptide_residue_info,
        'pdb_path': pdb_path
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Identify peptide residues in the binding site according to pipeline criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Criteria: Peptide residues that interact with protein residues meeting:
is_occluded OR passes_ca_filter OR is_interacting

Example:
    python scripts/identify_peptide_binding_site.py PDB_FILES/1cf4.pdb
    python scripts/identify_peptide_binding_site.py PDB_FILES/1cf4.pdb --output results.csv
        """
    )
    parser.add_argument('pdb_file', help='Path to PDB file')
    parser.add_argument('--receptor-chain', help='Receptor chain ID (auto-detect if not provided)')
    parser.add_argument('--ligand-chain', help='Ligand chain ID (auto-detect if not provided)')
    parser.add_argument('--output', help='Optional path to save results as CSV')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found: {args.pdb_file}")
        return
    
    # Run analysis
    try:
        results = identify_peptide_binding_site_residues(
            args.pdb_file,
            receptor_chain_id=getattr(args, 'receptor_chain', None),
            ligand_chain_id=getattr(args, 'ligand_chain', None)
        )
        
        # Print results
        print("\n" + "="*60)
        print("Peptide Binding Site Analysis")
        print("="*60)
        print(f"PDB file: {os.path.basename(results['pdb_path'])}")
        print(f"Receptor chain: {results['receptor_chain']}")
        print(f"Ligand chain: {results['ligand_chain']}")
        print()
        print(f"Protein residues in binding site: {len(results['protein_binding_site_residues'])}")
        print(f"Protein binding site residues: {sorted(results['protein_binding_site_residues'])}")
        print()
        print(f"Peptide residues in binding site: {len(results['peptide_binding_site_residues'])}")
        print(f"Peptide binding site residues: {sorted(results['peptide_binding_site_residues'])}")
        print("="*60)
        
        # Save to CSV if requested
        if args.output:
            with open(args.output, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'peptide_residue_number',
                    'residue_name',
                    'chain_id',
                    'interaction_types',
                    'protein_partners'
                ])
                for info in results['peptide_residue_info']:
                    writer.writerow([
                        info['peptide_residue_number'],
                        info['residue_name'],
                        info['chain_id'],
                        info['interaction_type_str'],
                        ','.join(str(p) for p in info['protein_partners'])
                    ])
            print(f"\nResults saved to: {args.output}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return


if __name__ == "__main__":
    main()
