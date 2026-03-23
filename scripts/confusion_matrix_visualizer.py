#!/usr/bin/env python3
"""
Confusion Matrix PyMOL Visualizer

This script creates PyMOL visualizations based on confusion matrix analysis results.
It takes a holo PDB file, positive strategy, and significance threshold as input,
then automatically finds the first model (for multi-model PDBs) and generates a 
temporary .pse file with classification coloring and automatically loads it in PyMOL.

Usage:
    python scripts/confusion_matrix_visualizer.py -p 1cf4 -s occ -t sig
    python scripts/confusion_matrix_visualizer.py -p 1d5g -s int -t 1sd
    python scripts/confusion_matrix_visualizer.py -p 1cf4 -s all -t 2sd -o my_viz.pse
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
import subprocess
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

try:
    from .confusion_matrix_analysis import ConfusionMatrix, compute_confusion_matrix, read_master_alignment
    from .config import classification_colors, hex_to_rgb01
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.confusion_matrix_analysis import ConfusionMatrix, compute_confusion_matrix, read_master_alignment
    from scripts.config import classification_colors, hex_to_rgb01


def find_medoid_model_simple(pdb_path: str) -> str:
    """
    Simple medoid model selection that extracts the first model from multi-model PDB files.
    This is a lightweight alternative that doesn't require mdtraj dependency.
    
    Args:
        pdb_path: Path to the PDB file
        
    Returns:
        Path to the medoid model PDB file (or original if single model)
    """
    try:
        # Check if it's a multi-model PDB
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        
        has_models = any(line.startswith('MODEL') for line in lines)
        
        if not has_models:
            print("[VISUALIZER] Single model PDB file, using as-is")
            return pdb_path
        
        print("[VISUALIZER] Multi-model PDB detected, extracting first model...")
        
        # Extract first model
        model_lines = []
        in_first_model = False
        
        for line in lines:
            if line.startswith('MODEL'):
                if 'MODEL        1' in line or 'MODEL       1' in line:
                    in_first_model = True
                    continue
                elif in_first_model:
                    break
            elif line.startswith('ENDMDL'):
                if in_first_model:
                    break
            elif in_first_model:
                model_lines.append(line)
        
        if not model_lines:
            print("[VISUALIZER] Could not extract first model, using original file")
            return pdb_path
        
        # Create temporary file for the first model
        import tempfile
        temp_fd, temp_path = tempfile.mkstemp(suffix='.pdb', prefix='medoid_')
        
        with os.fdopen(temp_fd, 'w') as f:
            f.writelines(model_lines)
        
        print(f"[VISUALIZER] Extracted first model to: {temp_path}")
        return temp_path
        
    except Exception as e:
        print(f"[VISUALIZER] Error extracting medoid model: {e}")
        print("[VISUALIZER] Falling back to original PDB file")
        return pdb_path


@dataclass
class ClassificationResult:
    """Result of classification for a single residue"""
    residue_number: int
    amino_acid: str
    csp_value: float
    is_significant: bool
    is_positive: bool
    classification: str  # 'TP', 'FP', 'TN', 'FN'


def find_system_directory(holo_pdb: str, outputs_dir: str) -> Optional[str]:
    """
    Find the system directory for a given holo PDB ID.
    
    Args:
        holo_pdb: PDB ID (e.g., '1cf4')
        outputs_dir: Path to outputs directory
        
    Returns:
        Path to system directory or None if not found
    """
    system_dir = os.path.join(outputs_dir, holo_pdb)
    if os.path.exists(system_dir):
        return system_dir
    
    # Try to find by searching subdirectories
    for item in os.listdir(outputs_dir):
        item_path = os.path.join(outputs_dir, item)
        if os.path.isdir(item_path) and item == holo_pdb:
            return item_path
    
    return None


def classify_residues_from_master_alignment(master_alignment_rows: List[Dict[str, str]], 
                                          significant_column: str = 'significant',
                                          positive_strategy: str = 'occluded_only') -> List[ClassificationResult]:
    """
    Classify residues based on master alignment data.
    
    Args:
        master_alignment_rows: Master alignment CSV rows
        significant_column: Column name for CSP significance
        positive_strategy: Strategy for defining positive residues
        
    Returns:
        List of ClassificationResult objects
    """
    results = []
    
    for row in master_alignment_rows:
        # Skip rows with missing essential data
        significant_str = row.get(significant_column, '').strip()
        holo_aa = row.get('holo_aa', '').strip()
        csp_str = row.get('csp_A', '').strip()
        
        if not significant_str or not holo_aa or not csp_str:
            continue
            
        try:
            significant = int(significant_str) == 1
            csp_value = float(csp_str)
        except ValueError:
            continue
        
        # Skip prolines (residues without NH groups)
        if holo_aa == 'P':
            continue
        
        # Get occlusion status
        is_occluded_str = row.get('is_occluded_occlusion', '').strip()
        is_occluded = is_occluded_str.lower() == 'true' if is_occluded_str else False
        
        # Get interaction status
        has_hbond_str = row.get('has_hbond_interaction', '').strip()
        has_charge_complement_str = row.get('has_charge_complement_interaction', '').strip()
        has_pi_contact_str = row.get('has_pi_contact_interaction', '').strip()
        has_hbond = has_hbond_str.lower() == 'true' if has_hbond_str else False
        has_charge_complement = has_charge_complement_str.lower() == 'true' if has_charge_complement_str else False
        has_pi_contact = has_pi_contact_str.lower() == 'true' if has_pi_contact_str else False
        is_interacting = has_hbond or has_charge_complement or has_pi_contact
        
        # Get CA distance status
        passes_ca_filter_str = row.get('passes_filter_distance', '').strip()
        passes_ca_filter = passes_ca_filter_str.lower() == 'true' if passes_ca_filter_str else False
        
        # Determine positive status based on strategy
        if positive_strategy == 'occluded_only':
            is_positive = is_occluded
        elif positive_strategy == 'occluded_or_ca':
            is_positive = is_occluded or passes_ca_filter
        elif positive_strategy == 'occluded_or_interaction':
            is_positive = is_occluded or is_interacting
        elif positive_strategy == 'occluded_or_ca_or_interaction':
            is_positive = is_occluded or passes_ca_filter or is_interacting
        else:
            # Default to occluded only
            is_positive = is_occluded
        
        # Classify based on significance and positive status
        if significant and is_positive:
            classification = 'TP'  # True Positive
        elif significant and not is_positive:
            classification = 'FP'  # False Positive
        elif not significant and not is_positive:
            classification = 'TN'  # True Negative
        else:  # not significant and is_positive
            classification = 'FN'  # False Negative
        
        # Get residue number (use sequential position as PDB residue number)
        try:
            residue_number = int(float(row.get('sequential_position', '0')))
        except ValueError:
            continue
        
        results.append(ClassificationResult(
            residue_number=residue_number,
            amino_acid=holo_aa,
            csp_value=csp_value,
            is_significant=significant,
            is_positive=is_positive,
            classification=classification
        ))
    
    return results




def create_simple_pymol_script(classifications: List[ClassificationResult], 
                              holo_pdb: str, 
                              output_path: str,
                              pdb_dir: str = "PDB_FILES",
                              print_to_stdout: bool = False) -> bool:
    """
    Create a simple PyMOL script without complex dependencies.
    Automatically finds and uses the first model for multi-model PDB files.
    
    Args:
        classifications: List of ClassificationResult objects
        holo_pdb: PDB ID
        output_path: Output script path
        pdb_dir: PDB files directory
        print_to_stdout: If True, print script to stdout instead of saving to file
        
    Returns:
        True if successful, False otherwise
    """
    pdb_file = os.path.join(pdb_dir, f"{holo_pdb}.pdb")
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file {pdb_file} does not exist")
        return False
    
    # Find medoid model if it's a multi-model PDB
    medoid_pdb_file = find_medoid_model_simple(pdb_file)
    print(f"[VISUALIZER] Using medoid model: {medoid_pdb_file}")
    
    # Count classifications
    tp_count = sum(1 for c in classifications if c.classification == 'TP')
    fp_count = sum(1 for c in classifications if c.classification == 'FP')
    tn_count = sum(1 for c in classifications if c.classification == 'TN')
    fn_count = sum(1 for c in classifications if c.classification == 'FN')
    
    # Generate PyMOL script
    lines = []
    lines.append("reinitialize")
    lines.append(f"load {medoid_pdb_file}, protein")
    lines.append("hide everything, protein")
    lines.append("show cartoon, protein")
    
    # Color peptide chain in gray (assume chain B)
    lines.append("color gray, protein and chain B")
    tp_rgb = hex_to_rgb01(classification_colors.TP)
    tn_rgb = hex_to_rgb01(classification_colors.TN)
    fp_rgb = hex_to_rgb01(classification_colors.FP)
    fn_rgb = hex_to_rgb01(classification_colors.FN)
    lines.append(f"set_color tp_color, [{tp_rgb[0]:.4f}, {tp_rgb[1]:.4f}, {tp_rgb[2]:.4f}]")
    lines.append(f"set_color tn_color, [{tn_rgb[0]:.4f}, {tn_rgb[1]:.4f}, {tn_rgb[2]:.4f}]")
    lines.append(f"set_color fp_color, [{fp_rgb[0]:.4f}, {fp_rgb[1]:.4f}, {fp_rgb[2]:.4f}]")
    lines.append(f"set_color fn_color, [{fn_rgb[0]:.4f}, {fn_rgb[1]:.4f}, {fn_rgb[2]:.4f}]")
    
    # Color protein chain residues by classification
    # TP: Bluish green (significant CSP in binding site)
    tp_residues = [str(c.residue_number) for c in classifications if c.classification == 'TP']
    if tp_residues:
        lines.append(f"# TP (Sig. CSP in Binding Site): {len(tp_residues)} residues")
        for res_num in tp_residues:
            lines.append(f"color tp_color, protein and chain A and resi {res_num}")
    
    # FP: Vermillion (significant CSP but not in binding site)
    fp_residues = [str(c.residue_number) for c in classifications if c.classification == 'FP']
    if fp_residues:
        lines.append(f"# FP (Sig. CSP -- Allosteric): {len(fp_residues)} residues")
        for res_num in fp_residues:
            lines.append(f"color fp_color, protein and chain A and resi {res_num}")
    
    # TN: Blue (low CSP and not in binding site)
    tn_residues = [str(c.residue_number) for c in classifications if c.classification == 'TN']
    if tn_residues:
        lines.append(f"# TN (low CSP -- allosteric): {len(tn_residues)} residues")
        for res_num in tn_residues:
            lines.append(f"color tn_color, protein and chain A and resi {res_num}")
    
    # FN: Reddish purple (low CSP but in binding site)
    fn_residues = [str(c.residue_number) for c in classifications if c.classification == 'FN']
    if fn_residues:
        lines.append(f"# FN (low CSP in Binding Site): {len(fn_residues)} residues")
        for res_num in fn_residues:
            lines.append(f"color fn_color, protein and chain A and resi {res_num}")
    
    # Set transparency and visual enhancements
    lines.append("set cartoon_transparency, 0.2, protein")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_ring_mode, 1")
    
    # Add colorbar/legend information
    lines.append("# CSP Classification Color Scheme:")
    lines.append(f"# TP ({classification_colors.TP}): Sig. CSP in Binding Site")
    lines.append(f"# FP ({classification_colors.FP}): Sig. CSP -- Allosteric")
    lines.append(f"# TN ({classification_colors.TN}): low CSP -- allosteric")
    lines.append(f"# FN ({classification_colors.FN}): low CSP in Binding Site")
    lines.append("# Gray: Peptide chain")
    
    # Add title and summary
    lines.append(f"# CSP Classification Analysis Summary:")
    lines.append(f"# Protein chain: A")
    lines.append(f"# Peptide chain: B")
    lines.append(f"# TP (Sig. CSP in Binding Site): {tp_count}")
    lines.append(f"# FP (Sig. CSP -- Allosteric): {fp_count}")
    lines.append(f"# TN (low CSP -- allosteric): {tn_count}")
    lines.append(f"# FN (low CSP in Binding Site): {fn_count}")
    lines.append(f"# Structure colored by CSP classification (TP/FP/TN/FN)")
    lines.append(f"# PDB file: {medoid_pdb_file}")
    
    # Write script file or print to stdout
    script_content = "\n".join(lines) + "\n"
    
    if print_to_stdout:
        print("\n" + "="*60)
        print("PYMOL SCRIPT CONTENT:")
        print("="*60)
        print(script_content)
        print("="*60)
    else:
        with open(output_path, "w") as f:
            f.write(script_content)
        print(f"[VISUALIZER] PyMOL script saved to: {output_path}")
    
    return True


def generate_pymol_session_file(script_path: str, output_path: str) -> bool:
    """
    Generate PyMOL session file from script.
    
    Args:
        script_path: Path to PyMOL script
        output_path: Path for output .pse file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Generate PyMOL session file
        pymol_cmd = [
            'pymol', '-c', '-q',  # Command line mode, quiet
            script_path,           # Run our script
            '-d', f'save {output_path}',  # Save session file
            '-d', 'quit'          # Quit PyMOL
        ]
        
        result = subprocess.run(pymol_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: PyMOL session generation failed: {result.stderr}")
            print("PyMOL script saved instead. You can run it manually in PyMOL.")
            return False
        else:
            print(f"✓ PyMOL session file saved: {output_path}")
            return True
            
    except FileNotFoundError:
        print("Warning: PyMOL not found. PyMOL script saved instead.")
        return False
    except Exception as e:
        print(f"Warning: Error generating PyMOL session: {e}")
        return False


def normalize_strategy_alias(strategy: str) -> str:
    """
    Convert short strategy aliases to full names.
    
    Args:
        strategy: Strategy string (may be alias)
        
    Returns:
        Full strategy name
    """
    alias_map = {
        'occ': 'occluded_only',
        'ca': 'occluded_or_ca', 
        'int': 'occluded_or_interaction',
        'all': 'occluded_or_ca_or_interaction'
    }
    return alias_map.get(strategy, strategy)


def normalize_threshold_alias(threshold: str) -> str:
    """
    Convert short threshold aliases to full names.
    
    Args:
        threshold: Threshold string (may be alias)
        
    Returns:
        Full threshold name
    """
    alias_map = {
        'sig': 'significant',
        '1sd': 'significant_1sd',
        '2sd': 'significant_2sd'
    }
    return alias_map.get(threshold, threshold)


def auto_load_pymol_session(pse_file_path: str) -> bool:
    """
    Automatically load the PyMOL session file in PyMOL GUI.
    
    Args:
        pse_file_path: Path to the .pse file to load
        
    Returns:
        True if successful, False otherwise
    """
    if not os.path.exists(pse_file_path):
        print(f"Error: Session file {pse_file_path} does not exist")
        return False
    
    try:
        # Try to open PyMOL GUI with the session file
        pymol_cmd = ['pymol', pse_file_path]
        
        # Run PyMOL in background (non-blocking)
        subprocess.Popen(pymol_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        print(f"✓ PyMOL session loaded: {pse_file_path}")
        print("PyMOL should open automatically with your visualization.")
        return True
        
    except FileNotFoundError:
        print("Warning: PyMOL not found in PATH. Please install PyMOL or add it to your PATH.")
        print(f"You can manually open the session file: {pse_file_path}")
        return False
    except Exception as e:
        print(f"Warning: Error loading PyMOL session: {e}")
        print(f"You can manually open the session file: {pse_file_path}")
        return False


def auto_run_pymol_script(script_path: str) -> bool:
    """
    Automatically run the PyMOL script and start PyMOL GUI.
    
    Args:
        script_path: Path to the .pml script file to run
        
    Returns:
        True if successful, False otherwise
    """
    if not os.path.exists(script_path):
        print(f"Error: Script file {script_path} does not exist")
        return False
    
    try:
        # Try to run PyMOL with the script file
        pymol_cmd = ['pymol', script_path]
        
        # Run PyMOL in background (non-blocking)
        subprocess.Popen(pymol_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        print(f"✓ PyMOL script executed: {script_path}")
        print("PyMOL should open automatically with your visualization.")
        return True
        
    except FileNotFoundError:
        print("Warning: PyMOL not found in PATH. Please install PyMOL or add it to your PATH.")
        print(f"You can manually run the script in PyMOL: {script_path}")
        return False
    except Exception as e:
        print(f"Warning: Error running PyMOL script: {e}")
        print(f"You can manually run the script in PyMOL: {script_path}")
        return False


def generate_pymol_visualization(holo_pdb: str, 
                                positive_strategy: str, 
                                significance_threshold: str,
                                output_path: str,
                                outputs_dir: str = "outputs",
                                pdb_dir: str = "PDB_FILES",
                                print_script: bool = False,
                                no_auto_run: bool = False) -> bool:
    """
    Generate PyMOL visualization based on confusion matrix analysis.
    
    Args:
        holo_pdb: PDB ID (e.g., '1cf4')
        positive_strategy: Strategy for defining positive residues
        significance_threshold: Significance threshold ('significant', 'significant_1sd', 'significant_2sd')
        output_path: Path for output .pse file
        outputs_dir: Path to outputs directory
        pdb_dir: Path to PDB files directory
        print_script: If True, print PyMOL script to stdout instead of saving
        no_auto_run: If True, don't automatically run PyMOL
        
    Returns:
        True if successful, False otherwise
    """
    print(f"[VISUALIZER] Generating PyMOL visualization for {holo_pdb}")
    print(f"[VISUALIZER] Strategy: {positive_strategy}, Threshold: {significance_threshold}")
    
    # Find system directory
    system_dir = find_system_directory(holo_pdb, outputs_dir)
    if not system_dir:
        print(f"Error: Could not find system directory for {holo_pdb}")
        return False
    
    # Read master alignment data
    master_alignment_path = os.path.join(system_dir, "master_alignment.csv")
    if not os.path.exists(master_alignment_path):
        print(f"Error: Could not find master_alignment.csv in {system_dir}")
        return False
    
    master_alignment_rows = read_master_alignment(master_alignment_path)
    if not master_alignment_rows:
        print(f"Error: Could not read master_alignment.csv")
        return False
    
    print(f"[VISUALIZER] Loaded {len(master_alignment_rows)} rows from master alignment")
    
    # Classify residues
    classifications = classify_residues_from_master_alignment(
        master_alignment_rows, 
        significance_threshold, 
        positive_strategy
    )
    
    if not classifications:
        print("Error: No valid classifications found")
        return False
    
    # Count classifications
    tp_count = sum(1 for c in classifications if c.classification == 'TP')
    fp_count = sum(1 for c in classifications if c.classification == 'FP')
    tn_count = sum(1 for c in classifications if c.classification == 'TN')
    fn_count = sum(1 for c in classifications if c.classification == 'FN')
    
    print(f"[VISUALIZER] Classification counts - TP: {tp_count}, FP: {fp_count}, TN: {tn_count}, FN: {fn_count}")
    
    # Create temporary script file
    script_path = output_path.replace('.pse', '.pml')
    
    # Generate PyMOL script
    success = create_simple_pymol_script(classifications, holo_pdb, script_path, pdb_dir, print_script)
    
    if not success:
        print(f"[VISUALIZER] ✗ Failed to generate PyMOL script")
        return False
    
    # If printing script, don't generate session file
    if print_script:
        print(f"[VISUALIZER] ✓ PyMOL script printed to stdout")
        return True
    
    # Try to generate session file
    session_success = generate_pymol_session_file(script_path, output_path)
    
    if session_success:
        print(f"[VISUALIZER] ✓ PyMOL session file saved: {output_path}")
        # Clean up script file
        try:
            os.unlink(script_path)
        except:
            pass
        # Automatically run the session file if not printing script and auto-run is enabled
        if not print_script and not no_auto_run:
            print(f"[VISUALIZER] Automatically running PyMOL session...")
            auto_load_pymol_session(output_path)
    else:
        print(f"[VISUALIZER] ✓ PyMOL script saved: {script_path}")
        # Automatically run the script if auto-run is enabled
        if not no_auto_run:
            print(f"[VISUALIZER] Automatically running PyMOL script...")
            auto_run_pymol_script(script_path)
    
    # Print summary metrics
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    print(f"[VISUALIZER] Summary metrics:")
    print(f"[VISUALIZER]   Precision: {precision:.3f}")
    print(f"[VISUALIZER]   Recall: {recall:.3f}")
    print(f"[VISUALIZER]   F1-Score: {f1_score:.3f}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate PyMOL visualization based on confusion matrix analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/confusion_matrix_visualizer.py -p 1cf4 -s occ -t sig
  python scripts/confusion_matrix_visualizer.py -p 1d5g -s int -t 1sd
  python scripts/confusion_matrix_visualizer.py -p 1cf4 -s all -t 2sd -o my_viz.pse
  python scripts/confusion_matrix_visualizer.py -p 1cf4 -s occ -t sig --print-script
  python scripts/confusion_matrix_visualizer.py -p 1cf4 -s occ -t sig --no-auto-run

Available strategies (with short aliases):
  - occ (occluded_only): Only occluded residues are positive
  - ca (occluded_or_ca): Occluded OR CA distance filtered residues are positive
  - int (occluded_or_interaction): Occluded OR interacting residues are positive
  - all (occluded_or_ca_or_interaction): Occluded OR CA distance OR interacting residues are positive

Available thresholds (with short aliases):
  - sig (significant): Standard significance threshold
  - 1sd (significant_1sd): 1 standard deviation threshold
  - 2sd (significant_2sd): 2 standard deviation threshold

Note: PyMOL will automatically open with your visualization unless --no-auto-run is specified.
        """
    )
    
    parser.add_argument("-p", "--holo-pdb", required=True,
                       help="Holo PDB ID (e.g., '1cf4')")
    parser.add_argument("-s", "--strategy", required=True,
                       choices=['occluded_only', 'occluded_or_ca', 'occluded_or_interaction', 'occluded_or_ca_or_interaction', 
                               'occ', 'ca', 'int', 'all'],
                       help="Strategy for defining positive residues (aliases: occ=occluded_only, ca=occluded_or_ca, int=occluded_or_interaction, all=occluded_or_ca_or_interaction)")
    parser.add_argument("-t", "--threshold", required=True,
                       choices=['significant', 'significant_1sd', 'significant_2sd', 'sig', '1sd', '2sd'],
                       help="Significance threshold for CSP classification (aliases: sig=significant, 1sd=significant_1sd, 2sd=significant_2sd)")
    parser.add_argument("-o", "--output", default="temp.pse",
                       help="Output path for .pse file (default: temp.pse)")
    parser.add_argument("--outputs-dir", default="outputs",
                       help="Path to outputs directory (default: outputs)")
    parser.add_argument("--pdb-dir", default="PDB_FILES",
                       help="Path to PDB files directory (default: PDB_FILES)")
    parser.add_argument("--no-auto-load", action="store_true",
                       help="Don't automatically load the .pse file in PyMOL")
    parser.add_argument("--print-script", action="store_true",
                       help="Print the PyMOL script content to stdout instead of saving to file")
    parser.add_argument("--no-auto-run", action="store_true",
                       help="Don't automatically run PyMOL (overrides --no-auto-load)")
    
    args = parser.parse_args()
    
    # Normalize aliases to full names
    strategy = normalize_strategy_alias(args.strategy)
    threshold = normalize_threshold_alias(args.threshold)
    
    # Validate inputs
    if not os.path.exists(args.outputs_dir):
        print(f"Error: Outputs directory {args.outputs_dir} does not exist")
        return 1
    
    if not os.path.exists(args.pdb_dir):
        print(f"Error: PDB directory {args.pdb_dir} does not exist")
        return 1
    
    pdb_file = os.path.join(args.pdb_dir, f"{args.holo_pdb}.pdb")
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file {pdb_file} does not exist")
        return 1
    
    # Generate visualization
    success = generate_pymol_visualization(
        args.holo_pdb,
        strategy,
        threshold,
        args.output,
        args.outputs_dir,
        args.pdb_dir,
        args.print_script,
        args.no_auto_run
    )
    
    # Note: Automatic PyMOL execution is now handled within generate_pymol_visualization
    # The --no-auto-run flag is passed to that function to control execution
    
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
