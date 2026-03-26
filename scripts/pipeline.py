"""
CLI orchestrator for the CSP pipeline.

Usage examples:
  python scripts/pipeline.py --input data/CSP_UBQ.csv --out outputs --workers 4
  python scripts/pipeline.py --input data/CSP_UBQ.csv --ids 18251,34688 --out outputs
  python scripts/pipeline.py --input data/CSP_UBQ.csv --holo-pdb 1cf4 --out outputs
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import List, Dict, Optional, Any, Callable, Iterator, TextIO
from pathlib import Path

# Support running as a script (python scripts/pipeline.py) or module (python -m scripts.pipeline)
try:
    from .config import paths, concurrency, ensure_directories, sasa_analysis, ca_distance_analysis, Referencing as _Referencing  # type: ignore
    from .bmrb_io import fetch_bmrb, parse_sequence_and_shifts, parse_sequence_and_shifts_from_saveframes
    from .align import align_global
    from .csp import compute_csp_A, compute_csp_multiple_saveframes, compute_csp_from_aligned_sequences, compute_csp_multiple_saveframes_ca
    from .rcsb_io import fetch_pdb, parse_pdb_sequences
    from .visualize import write_pymol_color_csp_mask_script, plot_csp_histogram, write_pymol_occlusion_script, write_pymol_combined_script, write_pymol_delta_sasa_script, write_pymol_session_file, create_pdb_with_delta_sasa_bfactors, write_pymol_csp_heatmap_script, write_pymol_csp_session_file, plot_csp_classification_bars, write_pymol_csp_classification_script, write_pymol_csp_classification_session_file, plot_per_atom_classification_panels
    from .HSQC_visualize import plot_hsqc_variants
    from .sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv, get_occlusion_summary
    from .interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv, get_interaction_summary, compute_ca_distance_filter, write_ca_distance_csv, get_ca_distance_summary, compute_nn_distance_filter, write_nn_distance_csv, get_nn_distance_summary, compute_min_atom_distance_filter, write_any_atom_distance_csv, get_any_atom_distance_summary
    from .merge_csv import merge_all_csv_files
    from .analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from .case_study import generate_case_study_figure
    from .case_study_2 import generate_case_study_2_figure
    from .annotate_csp_csv_metadata import annotate_csv_with_ec_and_scope
    from .confusion_matrix_analysis import generate_confusion_matrix_per_system
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import paths, concurrency, ensure_directories, sasa_analysis, ca_distance_analysis, Referencing as _Referencing  # type: ignore
    from scripts.bmrb_io import fetch_bmrb, parse_sequence_and_shifts, parse_sequence_and_shifts_from_saveframes
    from scripts.align import align_global
    from scripts.csp import compute_csp_A, compute_csp_multiple_saveframes, compute_csp_from_aligned_sequences, compute_csp_multiple_saveframes_ca
    from scripts.rcsb_io import fetch_pdb, parse_pdb_sequences
    from scripts.visualize import write_pymol_color_csp_mask_script, plot_csp_histogram, write_pymol_occlusion_script, write_pymol_combined_script, write_pymol_delta_sasa_script, write_pymol_session_file, create_pdb_with_delta_sasa_bfactors, write_pymol_csp_heatmap_script, write_pymol_csp_session_file, plot_csp_classification_bars, write_pymol_csp_classification_script, write_pymol_csp_classification_session_file, plot_per_atom_classification_panels
    from scripts.HSQC_visualize import plot_hsqc_variants
    from scripts.sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv, get_occlusion_summary
    from scripts.interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv, get_interaction_summary, compute_ca_distance_filter, write_ca_distance_csv, get_ca_distance_summary, compute_nn_distance_filter, write_nn_distance_csv, get_nn_distance_summary, compute_min_atom_distance_filter, write_any_atom_distance_csv, get_any_atom_distance_summary
    from scripts.merge_csv import merge_all_csv_files
    from scripts.analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from scripts.case_study import generate_case_study_figure
    from scripts.case_study_2 import generate_case_study_2_figure
    from scripts.annotate_csp_csv_metadata import annotate_csv_with_ec_and_scope
    from scripts.confusion_matrix_analysis import generate_confusion_matrix_per_system



def _console_line(message: str) -> None:
    """Write a concise progress line to real stdout."""
    print(message, file=sys.__stdout__, flush=True)


def _append_log(log_path: str, message: str) -> None:
    """Append one line to a log file."""
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(message + "\n")


@contextmanager
def _capture_to_log(log_path: str) -> Iterator[TextIO]:
    """Capture stdout/stderr to a per-step log file."""
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write(f"\n[{datetime.now().isoformat(timespec='seconds')}] capture_start\n")
        lf.flush()
        with redirect_stdout(lf), redirect_stderr(lf):
            yield lf
        lf.write(f"[{datetime.now().isoformat(timespec='seconds')}] capture_end\n")
        lf.flush()


def _run_logged(log_path: str, func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Run a callable while capturing stdout/stderr into log_path."""
    with _capture_to_log(log_path):
        return func(*args, **kwargs)


def _emit_warning(message: str, log_path: Optional[str] = None) -> None:
    """Emit warning/error text to stdout and optional log file."""
    _console_line(message)
    if log_path:
        _append_log(log_path, message)


def process_row(
    row: Dict[str, str],
    out_dir: str,
    sasa_args: Dict = None,
    csp_threshold_args: Dict = None,
    interaction_args: Dict = None,
    binary_mode: bool = False,
    include_alternative_thresholds: bool = False,
    directory_suffix: Optional[str] = None,
    generate_case_study: bool = True,
    include_numeric_residue_ticks: bool = False,
    force_case_study_view_reset: bool = False,
) -> None:
    apo_bmrb = (row.get("apo_bmrb") or "").strip()
    holo_bmrb = (row.get("holo_bmrb") or "").strip()
    holo_pdb = (row.get("holo_pdb") or "").strip() or (row.get("holo_pdb_id") or "").strip()
    if not apo_bmrb or not holo_bmrb or not holo_pdb:
        return

    # IO: ensure per-target output dir
    # If directory_suffix is provided, append it to create unique directories (e.g., 2mur_1, 2mur_2)
    # IMPORTANT: When duplicates exist, directory_suffix should ALWAYS be set by the caller
    # to prevent creating non-suffixed directories
    if directory_suffix:
        tgt_dir = os.path.join(out_dir, f"{holo_pdb}_{directory_suffix}")
    else:
        tgt_dir = os.path.join(out_dir, holo_pdb)
    os.makedirs(tgt_dir, exist_ok=True)
    logs_dir = os.path.join(tgt_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_files = {
        "fetch_parse_bmrb": os.path.join(logs_dir, "01_fetch_parse_bmrb.txt"),
        "compute_csp": os.path.join(logs_dir, "02_compute_csp.txt"),
        "structure_pdb_alignment": os.path.join(logs_dir, "03_structure_pdb_alignment.txt"),
        "sasa": os.path.join(logs_dir, "04_sasa_occlusion.txt"),
        "interaction": os.path.join(logs_dir, "05_interaction_filter.txt"),
        "ca_distance": os.path.join(logs_dir, "06_ca_distance_filter.txt"),
        "nn_distance": os.path.join(logs_dir, "07_nn_distance_filter.txt"),
        "any_atom_distance": os.path.join(logs_dir, "08_any_atom_distance_filter.txt"),
        "tables_outputs": os.path.join(logs_dir, "09_tables_and_visualizations.txt"),
        "master_csv": os.path.join(logs_dir, "10_master_csv.txt"),
        "case_study": os.path.join(logs_dir, "11_case_study.txt"),
    }
    target_label = os.path.basename(tgt_dir)
    
    # Normalize tgt_dir to be relative to project root for PyMOL scripts (e.g., "./outputs/2mur_1/")
    # Convert to relative path if it's absolute, ensure it starts with ./ and ends with /
    tgt_dir_for_pymol = tgt_dir
    if os.path.isabs(tgt_dir):
        # Get project root (parent of scripts directory)
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        tgt_dir_for_pymol = os.path.relpath(tgt_dir, project_root)
    # Normalize path separators and ensure it starts with ./ and ends with /
    tgt_dir_for_pymol = tgt_dir_for_pymol.replace('\\', '/')
    if not tgt_dir_for_pymol.startswith('./') and not tgt_dir_for_pymol.startswith('/'):
        tgt_dir_for_pymol = './' + tgt_dir_for_pymol
    if not tgt_dir_for_pymol.endswith('/'):
        tgt_dir_for_pymol += '/'
    
    _console_line(f"[PIPE] [{target_label}] Start")
    _append_log(
        log_files["fetch_parse_bmrb"],
        f"[PIPE] Start apo_bmrb={apo_bmrb} holo_bmrb={holo_bmrb} holo_pdb={holo_pdb}",
    )

    # Fetch and parse BMRB entries
    _console_line(f"[PIPE] [{target_label}] Step 1/5 Fetch+parse BMRB")
    apo_star = _run_logged(log_files["fetch_parse_bmrb"], fetch_bmrb, apo_bmrb)
    holo_star = _run_logged(log_files["fetch_parse_bmrb"], fetch_bmrb, holo_bmrb)
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Parsing BMRB apo={os.path.basename(apo_star)} holo={os.path.basename(holo_star)}")
    
    # Parse all sequences from multiple saveframes
    apo_sequences = _run_logged(log_files["fetch_parse_bmrb"], parse_sequence_and_shifts_from_saveframes, apo_star)
    holo_sequences = _run_logged(log_files["fetch_parse_bmrb"], parse_sequence_and_shifts_from_saveframes, holo_star)
    
    if not apo_sequences:
        _emit_warning(f"[PIPE] ERROR: No valid apo sequences found in {apo_bmrb}", log_files["fetch_parse_bmrb"])
        return
    if not holo_sequences:
        _emit_warning(f"[PIPE] ERROR: No valid holo sequences found in {holo_bmrb}", log_files["fetch_parse_bmrb"])
        return

    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Found {len(apo_sequences)} apo sequences, {len(holo_sequences)} holo sequences")

    # Compute CSPs using best alignment from all sequence pairs
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print("[PIPE] Computing CSP(A) for all sequence pairs")
    # Determine referencing method from config
    ref_method = None
    try:
        ref_method = _Referencing().method
    except Exception:
        ref_method = None
    _console_line(f"[PIPE] [{target_label}] Step 2/5 Compute CSP")
    results = _run_logged(
        log_files["compute_csp"],
        compute_csp_multiple_saveframes,
        apo_sequences,
        holo_sequences,
        apo_bmrb,
        holo_bmrb,
        holo_pdb,
        csp_threshold_args,
        referencing_method=ref_method,
        grid_params=None,
        target_id=target_label,
        output_root=out_dir,
    )
    
    if not results:
        _emit_warning(f"[PIPE] ERROR: No valid CSPs computed for {apo_bmrb} vs {holo_bmrb}", log_files["compute_csp"])
        return

    # Check if CA shifts are available in both apo and holo
    has_apo_ca = any(len(seq[3]) > 0 for seq in apo_sequences)  # CA_shifts is at index 3
    has_holo_ca = any(len(seq[3]) > 0 for seq in holo_sequences)  # CA_shifts is at index 3
    
    # Compute CA-inclusive CSPs (separate analysis) only if CA shifts are available in both
    results_ca = []
    if has_apo_ca and has_holo_ca:
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print("[PIPE] CA shifts found in both apo and holo, computing CA-inclusive CSPs")
        results_ca = _run_logged(
            log_files["compute_csp"],
            compute_csp_multiple_saveframes_ca,
            apo_sequences,
            holo_sequences,
            apo_bmrb,
            holo_bmrb,
            holo_pdb,
            csp_threshold_args,
            referencing_method=ref_method,
            grid_params=None,
            target_id=target_label,
            output_root=out_dir,
        )
    else:
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            if not has_apo_ca and not has_holo_ca:
                print("[PIPE] No CA shifts found in apo or holo, skipping CA-inclusive CSP analysis")
            elif not has_apo_ca:
                print("[PIPE] No CA shifts found in apo, skipping CA-inclusive CSP analysis")
            else:
                print("[PIPE] No CA shifts found in holo, skipping CA-inclusive CSP analysis")

    # Download holo PDB (in parallel with above in future; simple now)
    _console_line(f"[PIPE] [{target_label}] Step 3/5 Structure + interaction analyses")
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Fetching PDB {holo_pdb}")
    pdb_path = _run_logged(log_files["structure_pdb_alignment"], fetch_pdb, holo_pdb)
    chains = _run_logged(log_files["structure_pdb_alignment"], parse_pdb_sequences, pdb_path)

    # Find the best apo-holo alignment to determine which holo sequence was used
    # This matches the logic in compute_csp_multiple_saveframes
    best_alignment_score = float('-inf')
    best_holo_sequence = None
    best_apo_sequence = None
    best_aligned_apo = None
    best_aligned_holo = None
    best_mapping = None
    
    for apo_seq, _, _, _, apo_saveframe in apo_sequences:
        for holo_seq, _, _, _, holo_saveframe in holo_sequences:
            aligned_apo, aligned_holo, mapping, alignment_score = align_global(apo_seq, holo_seq)
            if alignment_score > best_alignment_score:
                best_alignment_score = alignment_score
                best_holo_sequence = holo_seq
                best_apo_sequence = apo_seq
                best_aligned_apo = aligned_apo
                best_aligned_holo = aligned_holo
                best_mapping = mapping
    
    # Save sequence alignment to text file
    if best_aligned_apo and best_aligned_holo:
        alignment_file_path = os.path.join(tgt_dir, "csp_sequence_alignment.txt")
        try:
            with open(alignment_file_path, "w") as f:
                f.write("Sequence Alignment: Apo vs Holo Chemical Shift Lists\n")
                f.write("=" * 80 + "\n\n")
                f.write(f"Apo BMRB ID: {apo_bmrb}\n")
                f.write(f"Holo BMRB ID: {holo_bmrb}\n")
                f.write(f"Alignment Score: {best_alignment_score:.2f}\n")
                f.write(f"Number of Mapped Pairs: {len(best_mapping) if best_mapping else 0}\n\n")
                
                # Calculate statistics
                matches = sum(1 for a, b in zip(best_aligned_apo, best_aligned_holo) if a == b and a != '-')
                gaps_apo = sum(1 for c in best_aligned_apo if c == '-')
                gaps_holo = sum(1 for c in best_aligned_holo if c == '-')
                mismatches = len(best_aligned_apo) - matches - gaps_apo - gaps_holo
                total_aligned = len(best_aligned_apo) - gaps_apo - gaps_holo
                identity = (matches / total_aligned * 100) if total_aligned > 0 else 0.0
                
                f.write("Alignment Statistics:\n")
                f.write(f"  Matches:    {matches:4d}\n")
                f.write(f"  Mismatches: {mismatches:4d}\n")
                f.write(f"  Gaps (apo): {gaps_apo:4d}\n")
                f.write(f"  Gaps (holo): {gaps_holo:4d}\n")
                f.write(f"  Identity:   {identity:.1f}%\n\n")
                f.write("=" * 80 + "\n\n")
                
                # Write alignment in blocks of 80 characters
                line_width = 80
                match_string = ''.join('|' if a == b and a != '-' else ' ' for a, b in zip(best_aligned_apo, best_aligned_holo))
                
                for i in range(0, len(best_aligned_apo), line_width):
                    chunk_apo = best_aligned_apo[i:i+line_width]
                    chunk_holo = best_aligned_holo[i:i+line_width]
                    chunk_match = match_string[i:i+line_width]
                    
                    start_pos = i + 1
                    end_pos = min(i + line_width, len(best_aligned_apo))
                    
                    f.write(f"Apo  {start_pos:4d}-{end_pos:4d}: {chunk_apo}\n")
                    f.write(f"      {' ' * 10} {chunk_match}\n")
                    f.write(f"Holo {start_pos:4d}-{end_pos:4d}: {chunk_holo}\n")
                    f.write("\n")
            
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ Sequence alignment saved: {alignment_file_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to save sequence alignment: {e}", log_files["structure_pdb_alignment"])
    
    # Identify receptor chain by aligning the apo sequence from best apo-holo alignment with each PDB chain.
    # The receptor is the chain that aligns most closely with the apo chemical shift list (user requirement).
    # Each entry computes independently; no cache, since duplicate holo_pdb entries can have different apo proteins.
    receptor_chain = None
    ligand_chain = None
    # Use apo sequence for chain selection (receptor = chain matching apo); fallback to holo if apo unavailable
    ref_sequence = best_apo_sequence or best_holo_sequence
    best_score = -1
    best_chain_length = 0
    best_length_match = float('inf')
    best_num_matches = -1
    chain_scores = {}
    if ref_sequence:
        ref_seq_length = len(ref_sequence)
        for ch, seq in chains.items():
            aligned_ref, aligned_ch, ch_map, _ = align_global(ref_sequence, seq)
            num_mapped = len(ch_map)
            chain_len = len(seq)
            length_diff = abs(chain_len - ref_seq_length)
            # Chain coverage: fraction of chain residues mapped to apo.
            match_fraction = num_mapped / chain_len if chain_len > 0 else 0
            # Ref coverage: fraction of apo (reference) covered by this chain. Receptor is the chain
            # that covers most of the apo; prefer over match_fraction to avoid short chains winning.
            ref_coverage = num_mapped / ref_seq_length if ref_seq_length > 0 else 0
            # Primary score: ref_coverage favors chains that cover the apo; length_diff penalizes size mismatch.
            score = ref_coverage * 1000 - length_diff
            # Count identity (matching residues) for tie-breaking when sequences are similar
            num_matches = sum(1 for a, b in zip(aligned_ref, aligned_ch) if a == b and a != "-")
            chain_scores[ch] = {"score": score, "num_mapped": num_mapped, "length_diff": length_diff, "chain_len": len(seq), "num_matches": num_matches, "ref_coverage": ref_coverage, "match_fraction": match_fraction}
            # Score includes length penalty; use num_matches for exact ties
            if score > best_score or (score >= best_score and num_matches > best_num_matches):
                best_score = score
                receptor_chain = ch
                best_chain_length = len(seq)
                best_length_match = length_diff
                best_num_matches = num_matches
    else:
        # Fallback: use any holo sequence if we couldn't determine the best one
        for ch, seq in chains.items():
            for holo_seq, _, _, _, _ in holo_sequences:
                _, _, ch_map, _ = align_global(holo_seq, seq)
                num_mapped = len(ch_map)
                if num_mapped > best_score or (num_mapped >= best_score * 0.95 and len(seq) > best_chain_length):
                    best_score = num_mapped
                    receptor_chain = ch
                    best_chain_length = len(seq)

    if receptor_chain:
        other_chains = [ch for ch in chains.keys() if ch != receptor_chain]
        if other_chains:
            ligand_chain = other_chains[0]
        elif len(chains) == 1:
            _emit_warning(
                "[PIPE] WARNING: Only one chain found in PDB, cannot identify ligand chain",
                log_files["structure_pdb_alignment"],
            )

    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Selected receptor chain={receptor_chain} (score={best_score:.2f}, length_match={best_length_match}) based on apo sequence alignment")
        if ligand_chain:
            print(f"[PIPE] Selected ligand chain={ligand_chain} (length={len(chains[ligand_chain])})")

    # Perform SASA occlusion analysis
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Performing SASA occlusion analysis")
    
    # Use SASA parameters from command line or defaults
    if sasa_args is None:
        sasa_args = {
            'probe_radius': 0.0,
            'nh_mode': 'sum',
            'mode': 'residue',
            'threshold': sasa_analysis.sasa_threshold
        }
    
    sasa_results = _run_logged(
        log_files["sasa"],
        compute_sasa_occlusion,
        pdb_path,
        sasa_threshold=sasa_args['threshold'],
        probe_radius=sasa_args['probe_radius'],
        nh_mode=sasa_args['nh_mode'],
        mode=sasa_args['mode'],
        receptor_chain_id=receptor_chain,
        ligand_chain_id=ligand_chain
    )
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(get_occlusion_summary(sasa_results))
    
    # Write occlusion analysis CSV
    occlusion_csv_path = os.path.join(tgt_dir, "occlusion_analysis.csv")
    _run_logged(log_files["sasa"], write_occlusion_analysis_csv, sasa_results['residue_info'], occlusion_csv_path)

    # Perform interaction analysis
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Performing interaction analysis")
    
    # Use interaction parameters from command line or defaults
    if interaction_args is None:
        interaction_args = {
            'distance_threshold': 4.5,
            'pi_distance_threshold': 6.0,
        }
    else:
        interaction_args.setdefault('pi_distance_threshold', 6.0)
    
    interaction_results = _run_logged(
        log_files["interaction"],
        compute_interaction_filter,
        pdb_path,
        distance_threshold=interaction_args['distance_threshold'],
        pi_distance_threshold=interaction_args['pi_distance_threshold'],
        receptor_chain_id=receptor_chain,
        ligand_chain_id=ligand_chain
    )
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(get_interaction_summary(interaction_results))
    # Write interaction analysis CSV
    interaction_csv_path = os.path.join(tgt_dir, "interaction_filter.csv")
    _run_logged(log_files["interaction"], write_interaction_analysis_csv, interaction_results['residue_info'], interaction_csv_path)

    # Perform CA distance filter analysis
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Performing CA distance filter analysis")
    
    ca_distance_results = _run_logged(
        log_files["ca_distance"],
        compute_ca_distance_filter,
        pdb_path,
        distance_threshold=ca_distance_analysis.ca_distance_threshold,
        receptor_chain_id=receptor_chain,
        ligand_chain_id=ligand_chain
    )
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(get_ca_distance_summary(ca_distance_results))
    

    # Write CA distance filter CSV
    ca_distance_csv_path = os.path.join(tgt_dir, "ca_distance_filter.csv")
    _run_logged(log_files["ca_distance"], write_ca_distance_csv, ca_distance_results['residue_info'], ca_distance_csv_path)

    # Perform N-N distance filter analysis
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Performing N-N distance filter analysis")
    nn_distance_results = _run_logged(
        log_files["nn_distance"],
        compute_nn_distance_filter,
        pdb_path,
        distance_threshold=ca_distance_analysis.ca_distance_threshold,
        receptor_chain_id=receptor_chain,
        ligand_chain_id=ligand_chain
    )
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(get_nn_distance_summary(nn_distance_results))
    nn_distance_csv_path = os.path.join(tgt_dir, "nn_distance_filter.csv")
    _run_logged(log_files["nn_distance"], write_nn_distance_csv, nn_distance_results['residue_info'], nn_distance_csv_path)

    # Perform any-atom distance filter analysis
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Performing any-atom distance filter analysis")
    any_atom_results = _run_logged(
        log_files["any_atom_distance"],
        compute_min_atom_distance_filter,
        pdb_path,
        distance_threshold=ca_distance_analysis.ca_distance_threshold,
        receptor_chain_id=receptor_chain,
        ligand_chain_id=ligand_chain
    )
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(get_any_atom_distance_summary(any_atom_results))
    any_atom_csv_path = os.path.join(tgt_dir, "any_atom_distance_filter.csv")
    _run_logged(log_files["any_atom_distance"], write_any_atom_distance_csv, any_atom_results['residue_info'], any_atom_csv_path)

    # Create mapping from residue numbers to occlusion data
    occlusion_map = {}
    for info in sasa_results['residue_info']:
        occlusion_map[info['residue_number']] = {
            'delta_sasa': info['delta_sasa'],
            'is_occluded': info['is_occluded']
        }
    
    # Create sequence alignment mapping between CSP sequential positions and PDB residue numbers
    # Extract sequences from CSP results (sequential positions)
    csp_sequence = []
    csp_positions = []
    for r in results:
        if r.holo_aa and r.holo_aa != 'P':  # Skip prolines
            csp_sequence.append(r.holo_aa)
            csp_positions.append(r.holo_index)
    
    # Extract sequences from SASA results (PDB residue numbers)
    occlusion_sequence = []
    occlusion_positions = []
    for info in sasa_results['residue_info']:
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
    csp_seq_str = ''.join(csp_sequence)
    occlusion_seq_str = ''.join(occlusion_sequence)
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] CSP sequence length: {len(csp_seq_str)}")
        print(f"[PIPE] Occlusion sequence length: {len(occlusion_seq_str)}")
        print(f"[PIPE] CSP sequence: {csp_seq_str[:50]}...")
        print(f"[PIPE] Occlusion sequence: {occlusion_seq_str[:50]}...")
    
    # Perform global alignment
    aligned_csp, aligned_occlusion, mapping, alignment_score = align_global(csp_seq_str, occlusion_seq_str)
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Alignment score: {alignment_score}")
        print(f"[PIPE] Mapped pairs: {len(mapping)}")
    
    # Create mapping from sequential position to PDB residue number
    sequential_to_pdb_map = {}
    for csp_pos, occlusion_pos in mapping:
        if csp_pos <= len(csp_positions) and occlusion_pos <= len(occlusion_positions):
            sequential_position = csp_positions[csp_pos - 1]  # Convert to 1-based
            pdb_residue_number = occlusion_positions[occlusion_pos - 1]  # Convert to 1-based
            sequential_to_pdb_map[sequential_position] = pdb_residue_number
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Created position mapping with {len(sequential_to_pdb_map)} entries")
        print(f"[PIPE] First 10 mappings: {list(sequential_to_pdb_map.items())[:10]}")

    # Build a quick lookup for CA data by holo_index (if available)
    ca_by_holo_index: Dict[int, Dict[str, Optional[float]]] = {}
    for r in results_ca or []:
        try:
            idx = int(r.holo_index)
        except Exception:
            continue
        ca_by_holo_index[idx] = {
            "CA_apo": getattr(r, "CA_apo", None),
            "CA_holo": getattr(r, "CA_holo", None),
            "CA_offset": getattr(r, "CA_offset", None),
        }

    _console_line(f"[PIPE] [{target_label}] Step 4/5 Write tables + visualizations")

    # Write outputs for H/N-only CSPs, now including CA columns (when available)
    table_path = os.path.join(tgt_dir, "csp_table.csv")
    with open(table_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "apo_bmrb","holo_bmrb","holo_pdb","chain","apo_resi","apo_aa","holo_resi","holo_aa",
            "H_apo","N_apo","CA_apo",
            "H_holo","N_holo","CA_holo",
            "H_holo_original","N_holo_original",
            "H_offset","N_offset","CA_offset",
            "dH","dN","csp_A","csp_z","significant","significant_1sd","significant_2sd",
            "delta_sasa","occluded",
        ])
        for r in results:
            # Get corresponding PDB residue number from alignment
            pdb_residue_number = sequential_to_pdb_map.get(r.holo_index)
            if pdb_residue_number is None:
                if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                    _emit_warning(
                        f"[PIPE WARNING] No alignment found for sequential position {r.holo_index}",
                        log_files["tables_outputs"],
                    )
                # Default to no occlusion data if no alignment found
                occlusion_data = {'delta_sasa': '', 'is_occluded': False}
            else:
                # Get occlusion data for this PDB residue
                occlusion_data = occlusion_map.get(pdb_residue_number, {'delta_sasa': '', 'is_occluded': False})

            ca_info = ca_by_holo_index.get(int(r.holo_index), {}) if r.holo_index is not None else {}
            CA_apo = ca_info.get("CA_apo")
            CA_holo = ca_info.get("CA_holo")
            CA_offset = ca_info.get("CA_offset")

            w.writerow([
                apo_bmrb, holo_bmrb, holo_pdb, receptor_chain or "?", r.apo_index, r.apo_aa, r.holo_index, r.holo_aa,
                f"{r.H_apo:.4f}" if r.H_apo is not None else "",
                f"{r.N_apo:.4f}" if r.N_apo is not None else "",
                f"{CA_apo:.4f}" if CA_apo is not None else "",
                f"{r.H_holo:.4f}" if r.H_holo is not None else "",
                f"{r.N_holo:.4f}" if r.N_holo is not None else "",
                f"{CA_holo:.4f}" if CA_holo is not None else "",
                f"{r.H_holo_original:.4f}" if r.H_holo_original is not None else "",
                f"{r.N_holo_original:.4f}" if r.N_holo_original is not None else "",
                f"{r.H_offset:.4f}" if r.H_offset is not None else "",
                f"{r.N_offset:.4f}" if r.N_offset is not None else "",
                f"{CA_offset:.4f}" if CA_offset is not None else "",
                f"{r.dH:.4f}" if r.dH is not None else "",
                f"{r.dN:.4f}" if r.dN is not None else "",
                f"{r.csp_A:.4f}" if r.csp_A is not None else "",
                f"{r.z_score:.4f}" if r.z_score is not None else "",
                int(bool(r.significant)) if r.significant is not None else "",
                int(bool(r.significant_1sd)) if r.significant_1sd is not None else "",
                int(bool(r.significant_2sd)) if r.significant_2sd is not None else "",
                f"{occlusion_data['delta_sasa']:.4f}" if occlusion_data['delta_sasa'] != '' else "",
                int(bool(occlusion_data['is_occluded'])) if occlusion_data['is_occluded'] is not None else "",
            ])

    # Write outputs for CA-inclusive CSPs (separate analysis) if available
    if results_ca:
        table_ca_path = os.path.join(tgt_dir, "csp_table_CA.csv")
        with open(table_ca_path, "w", newline="") as f_ca:
            w_ca = csv.writer(f_ca)
            w_ca.writerow([
                "apo_bmrb","holo_bmrb","holo_pdb","chain","apo_resi","apo_aa","holo_resi","holo_aa",
                "H_apo","N_apo","CA_apo",
                "H_holo","N_holo","CA_holo",
                "H_holo_original","N_holo_original","CA_holo_original",
                "H_offset","N_offset","CA_offset",
                "dH","dN","dCA",
                "csp_CA","csp_CA_z","csp_CA_significant","csp_CA_significant_1sd","csp_CA_significant_2sd",
                "delta_sasa","occluded",
            ])
            for r in results_ca:
                pdb_residue_number = sequential_to_pdb_map.get(r.holo_index)
                if pdb_residue_number is None:
                    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                        _emit_warning(
                            f"[PIPE WARNING] (CA) No alignment found for sequential position {r.holo_index}",
                            log_files["tables_outputs"],
                        )
                    occlusion_data = {'delta_sasa': '', 'is_occluded': False}
                else:
                    occlusion_data = occlusion_map.get(pdb_residue_number, {'delta_sasa': '', 'is_occluded': False})

                w_ca.writerow([
                    apo_bmrb, holo_bmrb, holo_pdb, receptor_chain or "?", r.apo_index, r.apo_aa, r.holo_index, r.holo_aa,
                    f"{r.H_apo:.4f}" if r.H_apo is not None else "",
                    f"{r.N_apo:.4f}" if r.N_apo is not None else "",
                    f"{r.CA_apo:.4f}" if getattr(r, "CA_apo", None) is not None else "",
                    f"{r.H_holo:.4f}" if r.H_holo is not None else "",
                    f"{r.N_holo:.4f}" if r.N_holo is not None else "",
                    f"{r.CA_holo:.4f}" if getattr(r, "CA_holo", None) is not None else "",
                    f"{r.H_holo_original:.4f}" if r.H_holo_original is not None else "",
                    f"{r.N_holo_original:.4f}" if r.N_holo_original is not None else "",
                    f"{r.CA_holo_original:.4f}" if getattr(r, "CA_holo_original", None) is not None else "",
                    f"{r.H_offset:.4f}" if r.H_offset is not None else "",
                    f"{r.N_offset:.4f}" if r.N_offset is not None else "",
                    f"{r.CA_offset:.4f}" if getattr(r, "CA_offset", None) is not None else "",
                    f"{r.dH:.4f}" if r.dH is not None else "",
                    f"{r.dN:.4f}" if r.dN is not None else "",
                    f"{r.dCA:.4f}" if getattr(r, "dCA", None) is not None else "",
                    f"{r.csp_A:.4f}" if r.csp_A is not None else "",
                    f"{r.z_score:.4f}" if r.z_score is not None else "",
                    int(bool(r.significant)) if r.significant is not None else "",
                    int(bool(r.significant_1sd)) if r.significant_1sd is not None else "",
                    int(bool(r.significant_2sd)) if r.significant_2sd is not None else "",
                    f"{occlusion_data['delta_sasa']:.4f}" if occlusion_data['delta_sasa'] != '' else "",
                    int(bool(occlusion_data['is_occluded'])) if occlusion_data['is_occluded'] is not None else "",
                ])

    # Visualizations (H/N CSPs)
    _run_logged(
        log_files["tables_outputs"],
        plot_hsqc_variants,
        results,
        os.path.join(tgt_dir, "hsqc_scatter.png"),
        title=f"{holo_pdb} HSQC comparison",
    )
    _run_logged(
        log_files["tables_outputs"],
        write_pymol_color_csp_mask_script,
        results,
        holo_pdb,
        os.path.join(tgt_dir, "color_csp_mask.pml"),
        sasa_results=sasa_results,
        receptor_chain=receptor_chain,
        ligand_chain=ligand_chain,
        output_dir=tgt_dir_for_pymol
    )

    # Visualizations (CA-inclusive CSPs) if available
    if results_ca:
        _run_logged(
            log_files["tables_outputs"],
            plot_hsqc_variants,
            results_ca,
            os.path.join(tgt_dir, "hsqc_scatter_CA.png"),
            title=f"{holo_pdb} HSQC comparison (CA-inclusive CSPs)",
        )
    
    # Generate CSP classification bar plots and PyMOL visualizations for each threshold
    significance_thresholds = [('significant', 'original')]
    if include_alternative_thresholds:
        significance_thresholds.extend([
            ('significant_1sd', '1sd'),
            ('significant_2sd', '2sd')
        ])
    
    # Generate occlusion visualizations
    _run_logged(
        log_files["tables_outputs"],
        write_pymol_occlusion_script,
        sasa_results,
        holo_pdb,
        os.path.join(tgt_dir, "color_occlusion.pml"),
        interaction_results=interaction_results,
        ca_distance_results=ca_distance_results,
        any_atom_results=any_atom_results,
        receptor_chain=receptor_chain,
        ligand_chain=ligand_chain,
        output_dir=tgt_dir_for_pymol
    )

    # Augment interaction results with SASA and CA distance annotations for union analysis
    interaction_info = list((interaction_results or {}).get('residue_info', []))
    union_map: Dict[int, Dict[str, object]] = {}

    def _get_or_create_union_entry(
        residue_number: int,
        *,
        residue_name: Optional[str] = None,
        chain_id: Optional[str] = None,
    ) -> Dict[str, object]:
        entry = union_map.get(residue_number)
        if entry is None:
            entry = {
                "residue_number": residue_number,
                "residue_name": residue_name,
                "chain_id": chain_id,
                "has_hbond": False,
                "has_charge_complement": False,
                "has_pi_contact": False,
                "interaction_category": "none",
                "partner_residues": [],
                "hbond_count": 0,
                "charge_pair_count": 0,
                "pi_contact_count": 0,
                "hbond_partner_residues": [],
                "charge_partner_residues": [],
                "pi_partner_residues": [],
                "is_charged": False,
                "charge_type": "neutral",
                "has_sasa_occlusion": False,
                "has_ca_distance": False,
                "has_any_atom_sub_2A": False,
                "delta_sasa": None,
                "min_ca_distance": None,
                "nearest_ligand_residue": None,
                "distance_threshold": None,
            }
            union_map[residue_number] = entry
        else:
            if residue_name and not entry.get("residue_name"):
                entry["residue_name"] = residue_name
            if chain_id and not entry.get("chain_id"):
                entry["chain_id"] = chain_id
        return entry

    for info in interaction_info:
        res_num = info.get('residue_number')
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            int(res_num),
            residue_name=info.get('residue_name'),
            chain_id=info.get('chain_id', interaction_results.get('receptor_chain')),
        )
        entry["has_hbond"] = bool(info.get('has_hbond'))
        entry["has_charge_complement"] = bool(info.get('has_charge_complement'))
        entry["has_pi_contact"] = bool(info.get('has_pi_contact'))
        entry["interaction_category"] = info.get('interaction_category', entry["interaction_category"])
        entry["partner_residues"] = info.get('partner_residues', entry["partner_residues"])
        entry["hbond_count"] = info.get('hbond_count', entry["hbond_count"])
        entry["charge_pair_count"] = info.get('charge_pair_count', entry["charge_pair_count"])
        entry["pi_contact_count"] = info.get('pi_contact_count', entry["pi_contact_count"])
        entry["hbond_partner_residues"] = info.get('hbond_partner_residues', entry["hbond_partner_residues"])
        entry["charge_partner_residues"] = info.get('charge_partner_residues', entry["charge_partner_residues"])
        entry["pi_partner_residues"] = info.get('pi_partner_residues', entry["pi_partner_residues"])
        entry["is_charged"] = info.get('is_charged', entry["is_charged"])
        entry["charge_type"] = info.get('charge_type', entry["charge_type"])

    for info in sasa_results.get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            int(res_num),
            residue_name=info.get('residue_name'),
            chain_id=info.get('chain_id', sasa_results.get('receptor_chain')),
        )
        entry["has_sasa_occlusion"] = entry["has_sasa_occlusion"] or bool(info.get('is_occluded'))
        if info.get('delta_sasa') is not None:
            entry["delta_sasa"] = info.get('delta_sasa')

    for info in (ca_distance_results or {}).get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            int(res_num),
            residue_name=info.get('residue_name'),
            chain_id=info.get('chain_id', ca_distance_results.get('receptor_chain')),
        )
        entry["has_ca_distance"] = entry["has_ca_distance"] or bool(info.get('passes_filter'))
        if info.get('min_ca_distance') is not None:
            entry["min_ca_distance"] = info.get('min_ca_distance')
        if info.get('nearest_ligand_residue') is not None:
            entry["nearest_ligand_residue"] = info.get('nearest_ligand_residue')
        if info.get('distance_threshold') is not None:
            entry["distance_threshold"] = info.get('distance_threshold')

    for info in (any_atom_results or {}).get('residue_info', []):
        res_num = info.get('residue_number')
        if res_num is None:
            continue
        entry = _get_or_create_union_entry(
            int(res_num),
            residue_name=info.get('residue_name'),
            chain_id=info.get('chain_id', any_atom_results.get('receptor_chain')),
        )
        entry["has_any_atom_sub_2A"] = entry["has_any_atom_sub_2A"] or bool(info.get('passes_sub_2A_filter'))

    union_residue_info = sorted(union_map.values(), key=lambda x: x["residue_number"])
    interaction_results['residue_info'] = union_residue_info
    interaction_results['dataset_type'] = 'union'

    interaction_results['n_hbond_residues'] = sum(1 for entry in union_residue_info if entry['has_hbond'])
    interaction_results['n_charge_residues'] = sum(1 for entry in union_residue_info if entry['has_charge_complement'])
    interaction_results['n_pi_residues'] = sum(1 for entry in union_residue_info if entry['has_pi_contact'])
    interaction_results['n_interacting_residues'] = sum(
        1
        for entry in union_residue_info
        if entry['has_hbond'] or entry['has_charge_complement'] or entry['has_pi_contact']
    )
    interaction_results['n_sasa_residues'] = sum(1 for entry in union_residue_info if entry['has_sasa_occlusion'])
    interaction_results['n_ca_residues'] = sum(1 for entry in union_residue_info if entry['has_ca_distance'])
    interaction_results['n_union_residues'] = sum(
        1
        for entry in union_residue_info
        if entry['has_hbond']
        or entry['has_charge_complement']
        or entry['has_pi_contact']
        or entry['has_sasa_occlusion']
        or entry['has_ca_distance']
        or entry['has_any_atom_sub_2A']
    )
    total_union_entries = len(union_residue_info)
    interaction_results['fraction_interacting'] = (
        interaction_results['n_interacting_residues'] / total_union_entries if total_union_entries else 0.0
    )
    interaction_results['fraction_union'] = (
        interaction_results['n_union_residues'] / total_union_entries if total_union_entries else 0.0
    )


    # Generate CSP classification visualizations only for original cutoff
    significance_field = 'significant'
    threshold_suffix = 'original'
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Generating CSP classification visualizations for {threshold_suffix} threshold")
    
    # Generate bar plot (H/N)
    csp_classification_plot_path = os.path.join(tgt_dir, f"csp_classification_bars_{threshold_suffix}.png")
    _run_logged(
        log_files["tables_outputs"],
        plot_csp_classification_bars,
        results,
        interaction_results,
        csp_classification_plot_path,
        title=f"{holo_pdb} CSP Classification ({threshold_suffix})",
        significance_field=significance_field,
        include_numeric_residue_ticks=include_numeric_residue_ticks,
    )
    
    # Generate PyMOL script (H/N)
    try:
        csp_classification_script_path = os.path.join(tgt_dir, f"csp_classification_{threshold_suffix}.pml")
        _run_logged(
            log_files["tables_outputs"],
            write_pymol_csp_classification_script,
            results,
            interaction_results,
            holo_pdb,
            csp_classification_script_path,
            significance_field,
            receptor_chain=receptor_chain,
            ligand_chain=ligand_chain,
            output_dir=tgt_dir_for_pymol,
        )
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ PyMOL CSP classification script saved: {csp_classification_script_path}")
    except Exception as e:
        _emit_warning(
            f"[PIPE] ✗ Failed to generate PyMOL CSP classification script for {threshold_suffix}: {e}",
            log_files["tables_outputs"],
        )
    
    # Generate PyMOL session file (H/N)
    try:
        csp_classification_session_path = os.path.join(tgt_dir, f"csp_classification_{threshold_suffix}.pse")
        success = _run_logged(
            log_files["tables_outputs"],
            write_pymol_csp_classification_session_file,
            results,
            interaction_results,
            holo_pdb,
            csp_classification_session_path,
            significance_field,
            receptor_chain=receptor_chain,
            ligand_chain=ligand_chain,
            output_dir=tgt_dir_for_pymol,
        )
        if success and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ PyMOL CSP classification session saved: {csp_classification_session_path}")
    except Exception as e:
        _emit_warning(
            f"[PIPE] ✗ Failed to generate PyMOL CSP classification session for {threshold_suffix}: {e}",
            log_files["tables_outputs"],
        )
    
    # Keep the loop for CA-inclusive CSP classification outputs (if available)
    for significance_field, threshold_suffix in significance_thresholds:
        
        # CA-inclusive CSP classification outputs (if available)
        if results_ca:
            try:
                csp_classification_plot_ca_path = os.path.join(
                    tgt_dir, f"csp_classification_bars_{threshold_suffix}_CA.png"
                )
                _run_logged(
                    log_files["tables_outputs"],
                    plot_csp_classification_bars,
                    results_ca,
                    interaction_results,
                    csp_classification_plot_ca_path,
                    title=f"{holo_pdb} CSP Classification (CA, {threshold_suffix})",
                    significance_field=significance_field,
                    include_numeric_residue_ticks=include_numeric_residue_ticks,
                )
            except Exception as e:
                _emit_warning(
                    f"[PIPE] ✗ Failed to generate CSP classification bar plot (CA) for {threshold_suffix}: {e}",
                    log_files["tables_outputs"],
                )

            try:
                csp_classification_script_ca_path = os.path.join(
                    tgt_dir, f"csp_classification_{threshold_suffix}_CA.pml"
                )
                _run_logged(
                    log_files["tables_outputs"],
                    write_pymol_csp_classification_script,
                    results_ca,
                    interaction_results,
                    holo_pdb,
                    csp_classification_script_ca_path,
                    significance_field,
                    receptor_chain=receptor_chain,
                    ligand_chain=ligand_chain,
                    output_dir=tgt_dir_for_pymol,
                )
                if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                    print(f"[PIPE] ✓ PyMOL CSP classification script (CA) saved: {csp_classification_script_ca_path}")
            except Exception as e:
                _emit_warning(
                    f"[PIPE] ✗ Failed to generate PyMOL CSP classification script (CA) for {threshold_suffix}: {e}",
                    log_files["tables_outputs"],
                )

            try:
                csp_classification_session_ca_path = os.path.join(
                    tgt_dir, f"csp_classification_{threshold_suffix}_CA.pse"
                )
                success = _run_logged(
                    log_files["tables_outputs"],
                    write_pymol_csp_classification_session_file,
                    results_ca,
                    interaction_results,
                    holo_pdb,
                    csp_classification_session_ca_path,
                    significance_field,
                    receptor_chain=receptor_chain,
                    ligand_chain=ligand_chain,
                    output_dir=tgt_dir_for_pymol,
                )
                if success and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                    print(f"[PIPE] ✓ PyMOL CSP classification session (CA) saved: {csp_classification_session_ca_path}")
            except Exception as e:
                _emit_warning(
                    f"[PIPE] ✗ Failed to generate PyMOL CSP classification session (CA) for {threshold_suffix}: {e}",
                    log_files["tables_outputs"],
                )
    
    # Per-atom (H, N, CA) perturbation classification panels
    try:
        per_atom_panels_path = os.path.join(tgt_dir, "per_atom_classification_panels.png")
        _run_logged(
            log_files["tables_outputs"],
            plot_per_atom_classification_panels,
            results_hn=results,
            results_ca=results_ca if results_ca else None,
            binding_results=interaction_results,
            out_png=per_atom_panels_path,
            title=f"{holo_pdb} per-atom perturbations",
            threshold_config=csp_threshold_args,
            pdb_id=holo_pdb,
            structure_pdb_path=os.path.join(tgt_dir_for_pymol, f"{holo_pdb}_csp.pdb"),
            receptor_chain=receptor_chain,
        )
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ Per-atom classification panels saved: {per_atom_panels_path}")
    except Exception as e:
        _emit_warning(f"[PIPE] ✗ Failed to generate per-atom classification panels: {e}", log_files["tables_outputs"])
    
    # Generate delta SASA heatmap visualizations (same as test_sasa.py)
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Generating delta SASA heatmap visualizations")
    
    try:
        # Generate PyMOL script with delta SASA heatmap
        delta_sasa_script_path = os.path.join(tgt_dir, "delta_sasa_coloring.pml")
        _run_logged(
            log_files["tables_outputs"],
            write_pymol_delta_sasa_script,
            sasa_results,
            holo_pdb,
            delta_sasa_script_path,
            output_dir=tgt_dir_for_pymol,
        )
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ PyMOL delta SASA script saved: {delta_sasa_script_path}")
    except Exception as e:
        _emit_warning(f"[PIPE] ✗ Failed to generate PyMOL delta SASA script: {e}", log_files["tables_outputs"])
    
    try:
        # Generate PyMOL session file
        delta_sasa_session_path = os.path.join(tgt_dir, "delta_sasa_coloring.pse")
        success = _run_logged(
            log_files["tables_outputs"],
            write_pymol_session_file,
            sasa_results,
            holo_pdb,
            delta_sasa_session_path,
            output_dir=tgt_dir_for_pymol,
        )
        if success and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ PyMOL delta SASA session saved: {delta_sasa_session_path}")
    except Exception as e:
        _emit_warning(f"[PIPE] ✗ Failed to generate PyMOL delta SASA session: {e}", log_files["tables_outputs"])
    
    if binary_mode:
        binary_dir = os.path.join(tgt_dir, "binary_pymol")
        os.makedirs(binary_dir, exist_ok=True)
        # Normalize binary_dir for PyMOL scripts (relative to project root)
        binary_dir_for_pymol = binary_dir
        if os.path.isabs(binary_dir):
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            binary_dir_for_pymol = os.path.relpath(binary_dir, project_root)
        binary_dir_for_pymol = binary_dir_for_pymol.replace('\\', '/')
        if not binary_dir_for_pymol.startswith('./') and not binary_dir_for_pymol.startswith('/'):
            binary_dir_for_pymol = './' + binary_dir_for_pymol
        if not binary_dir_for_pymol.endswith('/'):
            binary_dir_for_pymol += '/'
        
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] Generating binary PyMOL visualizations in {binary_dir}")

        try:
            delta_sasa_binary_path = os.path.join(binary_dir, "delta_sasa_coloring_binary.pml")
            _run_logged(
                log_files["tables_outputs"],
                write_pymol_delta_sasa_script,
                sasa_results,
                holo_pdb,
                delta_sasa_binary_path,
                binary_mode=True,
                output_dir=binary_dir_for_pymol,
            )
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ PyMOL binary delta SASA script saved: {delta_sasa_binary_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate binary delta SASA script: {e}", log_files["tables_outputs"])

        try:
            delta_sasa_binary_session = os.path.join(binary_dir, "delta_sasa_coloring_binary.pse")
            success = _run_logged(
                log_files["tables_outputs"],
                write_pymol_session_file,
                sasa_results,
                holo_pdb,
                delta_sasa_binary_session,
                binary_mode=True,
                output_dir=binary_dir_for_pymol,
            )
            if success and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ PyMOL binary delta SASA session saved: {delta_sasa_binary_session}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate binary delta SASA session: {e}", log_files["tables_outputs"])

        try:
            csp_heatmap_binary_path = os.path.join(binary_dir, "csp_heatmap_binary.pml")
            _run_logged(
                log_files["tables_outputs"],
                write_pymol_csp_heatmap_script,
                results,
                sasa_results,
                holo_pdb,
                csp_heatmap_binary_path,
                binary_mode=True,
                output_dir=binary_dir_for_pymol,
            )
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ PyMOL binary CSP heatmap script saved: {csp_heatmap_binary_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate binary CSP heatmap script: {e}", log_files["tables_outputs"])

        try:
            csp_binary_session_path = os.path.join(binary_dir, "csp_heatmap_binary.pse")
            success = _run_logged(
                log_files["tables_outputs"],
                write_pymol_csp_session_file,
                results,
                sasa_results,
                holo_pdb,
                csp_binary_session_path,
                binary_mode=True,
            )
            if success and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ PyMOL binary CSP heatmap session saved: {csp_binary_session_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate binary CSP heatmap session: {e}", log_files["tables_outputs"])

        try:
            combined_binary_path = os.path.join(binary_dir, "color_combined_binary.pml")
            _run_logged(
                log_files["tables_outputs"],
                write_pymol_combined_script,
                results,
                interaction_results,
                holo_pdb,
                combined_binary_path,
                binary_mode=True,
                receptor_chain=receptor_chain,
                ligand_chain=ligand_chain,
                sasa_results=sasa_results,
                output_dir=binary_dir_for_pymol,
            )
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ PyMOL binary combined script saved: {combined_binary_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate binary combined PyMOL script: {e}", log_files["tables_outputs"])
    
    # Generate 1D single-atom shift analysis for this target
    try:
        _run_logged(log_files["tables_outputs"], compute_1d_metrics_for_target, Path(tgt_dir))
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ 1D single-atom analysis saved: {os.path.join(tgt_dir, '1d_analysis.csv')}")
    except Exception as e:
        _emit_warning(f"[PIPE] ✗ Failed to generate 1D single-atom analysis: {e}", log_files["tables_outputs"])

    _console_line(f"[PIPE] [{target_label}] Step 5/5 Master CSV + case study")
    # Generate master CSV file with alignment
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] Generating master CSV with sequence alignment")
    
    try:
        master_csv_path = os.path.join(tgt_dir, "master_alignment.csv")
        _run_logged(log_files["master_csv"], merge_all_csv_files, tgt_dir, master_csv_path)
        if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] ✓ Master CSV saved: {master_csv_path}")
    except Exception as e:
        _emit_warning(f"[PIPE] ✗ Failed to generate master CSV: {e}", log_files["master_csv"])

    # Generate case-study figures (v1 and v2) for each target
    if generate_case_study:
        # If forcing view reset, capture once in case_study v1 and let v2 reuse it.
        force_reset_for_case_study_1 = force_case_study_view_reset
        force_reset_for_case_study_2 = False if force_case_study_view_reset else force_case_study_view_reset
        try:
            case_study_path = _run_logged(
                log_files["case_study"],
                generate_case_study_figure,
                target_dir=tgt_dir,
                pdb_id=holo_pdb,
                apo_bmrb=apo_bmrb,
                holo_bmrb=holo_bmrb,
                force_view_reset=force_reset_for_case_study_1,
                view_key=target_label,
            )
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ Case-study figure saved: {case_study_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate case-study figure: {e}", log_files["case_study"])
        try:
            case_study_2_path = _run_logged(
                log_files["case_study"],
                generate_case_study_2_figure,
                target_dir=tgt_dir,
                pdb_id=holo_pdb,
                apo_bmrb=apo_bmrb,
                holo_bmrb=holo_bmrb,
                force_view_reset=force_reset_for_case_study_2,
                view_key=target_label,
            )
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] ✓ Case-study v2 figure saved: {case_study_2_path}")
        except Exception as e:
            _emit_warning(f"[PIPE] ✗ Failed to generate case-study v2 figure: {e}", log_files["case_study"])
    
    _console_line(f"[PIPE] [{target_label}] Complete -> {tgt_dir}")


def lookup_all_rows_from_holo_pdb(csv_path: str, holo_pdb: str) -> List[Dict[str, str]]:
    """Look up all rows with matching holo_pdb from CSP_UBQ.csv."""
    rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("holo_pdb", "").strip() == holo_pdb:
                rows.append({
                    "apo_bmrb": row.get("apo_bmrb", "").strip(),
                    "holo_bmrb": row.get("holo_bmrb", "").strip(),
                    "holo_pdb": row.get("holo_pdb", "").strip()
                })
    if not rows:
        raise ValueError(f"No entry found for holo_pdb ID: {holo_pdb}")
    return rows


def lookup_bmrb_ids_from_holo_pdb(csv_path: str, holo_pdb: str) -> Dict[str, str]:
    """Look up apo_bmrb and holo_bmrb from CSP_UBQ.csv using holo_pdb ID.
    
    DEPRECATED: Use lookup_all_rows_from_holo_pdb for duplicate handling.
    This function returns only the first match for backward compatibility.
    """
    rows = lookup_all_rows_from_holo_pdb(csv_path, holo_pdb)
    return rows[0]


def get_duplicate_index_for_row(csv_path: str, target_row: Dict[str, str]) -> Optional[int]:
    """Determine the index (1-based) of a specific row among duplicates with the same holo_pdb.
    
    Returns None if there are no duplicates, or the 1-based index if duplicates exist.
    """
    holo_pdb = (target_row.get("holo_pdb") or "").strip()
    if not holo_pdb:
        return None
    
    apo_bmrb = (target_row.get("apo_bmrb") or "").strip()
    holo_bmrb = (target_row.get("holo_bmrb") or "").strip()
    
    # Find all rows with the same holo_pdb
    matching_rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("holo_pdb", "").strip() == holo_pdb:
                matching_rows.append({
                    "apo_bmrb": row.get("apo_bmrb", "").strip(),
                    "holo_bmrb": row.get("holo_bmrb", "").strip(),
                    "holo_pdb": row.get("holo_pdb", "").strip()
                })
    
    # If only one match, no duplicates
    if len(matching_rows) <= 1:
        return None
    
    # Find the index of the current row among duplicates
    for idx, row in enumerate(matching_rows, start=1):
        if (row.get("apo_bmrb") == apo_bmrb and 
            row.get("holo_bmrb") == holo_bmrb):
            return idx
    
    # If exact match not found, return None (shouldn't happen in normal usage)
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="CSP Pipeline v1")
    parser.add_argument("--input", default=paths.input_csv)
    parser.add_argument("--out", default=paths.outputs_dir)
    
    # Create mutually exclusive group for filtering
    filter_group = parser.add_mutually_exclusive_group()
    filter_group.add_argument("--ids", help="Comma-separated apo_bmrb IDs to process (filters rows)")
    filter_group.add_argument("--holo-pdb", help="Single holo_pdb ID to process (looks up corresponding apo_bmrb and holo_bmrb)")
    
    parser.add_argument("--workers", type=int, default=concurrency.workers)
    
    # SASA analysis parameters (same as test_sasa.py)
    parser.add_argument("--probe-radius", type=float, default=0.0, 
                       help="Probe radius for SASA calculation (default: 0.0)")
    parser.add_argument("--nh-mode", choices=['sum', 'n_only'], default='sum',
                       help="NH calculation mode (default: sum)")
    parser.add_argument("--mode", choices=['standard', 'fast', 'residue'], default='residue',
                       help="SASA computation mode (default: residue)")
    parser.add_argument("--threshold", type=float, default=sasa_analysis.sasa_threshold,
                       help="SASA threshold for occlusion (default: from config)")
    
    # Interaction analysis parameters
    parser.add_argument("--interaction-distance", type=float, default=4.5,
                       help="Distance threshold for charge complementarity (default: 4.5 Å)")
    
    # CSP threshold arguments
    parser.add_argument("--outlier-z", type=float, 
                       help="Z-score for outlier detection (default: 3.0)")
    parser.add_argument("--significance-z", type=float,
                       help="Z-score for final threshold (default: 0.0)")
    parser.add_argument("--max-outlier-iterations", type=int,
                       help="Max iterations for outlier removal (default: 10)")
    parser.add_argument("--max-outlier-fraction", type=float,
                       help="Max fraction of residues to remove (default: 0.2)")
    parser.add_argument("--absolute-cutoff", type=float,
                       help="Absolute CSP cutoff (overrides all other methods)")
    parser.add_argument("--binary-visualizations", action="store_true",
                       help="Also emit binary red/blue PyMOL visualizations into a separate subdirectory.")
    parser.add_argument("--include-alternative-thresholds", action="store_true",
                       help="Include 1sd and 2sd significance cutoff analyses (default: False, only original cutoff)")
    parser.add_argument(
        "--run-csv-metadata-annotation",
        action="store_true",
        help="Run ec_classes and scope_fold_type metadata annotation on input CSV before processing (opt-in).",
    )
    parser.add_argument(
        "--no-case-study",
        action="store_true",
        help="Skip per-target case-study figure generation (figures are generated by default).",
    )
    parser.add_argument(
        "--force-case-study-view-reset",
        action="store_true",
        help="Always recapture PyMOL perspective (F5) even if a saved view already exists.",
    )
    parser.add_argument(
        "--include-numeric-residue-ticks",
        action="store_true",
        help="Include numeric residue index rows under CSP classification bar plots (default: off, amino-acid row only).",
    )
    
    args = parser.parse_args()

    generate_case_study = not args.no_case_study
    if generate_case_study and not args.holo_pdb and args.workers and args.workers > 1:
        print("[PIPE] WARNING: Case-study in batch mode may use interactive PyMOL for targets without saved views; forcing --workers=1.")
        args.workers = 1

    should_annotate_csv = args.run_csv_metadata_annotation and (args.holo_pdb is None)
    if should_annotate_csv:
        try:
            annotate_csv_with_ec_and_scope(
                Path(args.input),
                verbose=(os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")),
            )
        except Exception as exc:
            print(f"[PIPE] WARNING: Failed to annotate CSV metadata columns: {exc}")
    elif args.run_csv_metadata_annotation and args.holo_pdb is not None:
        print("[PIPE] WARNING: --run-csv-metadata-annotation is ignored when --holo-pdb is set.")

    ensure_directories(args.out)

    rows: List[Dict[str, str]] = []
    row_suffixes: Dict[int, Optional[str]] = {}  # Map row index to directory suffix
    
    if args.holo_pdb:
        # Look up all entries with matching holo_pdb from CSP_UBQ.csv
        try:
            matching_rows = lookup_all_rows_from_holo_pdb(args.input, args.holo_pdb)
            
            # Check if there are duplicates
            if len(matching_rows) > 1:
                # Multiple entries found - create unique subdirectories for ALL entries
                # IMPORTANT: When duplicates exist, we MUST use suffixes for all entries
                print(f"[PIPE] Found {len(matching_rows)} entries for holo_pdb {args.holo_pdb}, creating unique subdirectories")
                for idx, row in enumerate(matching_rows, start=1):
                    rows.append(row)
                    row_suffixes[len(rows) - 1] = str(idx)  # Always set suffix when duplicates exist
                    print(f"[PIPE] Entry {idx}/{len(matching_rows)}: apo_bmrb={row['apo_bmrb']}, holo_bmrb={row['holo_bmrb']} -> {args.holo_pdb}_{idx}")
            else:
                # Single entry - use current behavior (no suffix)
                rows.append(matching_rows[0])
                row_suffixes[0] = None
                if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                    print(f"[PIPE] Found single entry for holo_pdb {args.holo_pdb}: apo_bmrb={matching_rows[0]['apo_bmrb']}, holo_bmrb={matching_rows[0]['holo_bmrb']}")
        except ValueError as e:
            print(f"[PIPE] ERROR: {e}")
            return
    else:
        # Original logic for --ids filtering
        filter_ids = None
        if args.ids:
            filter_ids = {s.strip() for s in args.ids.split(",") if s.strip()}

        # First pass: collect all rows and build a mapping of holo_pdb to all matching rows
        all_rows = []
        holo_pdb_to_rows: Dict[str, List[Dict[str, str]]] = {}
        
        with open(args.input, "r", newline="") as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                if filter_ids and row.get("apo_bmrb") not in filter_ids:
                    continue
                all_rows.append(row)
                holo_pdb = (row.get("holo_pdb") or "").strip()
                if holo_pdb:
                    if holo_pdb not in holo_pdb_to_rows:
                        holo_pdb_to_rows[holo_pdb] = []
                    holo_pdb_to_rows[holo_pdb].append(row)
        
        # Second pass: determine directory suffixes for rows with duplicates
        # For each row, find its position among duplicates with the same holo_pdb
        for row in all_rows:
            holo_pdb = (row.get("holo_pdb") or "").strip()
            apo_bmrb = (row.get("apo_bmrb") or "").strip()
            holo_bmrb = (row.get("holo_bmrb") or "").strip()
            
            # Check if this holo_pdb has duplicates
            matching_rows = holo_pdb_to_rows.get(holo_pdb, [])
            if len(matching_rows) > 1:
                # Has duplicates - find which duplicate this row is
                # Match based on apo_bmrb and holo_bmrb to ensure correct identification
                duplicate_index = None
                for idx, match_row in enumerate(matching_rows, start=1):
                    match_apo = (match_row.get("apo_bmrb") or "").strip()
                    match_holo = (match_row.get("holo_bmrb") or "").strip()
                    if match_apo == apo_bmrb and match_holo == holo_bmrb:
                        duplicate_index = idx
                        break
                
                if duplicate_index:
                    rows.append(row)
                    row_suffixes[len(rows) - 1] = str(duplicate_index)
                else:
                    # Shouldn't happen, but fallback to position-based matching
                    rows.append(row)
                    # Count how many rows with same holo_pdb we've already processed
                    seen_count = sum(1 for i in range(len(rows) - 1) 
                                   if (rows[i].get("holo_pdb") or "").strip() == holo_pdb)
                    row_suffixes[len(rows) - 1] = str(seen_count + 1)
            else:
                # No duplicates - use current behavior (no suffix)
                rows.append(row)
                row_suffixes[len(rows) - 1] = None

    # Prepare SASA arguments
    sasa_args = {
        'probe_radius': args.probe_radius,
        'nh_mode': args.nh_mode,
        'mode': args.mode,
        'threshold': args.threshold
    }
    
    # Prepare interaction arguments
    interaction_args = {
        'distance_threshold': args.interaction_distance
    }
    
    # Prepare CSP threshold arguments
    csp_threshold_args = {}
    if args.outlier_z is not None:
        csp_threshold_args['outlier_z_score'] = args.outlier_z
    if args.significance_z is not None:
        csp_threshold_args['significance_z_score'] = args.significance_z
    if args.max_outlier_iterations is not None:
        csp_threshold_args['max_outlier_iterations'] = args.max_outlier_iterations
    if args.max_outlier_fraction is not None:
        csp_threshold_args['max_outlier_fraction'] = args.max_outlier_fraction
    if args.absolute_cutoff is not None:
        csp_threshold_args['absolute_cutoff'] = args.absolute_cutoff
    
    if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(f"[PIPE] SASA parameters: probe_radius={sasa_args['probe_radius']}, nh_mode={sasa_args['nh_mode']}, mode={sasa_args['mode']}, threshold={sasa_args['threshold']}")
        if csp_threshold_args:
            print(f"[PIPE] CSP threshold overrides: {csp_threshold_args}")
        else:
            print(f"[PIPE] Using default CSP thresholds")
    
    # Matplotlib is not thread-safe. Use parallel only for single-target mode (--holo-pdb) with no duplicates.
    # Batch mode (many targets) and duplicate holo_pdb entries must run serially.
    has_duplicates = args.holo_pdb and len(rows) > 1
    use_parallel = args.workers and args.workers > 1 and args.holo_pdb and not has_duplicates

    if use_parallel:
        # Parallel processing for non-duplicate cases
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] Scheduling {len(rows)} rows across {args.workers} workers")
            futures = [
                ex.submit(
                    process_row,
                    row,
                    args.out,
                    sasa_args,
                    csp_threshold_args,
                    interaction_args,
                    args.binary_visualizations,
                    args.include_alternative_thresholds,
                    row_suffixes.get(idx, None),
                    generate_case_study,
                    args.include_numeric_residue_ticks,
                    args.force_case_study_view_reset,
                )
                for idx, row in enumerate(rows)
            ]
            for f in as_completed(futures):
                f.result()  # Surface any exception from worker
    else:
        # Serial processing (batch mode, duplicates, or workers=1)
        if has_duplicates:
            if (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
                print(f"[PIPE] Processing {len(rows)} duplicate rows serially to avoid matplotlib threading issues")
        elif not args.holo_pdb and (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] Processing {len(rows)} rows serially (batch mode; matplotlib is not thread-safe)")
        elif (os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
            print(f"[PIPE] Processing {len(rows)} rows serially")
        for idx, row in enumerate(rows):
            process_row(
                row,
                args.out,
                sasa_args,
                csp_threshold_args,
                interaction_args,
                args.binary_visualizations,
                args.include_alternative_thresholds,
                row_suffixes.get(idx, None),
                generate_case_study,
                args.include_numeric_residue_ticks,
                args.force_case_study_view_reset,
            )

    # Generate confusion_matrix_per_system.csv for downstream scripts (create_si_fig_s20, etc.)
    _console_line("[PIPE] Finalize: refresh confusion-matrix summary")
    run_logs_dir = os.path.join(args.out, "logs")
    os.makedirs(run_logs_dir, exist_ok=True)
    run_summary_log = os.path.join(run_logs_dir, "confusion_matrix_summary.txt")
    _run_logged(
        run_summary_log,
        generate_confusion_matrix_per_system,
        args.out,
        verbose=(os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")),
    )


if __name__ == "__main__":
    main()


