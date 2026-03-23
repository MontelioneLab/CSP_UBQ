"""
Utilities to download and parse BMRB chemical shift files (NMR-STAR).

Outputs per entry:
- amino acid sequence string as reported alongside chemical shift assignments
- per-residue dictionaries for backbone amide H and N chemical shifts

Notes:
- We accept H atom names of {"H", "HN", "H1"} and N as {"N"}
- Residue indexing is returned as 1-based sequence position (not BMRB seq_id)
- Multiple atoms per residue are aggregated by median to reduce outliers
"""

from __future__ import annotations

import io
import os
import re
import statistics
from typing import Dict, Tuple, Optional, List

import requests

# Verbose printing controlled by env var CSP_VERBOSE
def _vprint(*args, **kwargs):
    import os as _os
    if (_os.environ.get("CSP_VERBOSE", "").lower() in ("1", "true", "yes")):
        print(*args, **kwargs)

try:
    from .config import network, paths
except Exception:
    import os as _os, sys as _sys
    _sys.path.append(_os.path.dirname(_os.path.dirname(__file__)))
    from scripts.config import network, paths


class DownloadError(RuntimeError):
    pass


def _download_bmrb_nmr_star(entry_id: str) -> Tuple[bytes, str]:
    """Download BMRB NMR-STAR file.

    Returns (content, version_str), where version_str is '21' or '3'.
    Tries version 21 first, then falls back to version 3.
    """
    # Try primary URL (version 21) first
    primary_url = network.bmrb_ftp_url_template.format(id=entry_id)
    fallback_url = network.bmrb_ftp_fallback_url_template.format(id=entry_id)
    
    last_err: Optional[Exception] = None
    
    # Try primary URL with retries
    for attempt in range(network.retries):
        try:
            _vprint(f"[BMRB] GET {primary_url}")
            resp = requests.get(primary_url, timeout=(network.connect_timeout, network.read_timeout))
            if resp.status_code == 200 and resp.content:
                _vprint(f"[BMRB] OK {entry_id} ({len(resp.content)} bytes)")
                return resp.content, "21"
            last_err = DownloadError(f"HTTP {resp.status_code} for {primary_url}")
        except Exception as e:  # network errors
            last_err = e
        # simple backoff without sleep (no async); retries will try templates again
    
    # If primary URL failed, try fallback URL
    _vprint(f"[BMRB] Primary URL failed, trying fallback: {fallback_url}")
    for attempt in range(network.retries):
        try:
            _vprint(f"[BMRB] GET {fallback_url}")
            resp = requests.get(fallback_url, timeout=(network.connect_timeout, network.read_timeout))
            if resp.status_code == 200 and resp.content:
                _vprint(f"[BMRB] OK {entry_id} (fallback, {len(resp.content)} bytes)")
                return resp.content, "3"
            last_err = DownloadError(f"HTTP {resp.status_code} for {fallback_url}")
        except Exception as e:  # network errors
            last_err = e
    
    # Both URLs failed
    raise DownloadError(f"Failed to download BMRB {entry_id}: Neither version 21 nor version 3 available. Last error: {last_err}")


def fetch_bmrb(entry_id: str, cache_dir: Optional[str] = None, force: bool = False) -> str:
    """Fetch an NMR-STAR file and cache it locally.

    Returns the absolute path to the cached file.
    """
    cache_dir = cache_dir or paths.cs_cache_dir
    os.makedirs(cache_dir, exist_ok=True)
    
    # Try primary URL (version 21) first
    primary_url = network.bmrb_ftp_url_template.format(id=entry_id)
    fallback_url = network.bmrb_ftp_fallback_url_template.format(id=entry_id)
    
    # Determine which URL to use and what suffix to preserve
    primary_path = os.path.join(cache_dir, f"{entry_id}_21.str")
    fallback_path = os.path.join(cache_dir, f"{entry_id}_3.str")
    
    # Check if we already have a cached file
    if not force:
        if os.path.exists(primary_path) and os.path.getsize(primary_path) > 0:
            _vprint(f"[BMRB] Cache hit {primary_path}")
            return primary_path
        elif os.path.exists(fallback_path) and os.path.getsize(fallback_path) > 0:
            _vprint(f"[BMRB] Cache hit {fallback_path}")
            return fallback_path
    
    # Try download with internal fallback and save using the detected version suffix
    try:
        content, version = _download_bmrb_nmr_star(entry_id)
        if version == "21":
            out_path = primary_path
        else:
            out_path = fallback_path
        with open(out_path, "wb") as f:
            f.write(content)
        _vprint(f"[BMRB] Wrote {out_path}")
        return out_path
    except DownloadError as e:
        raise DownloadError(f"Failed to download BMRB {entry_id}: Neither version 21 nor version 3 available. Last error: {e}")


def _detect_bmrb_format(star_path: str) -> str:
    """Detect whether the STR file is version 21 or version 3 format.
    
    Returns '21' or '3' based on the filename suffix or content analysis.
    """
    filename = os.path.basename(star_path)
    if filename.endswith('_21.str'):
        return '21'
    elif filename.endswith('_3.str'):
        return '3'
    else:
        # Analyze content to determine format
        with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
            if "_Atom_chem_shift.Val" in content:
                return '3'
            elif "_Chem_shift_value" in content:
                return '21'
            else:
                # Default to 21 if unclear
                return '21'


def _parse_chem_shifts_v3(lines: List[str]) -> List[Tuple[str, List[Dict[str, str]]]]:
    """Parse chemical shift saveframes for version 3 format.
    
    Version 3 uses _Atom_chem_shift.Val instead of _Chem_shift_value.
    """
    saveframes: List[Tuple[str, List[Dict[str, str]]]] = []
    i = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Look for any saveframe
        if line.startswith("save_"):
            saveframe_name = line
            saveframe_start = i
            
            # Check if this saveframe contains chemical shift data
            # First check the saveframe name pattern
            name_matches = (
                "assigned_chemical_shifts" in line or 
                "_Assigned_chem_shift_list" in line or
                "shift_set" in line or
                "chem_shift" in line
            )
            
            # Also check the _Assigned_chem_shift_list.Sf_category tag within the saveframe
            category_matches = False
            j = i + 1
            while j < len(lines) and lines[j].strip() != "save_":
                cat_line = lines[j].strip()
                if cat_line.startswith("_Assigned_chem_shift_list.Sf_category") or cat_line.startswith("_Saveframe_category"):
                    # Extract the category value
                    if "assigned_chemical_shifts" in cat_line or "assigned_chem_shift" in cat_line.lower():
                        category_matches = True
                    break
                j += 1
            
            if name_matches or category_matches:
                _vprint(f"[BMRB] Found chemical shift saveframe (v3): {saveframe_name}")
                
                # Parse this saveframe
                rows = _parse_single_chem_shift_saveframe_v3(lines, saveframe_start)
                if rows:
                    saveframes.append((saveframe_name, rows))
                    _vprint(f"[BMRB] Parsed {len(rows)} rows from {saveframe_name}")
            
            # Skip to end of this saveframe
            while i < len(lines) and lines[i].strip() != "save_":
                i += 1
            i += 1
        else:
            i += 1
    
    return saveframes


def _parse_single_chem_shift_saveframe_v3(lines: List[str], start_idx: int) -> List[Dict[str, str]]:
    """Parse a single chemical shift saveframe for version 3 format."""
    rows: List[Dict[str, str]] = []
    i = start_idx + 1
    in_saveframe = True
    
    while i < len(lines) and in_saveframe:
        line = lines[i].strip()
        
        # End of saveframe
        if line == "save_":
            in_saveframe = False
            break
            
        # Look for loop_ within the saveframe
        if line == "loop_":
            i += 1
            tags: List[str] = []
            
            # Collect tags
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("_"):
                    tags.append(line)
                    i += 1
                else:
                    break
                    
            if not tags:
                continue
                
            # Check if this is the chemical shift data loop (v3 format)
            chem_shift_tags = any(
                "_Atom_chem_shift.Val" in tag or 
                "_Atom_chem_shift.Seq_ID" in tag or
                "_Atom_chem_shift.Atom_ID" in tag or
                "_Atom_chem_shift.Comp_ID" in tag
                for tag in tags
            )
            if not chem_shift_tags:
                _vprint(f"[BMRB] Skipping non-chemical shift loop with tags: {tags[:3]}...")
                # Skip this loop, continue to next
                while i < len(lines) and lines[i].strip() != "stop_":
                    i += 1
                i += 1
                continue
                
            _vprint(f"[BMRB] Found chemical shift loop (v3) with {len(tags)} tags")
                
            # Collect data rows until stop_ or end of saveframe
            while i < len(lines):
                line = lines[i].strip()
                if line == "stop_" or line == "save_":
                    break
                if not line:
                    i += 1
                    continue
                    
                # Parse data line
                tokens = _split_star_line(line)
                if len(tokens) == len(tags):
                    row = {tags[j]: tokens[j] for j in range(len(tags))}
                    # Preserve tag order to allow positional value inference later
                    row["__tags_order__"] = "\t".join(tags)
                    rows.append(row)
                elif len(tokens) > 0:  # Handle multi-line data
                    # Accumulate tokens until we have enough
                    buf = tokens[:]
                    i += 1
                    while i < len(lines) and len(buf) < len(tags):
                        next_line = lines[i].strip()
                        if next_line == "stop_" or next_line == "save_":
                            break
                        if next_line:
                            next_tokens = _split_star_line(next_line)
                            buf.extend(next_tokens)
                        i += 1
                    if len(buf) >= len(tags):
                        row = {tags[j]: buf[j] for j in range(len(tags))}
                        row["__tags_order__"] = "\t".join(tags)
                        rows.append(row)
                i += 1
        else:
            i += 1
            
    return rows


_loop_start_re = re.compile(r"^\s*loop_\s*$", re.IGNORECASE)
_tag_re = re.compile(r"^\s*_([\w\.:]+)\s*$")
_data_stop_re = re.compile(r"^\s*stop_\s*$", re.IGNORECASE)


def _tokenize_star_lines(text: str) -> List[str]:
    # Keep simple: split into lines and strip trailing newlines; do not strip leading spaces
    return text.splitlines()


def _parse_all_chem_shift_saveframes(lines: List[str]) -> List[Tuple[str, List[Dict[str, str]]]]:
    """Parse all chemical shift saveframes and return list of (saveframe_name, rows) tuples.
    
    Handles multiple chemical shift saveframes in a single file.
    """
    saveframes: List[Tuple[str, List[Dict[str, str]]]] = []
    i = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Look for any saveframe
        if line.startswith("save_"):
            saveframe_name = line
            saveframe_start = i
            
            # Check if this saveframe contains chemical shift data
            # First check the saveframe name pattern
            name_matches = (
                "assigned_chem_shift_list" in line or 
                "_Assigned_chem_shift_list" in line or
                "shift_set" in line or
                "_shift_set" in line or
                "chem_shift_list" in line  # Added for save_chem_shift_list_1 pattern
            )
            
            # Also check the _Saveframe_category tag within the saveframe
            category_matches = False
            j = i + 1
            while j < len(lines) and lines[j].strip() != "save_":
                cat_line = lines[j].strip()
                if cat_line.startswith("_Saveframe_category"):
                    # Extract the category value
                    if "assigned_chemical_shifts" in cat_line or "assigned_chem_shift" in cat_line.lower():
                        category_matches = True
                    break
                j += 1
            
            if name_matches or category_matches:
                _vprint(f"[BMRB] Found chemical shift saveframe: {saveframe_name}")
                
                # Parse this saveframe
                rows = _parse_single_chem_shift_saveframe(lines, saveframe_start)
                if rows:
                    saveframes.append((saveframe_name, rows))
                    _vprint(f"[BMRB] Parsed {len(rows)} rows from {saveframe_name}")
            
            # Skip to end of this saveframe
            while i < len(lines) and lines[i].strip() != "save_":
                i += 1
            i += 1
        else:
            i += 1
            
    return saveframes


def _parse_single_chem_shift_saveframe(lines: List[str], start_idx: int) -> List[Dict[str, str]]:
    """Parse a single chemical shift saveframe starting at start_idx."""
    rows: List[Dict[str, str]] = []
    i = start_idx + 1
    in_saveframe = True
    
    while i < len(lines) and in_saveframe:
        line = lines[i].strip()
        
        # End of saveframe
        if line == "save_":
            in_saveframe = False
            break
            
        # Look for loop_ within the saveframe
        if line == "loop_":
            i += 1
            tags: List[str] = []
            
            # Collect tags
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("_"):
                    tags.append(line)
                    i += 1
                else:
                    break
                    
            if not tags:
                continue
                
            # Check if this is the chemical shift data loop
            chem_shift_tags = any(
                "Atom_shift_assign_ID" in tag or 
                tag == "_Chem_shift_value" or
                "Residue_seq_code" in tag or
                "Atom_name" in tag
                for tag in tags
            )
            if not chem_shift_tags:
                _vprint(f"[BMRB] Skipping non-chemical shift loop with tags: {tags[:3]}...")
                # Skip this loop, continue to next
                while i < len(lines) and lines[i].strip() != "stop_":
                    i += 1
                i += 1
                continue
                
            _vprint(f"[BMRB] Found chemical shift loop with {len(tags)} tags")
                
            # Collect data rows until stop_ or end of saveframe
            while i < len(lines):
                line = lines[i].strip()
                if line == "stop_" or line == "save_":
                    break
                if not line:
                    i += 1
                    continue
                    
                # Parse data line
                tokens = _split_star_line(line)
                if len(tokens) == len(tags):
                    row = {tags[j]: tokens[j] for j in range(len(tags))}
                    # Preserve tag order to allow positional value inference later
                    row["__tags_order__"] = "\t".join(tags)
                    rows.append(row)
                elif len(tokens) > 0:  # Handle multi-line data
                    # Accumulate tokens until we have enough
                    buf = tokens[:]
                    i += 1
                    while i < len(lines) and len(buf) < len(tags):
                        next_line = lines[i].strip()
                        if next_line == "stop_" or next_line == "save_":
                            break
                        if next_line:
                            next_tokens = _split_star_line(next_line)
                            buf.extend(next_tokens)
                        i += 1
                    if len(buf) >= len(tags):
                        row = {tags[j]: buf[j] for j in range(len(tags))}
                        row["__tags_order__"] = "\t".join(tags)
                        rows.append(row)
                i += 1
        else:
            i += 1
            
    return rows


def _split_star_line(line: str) -> List[str]:
    tokens: List[str] = []
    buf = []
    in_single = False
    in_double = False
    i = 0
    while i < len(line):
        c = line[i]
        if c == "'" and not in_double:
            in_single = not in_single
            i += 1
            continue
        if c == '"' and not in_single:
            in_double = not in_double
            i += 1
            continue
        if not in_single and not in_double and c.isspace():
            if buf:
                tokens.append("".join(buf))
                buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    if buf:
        tokens.append("".join(buf))
    return tokens


def _export_chem_shifts_to_csv(star_path: str, all_rows: List[Dict[str, str]], saveframes: List[Tuple[str, List[Dict[str, str]]]]) -> None:
    """Export parsed chemical shift data to CSV files."""
    import csv
    
    # Create output directory
    base_name = os.path.splitext(os.path.basename(star_path))[0]
    output_dir = os.path.join(os.path.dirname(star_path), "parsed")
    os.makedirs(output_dir, exist_ok=True)
    
    # Export combined data
    combined_csv_path = os.path.join(output_dir, f"{base_name}_combined.csv")
    if all_rows:
        # Get all unique field names
        all_fields = set()
        for row in all_rows:
            all_fields.update(row.keys())
        fieldnames = sorted(all_fields)
        
        with open(combined_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        _vprint(f"[BMRB] Exported combined data to {combined_csv_path}")
    
    # Export individual saveframes
    for i, (saveframe_name, rows) in enumerate(saveframes):
        if not rows:
            continue
            
        # Clean saveframe name for filename
        clean_name = saveframe_name.replace("save_", "").replace("_", "-")
        saveframe_csv_path = os.path.join(output_dir, f"{base_name}_{clean_name}.csv")
        
        # Get field names for this saveframe
        fieldnames = sorted(set().union(*(row.keys() for row in rows)))
        
        with open(saveframe_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        _vprint(f"[BMRB] Exported {saveframe_name} to {saveframe_csv_path}")


def _parse_sequence_and_shifts_v3(star_path: str) -> Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float]]:
    """Extract sequence and per-residue H/N/CA chemical shifts from NMR-STAR version 3.
    
    Version 3 format uses _Atom_chem_shift.Val instead of _Chem_shift_value.
    """
    with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    lines = _tokenize_star_lines(text)
    
    # Parse all chemical shift saveframes using v3 parser
    saveframes = _parse_chem_shifts_v3(lines)
    _vprint(f"[BMRB] Found {len(saveframes)} chemical shift saveframes (v3) in {os.path.basename(star_path)}")
    
    # Combine data from all saveframes
    all_rows = []
    for saveframe_name, rows in saveframes:
        all_rows.extend(rows)
        _vprint(f"[BMRB] Parsed {len(rows)} rows from {saveframe_name}")
    
    _vprint(f"[BMRB] Total parsed {len(all_rows)} chem shift rows (v3) from {os.path.basename(star_path)}")

    # Extract sequence from monomeric polymer saveframe first
    sequence = _extract_sequence_from_saveframe(lines)
    if not sequence:
        _vprint(f"[BMRB] Warning: Could not extract sequence from saveframe, using CS data")
        sequence = ""

    # Collect per residue atom shifts from chemical shift data (v3 format)
    residue_to_atoms: Dict[Tuple[int, str], Dict[str, List[float]]] = {}
    for r in all_rows:
        try:
            # Extract key fields - handle the v3 BMRB tag names
            seq_id = None
            comp_id = ""
            atom_id = ""
            value = None
            
            # Try different tag name patterns based on v3 BMRB format
            for tag, val in r.items():
                if "_Atom_chem_shift.Seq_ID" in tag and not seq_id:
                    try:
                        seq_id = int(val)
                    except:
                        pass
                elif "_Atom_chem_shift.Comp_ID" in tag:
                    comp_id = val.strip()
                elif "_Atom_chem_shift.Atom_ID" in tag:
                    atom_id = val.strip().upper()
                elif tag == "_Atom_chem_shift.Val":
                    try:
                        value = float(val)
                    except:
                        pass

            # If explicit Val not present, infer value as column after Atom_ID
            if value is None and atom_id:
                try:
                    # Recover ordered tags for this loop
                    order_s = r.get("__tags_order__", "")
                    if order_s:
                        order = order_s.split("\t")
                        # Find index of Atom_ID tag in order
                        atom_idx = -1
                        for j, t in enumerate(order):
                            if "_Atom_chem_shift.Atom_ID" in t:
                                atom_idx = j
                                break
                        if atom_idx != -1 and (atom_idx + 1) < len(order):
                            # The next column is the value by convention
                            next_tag = order[atom_idx + 1]
                            raw_val = r.get(next_tag)
                            if raw_val is not None and raw_val != "." and raw_val != "?":
                                value = float(raw_val)
                except Exception:
                    # best-effort; skip if cannot infer
                    pass
                        
            if seq_id is None or not comp_id or not atom_id or value is None:
                continue
                
            key = (seq_id, comp_id)
            bank = residue_to_atoms.setdefault(key, {})
            bank.setdefault(atom_id, []).append(value)
            
        except Exception as e:
            _vprint(f"[BMRB] Error parsing row (v3): {e}")
            continue

    # Build sequence from chemical shift data if not found in saveframe
    if not sequence:
        sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
        seq_letters: List[str] = []
        
        # minimal 3-letter to 1-letter mapping
        aa3_to1 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        
        for idx1, comp_id in sorted_keys:
            aa1 = aa3_to1.get(comp_id.upper(), "X")
            seq_letters.append(aa1)
        sequence = "".join(seq_letters)

    # Export raw chemical shift data to CSV
    _export_chem_shifts_to_csv(star_path, all_rows, saveframes)

    # Extract H, N, and CA shifts by sequence position
    H_by_pos: Dict[int, float] = {}
    N_by_pos: Dict[int, float] = {}
    CA_by_pos: Dict[int, float] = {}
    
    # Map residue numbers to sequence positions
    sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
    for seq_pos, (seq_id, comp_id) in enumerate(sorted_keys, 1):
        atoms = residue_to_atoms[(seq_id, comp_id)]
        
        # Aggregate H atoms (H, HN, H1)
        h_candidates: List[float] = []
        for name in ("H", "HN", "H1"):
            if name in atoms:
                h_candidates.extend(atoms[name])
        if h_candidates:
            H_by_pos[seq_pos] = statistics.median(h_candidates)
            
        # Aggregate N atoms
        n_candidates = atoms.get("N")
        if n_candidates:
            N_by_pos[seq_pos] = statistics.median(n_candidates)
        # Aggregate CA atoms
        ca_candidates = atoms.get("CA")
        if ca_candidates:
            CA_by_pos[seq_pos] = statistics.median(ca_candidates)

    _vprint(
        f"[BMRB] Sequence length {len(sequence)}; "
        f"H entries {len(H_by_pos)}, N entries {len(N_by_pos)}, CA entries {len(CA_by_pos)} (v3)"
    )
    return sequence, H_by_pos, N_by_pos, CA_by_pos


def parse_sequence_and_shifts_from_saveframes(star_path: str) -> List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]]:
    """Parse NMR-STAR file and return list of (sequence, H_shifts, N_shifts, saveframe_name) tuples.
    
    Returns all sequences from chemical shift saveframes that have both H and N data.
    """
    # Detect format and use appropriate parser
    format_version = _detect_bmrb_format(star_path)
    _vprint(f"[BMRB] Detected format version: {format_version}")
    
    if format_version == '3':
        sequence, H_shifts, N_shifts, CA_shifts = _parse_sequence_and_shifts_v3(star_path)
        if sequence and H_shifts and N_shifts:
            return [(sequence, H_shifts, N_shifts, CA_shifts, "assigned_chemical_shifts_1")]
        else:
            return []
    else:
        # Use original v21 parser
        with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
        lines = _tokenize_star_lines(text)
        
        # Parse all chemical shift saveframes
        saveframes = _parse_all_chem_shift_saveframes(lines)
        _vprint(f"[BMRB] Found {len(saveframes)} chemical shift saveframes in {os.path.basename(star_path)}")
        
        results: List[Tuple[str, Dict[int, float], Dict[int, float], Dict[int, float], str]] = []
        
        for saveframe_name, rows in saveframes:
            if not rows:
                continue
                
            _vprint(f"[BMRB] Processing {saveframe_name} with {len(rows)} rows")
        
            # Collect per residue atom shifts from this saveframe
            residue_to_atoms: Dict[Tuple[int, str], Dict[str, List[float]]] = {}
            for r in rows:
                try:
                    # Extract key fields - handle the actual BMRB tag names
                    seq_id = None
                    comp_id = ""
                    atom_id = ""
                    value = None
                    
                    # Try different tag name patterns based on actual BMRB format
                    for tag, val in r.items():
                        if "Residue_seq_code" in tag and not seq_id:
                            try:
                                seq_id = int(val)
                            except:
                                pass
                        elif "Residue_label" in tag:
                            comp_id = val.strip()
                        elif "Atom_name" in tag:
                            atom_id = val.strip().upper()
                        elif tag == "_Chem_shift_value":
                            try:
                                value = float(val)
                            except:
                                pass

                    # If explicit Chem_shift_value not present, infer value as column after Atom_name
                    if value is None and atom_id:
                        try:
                            # Recover ordered tags for this loop
                            order_s = r.get("__tags_order__", "")
                            if order_s:
                                order = order_s.split("\t")
                                # Find index of Atom_name tag in order
                                atom_idx = -1
                                for j, t in enumerate(order):
                                    if "Atom_name" in t:
                                        atom_idx = j
                                        break
                                if atom_idx != -1 and (atom_idx + 1) < len(order):
                                    # The next column is the value by convention
                                    next_tag = order[atom_idx + 1]
                                    raw_val = r.get(next_tag)
                                    if raw_val is not None and raw_val != "." and raw_val != "?":
                                        value = float(raw_val)
                        except Exception:
                            # best-effort; skip if cannot infer
                            pass
                                
                    if seq_id is None or not comp_id or not atom_id or value is None:
                        continue
                        
                    key = (seq_id, comp_id)
                    bank = residue_to_atoms.setdefault(key, {})
                    bank.setdefault(atom_id, []).append(value)
                    
                except Exception as e:
                    _vprint(f"[BMRB] Error parsing row: {e}")
                    continue

            # Build sequence from chemical shift data
            sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
            seq_letters: List[str] = []
            
            # minimal 3-letter to 1-letter mapping
            aa3_to1 = {
                "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            }
            
            for idx1, comp_id in sorted_keys:
                aa1 = aa3_to1.get(comp_id.upper(), "X")
                seq_letters.append(aa1)
            sequence = "".join(seq_letters)

            # Extract H, N, and CA shifts by sequence position
            H_by_pos: Dict[int, float] = {}
            N_by_pos: Dict[int, float] = {}
            CA_by_pos: Dict[int, float] = {}
            
            # Map residue numbers to sequence positions
            for seq_pos, (seq_id, comp_id) in enumerate(sorted_keys, 1):
                atoms = residue_to_atoms[(seq_id, comp_id)]
                
                # Aggregate H atoms (H, HN, H1)
                h_candidates: List[float] = []
                for name in ("H", "HN", "H1"):
                    if name in atoms:
                        h_candidates.extend(atoms[name])
                if h_candidates:
                    H_by_pos[seq_pos] = statistics.median(h_candidates)
                    
                # Aggregate N atoms
                n_candidates = atoms.get("N")
                if n_candidates:
                    N_by_pos[seq_pos] = statistics.median(n_candidates)
                # Aggregate CA atoms
                ca_candidates = atoms.get("CA")
                if ca_candidates:
                    CA_by_pos[seq_pos] = statistics.median(ca_candidates)

            # Only include sequences that have both H and N shifts
            if H_by_pos and N_by_pos:
                results.append((sequence, H_by_pos, N_by_pos, CA_by_pos, saveframe_name))
                _vprint(
                    f"[BMRB] Added sequence from {saveframe_name}: "
                    f"len={len(sequence)}, H={len(H_by_pos)}, N={len(N_by_pos)}, CA={len(CA_by_pos)}"
                )
            else:
                _vprint(
                    f"[BMRB] Skipping {saveframe_name}: missing H or N shifts "
                    f"(H={len(H_by_pos)}, N={len(N_by_pos)})"
                )
        
        return results


def parse_sequence_and_shifts(star_path: str) -> Tuple[str, Dict[int, float], Dict[int, float]]:
    """Extract sequence and per-residue H/N chemical shifts from NMR-STAR.

    Returns
    -------
    (sequence, H_by_pos, N_by_pos)
        sequence: str of one-letter AAs
        H_by_pos: map of 1-based sequence position → H shift (median if multiple)
        N_by_pos: map of 1-based sequence position → N shift (median if multiple)
    """
    # Detect format and use appropriate parser
    format_version = _detect_bmrb_format(star_path)
    _vprint(f"[BMRB] Detected format version: {format_version}")
    
    if format_version == '3':
        return _parse_sequence_and_shifts_v3(star_path)
    else:
        # Use original v21 parser
        with open(star_path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
        lines = _tokenize_star_lines(text)
        
        # Parse all chemical shift saveframes
        saveframes = _parse_all_chem_shift_saveframes(lines)
        _vprint(f"[BMRB] Found {len(saveframes)} chemical shift saveframes in {os.path.basename(star_path)}")
        
        # Combine data from all saveframes
        all_rows = []
        for saveframe_name, rows in saveframes:
            all_rows.extend(rows)
            _vprint(f"[BMRB] Parsed {len(rows)} rows from {saveframe_name}")
        
        _vprint(f"[BMRB] Total parsed {len(all_rows)} chem shift rows from {os.path.basename(star_path)}")

        # Extract sequence from monomeric polymer saveframe first
        sequence = _extract_sequence_from_saveframe(lines)
        if not sequence:
            _vprint(f"[BMRB] Warning: Could not extract sequence from saveframe, using CS data")
            sequence = ""

    # Collect per residue atom shifts from chemical shift data
    residue_to_atoms: Dict[Tuple[int, str], Dict[str, List[float]]] = {}
    for r in all_rows:
        try:
            # Extract key fields - handle the actual BMRB tag names
            seq_id = None
            comp_id = ""
            atom_id = ""
            value = None
            
            # Try different tag name patterns based on actual BMRB format
            for tag, val in r.items():
                if "Residue_seq_code" in tag and not seq_id:
                    try:
                        seq_id = int(val)
                    except:
                        pass
                elif "Residue_label" in tag:
                    comp_id = val.strip()
                elif "Atom_name" in tag:
                    atom_id = val.strip().upper()
                elif tag == "_Chem_shift_value":
                    try:
                        value = float(val)
                    except:
                        pass

            # If explicit Chem_shift_value not present, infer value as column after Atom_name
            if value is None and atom_id:
                try:
                    # Recover ordered tags for this loop
                    order_s = r.get("__tags_order__", "")
                    if order_s:
                        order = order_s.split("\t")
                        # Find index of Atom_name tag in order
                        atom_idx = -1
                        for j, t in enumerate(order):
                            if "Atom_name" in t:
                                atom_idx = j
                                break
                        if atom_idx != -1 and (atom_idx + 1) < len(order):
                            # The next column is the value by convention
                            next_tag = order[atom_idx + 1]
                            raw_val = r.get(next_tag)
                            if raw_val is not None and raw_val != "." and raw_val != "?":
                                value = float(raw_val)
                except Exception:
                    # best-effort; skip if cannot infer
                    pass
                        
            if seq_id is None or not comp_id or not atom_id or value is None:
                continue
                
            key = (seq_id, comp_id)
            bank = residue_to_atoms.setdefault(key, {})
            bank.setdefault(atom_id, []).append(value)
            
        except Exception as e:
            _vprint(f"[BMRB] Error parsing row: {e}")
            continue

    # Build sequence from chemical shift data if not found in saveframe
    if not sequence:
        sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
        seq_letters: List[str] = []
        
        # minimal 3-letter to 1-letter mapping
        aa3_to1 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        
        for idx1, comp_id in sorted_keys:
            aa1 = aa3_to1.get(comp_id.upper(), "X")
            seq_letters.append(aa1)
        sequence = "".join(seq_letters)

    # Export raw chemical shift data to CSV
    _export_chem_shifts_to_csv(star_path, all_rows, saveframes)

    # Extract H and N shifts by sequence position
    H_by_pos: Dict[int, float] = {}
    N_by_pos: Dict[int, float] = {}
    
    # Map residue numbers to sequence positions
    sorted_keys = sorted(residue_to_atoms.keys(), key=lambda k: k[0])
    for seq_pos, (seq_id, comp_id) in enumerate(sorted_keys, 1):
        atoms = residue_to_atoms[(seq_id, comp_id)]
        
        # Aggregate H atoms (H, HN, H1)
        h_candidates: List[float] = []
        for name in ("H", "HN", "H1"):
            if name in atoms:
                h_candidates.extend(atoms[name])
        if h_candidates:
            H_by_pos[seq_pos] = statistics.median(h_candidates)
            
        # Aggregate N atoms
        n_candidates = atoms.get("N")
        if n_candidates:
            N_by_pos[seq_pos] = statistics.median(n_candidates)

    _vprint(f"[BMRB] Sequence length {len(sequence)}; H entries {len(H_by_pos)}, N entries {len(N_by_pos)}")
    return sequence, H_by_pos, N_by_pos


def _extract_sequence_from_saveframe(lines: List[str]) -> str:
    """Extract sequence from monomeric polymer saveframe."""
    sequence = ""
    in_polymer_saveframe = False
    found_sequence_tag = False
    collecting_sequence = False
    seq_lines = []
    
    for i, line in enumerate(lines):
        line_stripped = line.strip()
        line_orig = line  # Keep original for multi-line sequence extraction
        
        # Look for monomeric polymer saveframe
        if line_stripped.startswith("save_") and "monomeric_polymer" in line_stripped:
            in_polymer_saveframe = True
            found_sequence_tag = False
            collecting_sequence = False
            seq_lines = []
            continue
            
        if in_polymer_saveframe and line_stripped == "save_":
            break
            
        if in_polymer_saveframe and line_stripped.startswith("_Mol_residue_sequence"):
            found_sequence_tag = True
            # Check if sequence is on the same line (after semicolon)
            if ";" in line_stripped:
                seq_part = line_stripped.split(";", 1)[1]
                if seq_part.strip():
                    seq_lines.append(seq_part)
                    collecting_sequence = True
            # Otherwise, sequence will be on following lines between semicolons
            continue
        
        # If we found the sequence tag, look for sequence lines
        if in_polymer_saveframe and found_sequence_tag:
            # Check if this line starts a semicolon block (opening semicolon)
            if line_stripped == ";" and not collecting_sequence:
                collecting_sequence = True
                continue
            # Check if this line ends the semicolon block (closing semicolon)
            elif line_stripped == ";" and collecting_sequence:
                # We've collected all sequence lines
                break
            # Collect sequence lines
            elif collecting_sequence:
                seq_lines.append(line_stripped)
    
    # Process collected sequence lines
    if seq_lines:
        # Join all sequence lines and remove whitespace
        seq_clean = "".join(seq_lines).replace(" ", "").replace("\n", "")
        
        # Convert 3-letter codes to 1-letter (basic mapping)
        aa3_to1 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        
        # Try to convert 3-letter codes, but also handle 1-letter codes
        result = []
        i = 0
        while i < len(seq_clean):
            if i + 3 <= len(seq_clean):
                three_letter = seq_clean[i:i+3].upper()
                if three_letter in aa3_to1:
                    result.append(aa3_to1[three_letter])
                    i += 3
                    continue
            # If not a 3-letter code, check if it's a valid 1-letter code
            one_letter = seq_clean[i].upper()
            if one_letter in "ACDEFGHIKLMNPQRSTVWY":
                result.append(one_letter)
            i += 1
        sequence = "".join(result)
                
    return sequence


