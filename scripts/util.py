# util.py
from Bio.Seq import Seq
from paths import *
import pickle
import urllib.request
import json
from tqdm import tqdm
import requests
from os.path import exists
from os import listdir
import os
import csv
import math
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser, PDBIO, Select
#import pymol
#from pymol import cmd, stored
from scipy.optimize import minimize
import re
from collections import defaultdict
import os
from os import listdir
from os.path import exists, isdir
import mdtraj as md
import pandas as pd
import ast
from scipy.stats import linregress
import statistics
from scipy import stats
import subprocess
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.Polypeptide import is_aa
import zipfile
from Bio.PDB import *
from typing import List

def convert_aa_name(residue_name: str) -> str:
    """
    Convert a 3-letter residue name to a one-letter amino-acid code.

    Handles common variants and non-standard residues; falls back to 'X' when unknown.
    """
    if residue_name is None:
        return 'X'
    name = str(residue_name).strip().upper()

    # Common non-standard/variant mappings
    variant_map = {
        'HIS': 'H', 'HSE': 'H', 'HSD': 'H', 'HSP': 'H',  # histidine protonation variants
        'MSE': 'M',  # selenomethionine -> methionine
        'SEC': 'U',  # selenocysteine
        'PYL': 'O',  # pyrrolysine
        'ASX': 'B',  # Asn/Asp ambiguous
        'GLX': 'Z',  # Gln/Glu ambiguous
        'UNK': 'X', 'XAA': 'X'
    }

    try:
        # Try Biopython canonical conversion first
        return Polypeptide.three_to_one(name)
    except Exception:
        pass

    # Try variant map
    if name in variant_map:
        return variant_map[name]

    # Strip common suffix/prefix decorations and retry
    # e.g., "ALA" with trailing digits or modified forms like "pSER" -> "SER"
    # Keep only last 3 uppercase letters if that looks like a standard code
    import re as _re
    m = _re.search(r'([A-Z]{3})$', name)
    if m:
        core = m.group(1)
        try:
            return Polypeptide.three_to_one(core)
        except Exception:
            if core in variant_map:
                return variant_map[core]

    return 'X'

def continue_prompt():
    while True:
        user_input = input("Do you want to continue? (y/n): ").lower()
        if user_input in ['n', 'no']:
            print("Terminating script.")
            return False
            exit()
        elif user_input in ['y', 'yes']:
            print("Continuing execution.")
            return True
            break
        else:
            print("Invalid input. Please enter 'y' for yes or 'n' for no.")

def get_pdb_file(pdb_id, output_dir = './bmrb_dat/'):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"{output_dir}{pdb_id}.pdb", 'w') as out_file:
                out_file.write(response.text)
        else:
            print(f"Couldn't download PDB file {pdb_id}")
            return None
        return f"{output_dir}{pdb_id}.pdb"
    except:
        return None


def execute_command_and_get_output(command):
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = result.stdout.decode().strip()  # decode stdout to string and remove trailing whitespace
    return output


def merge_chains(input_pdb, output_pdb):
    new_chain_id = "A"
    last_residue_id = 0
    current_chain_id = None

    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                residue_id = int(line[22:26].strip())

                if current_chain_id != chain_id:
                    residue_id_offset = last_residue_id
                    current_chain_id = chain_id

                new_residue_id = residue_id + residue_id_offset
                new_line = line[:21] + new_chain_id + f"{new_residue_id:>4}" + line[26:]
                f_out.write(new_line)
                last_residue_id = new_residue_id
            else:
                f_out.write(line)
    #print("wrote to " + output_pdb)

def get_pdb_sequence(pdb_path):

    sequence = []
    current_chain = ""
    prev_chain = ""
    prev_res_id = None

    with open(pdb_path, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                chain = line[21]
                res_name = line[17:20].strip()
                res_id = int(line[22:26])

                if prev_chain != "" and prev_chain != chain:
                    sequence.append(':')
                if prev_res_id != res_id:
                    sequence.append(convert_aa_name(res_name))

                prev_chain = chain
                prev_res_id = res_id

    return "".join(sequence)

def calculate_rmsd_matrix(ensemble):
    n_structures = len(ensemble)
    rmsd_matrix = np.zeros((n_structures, n_structures))

    for i in range(n_structures):
        for j in range(i+1, n_structures):
            rmsd = md.rmsd(ensemble[i], ensemble[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    return rmsd_matrix

def find_medoid_structure(pdb_files):

    if not pdb_files:
        print(pdb_files)
        raise ValueError("No PDB files found in the given directory")

    pdb_paths = [pdb_file for pdb_file in pdb_files if pdb_file.find('_A') == -1 and pdb_file.find('_B') == -1 and pdb_file.endswith('.pdb')]
    # print(pdb_paths)
    ensemble = []
    for pdb_path in pdb_paths:
        # print(pdb_path)
        ensemble.append(md.load(pdb_path))
    
    #ensemble = [md.load(pdb_path) for pdb_path in pdb_paths]

    rmsd_matrix = calculate_rmsd_matrix(ensemble)
    sum_rmsd = np.sum(rmsd_matrix, axis=1)
    medoid_index = np.argmin(sum_rmsd)

    return pdb_paths[medoid_index]

def get_bfactors(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('PDB', pdb_file)
    largest_chain = max(structure.get_chains(), key=lambda chain: len(list(chain.get_residues())))
    bfactors = []

    # Calculate average B-factor per residue
    for residue in largest_chain.get_residues():
        bfactor_sum = 0
        atom_count = 0
        for atom in residue.get_atoms():
            bfactor_sum += atom.get_bfactor()
            atom_count += 1
        average_bfactor = bfactor_sum / atom_count if atom_count > 0 else None
        bfactors.append(average_bfactor)

    return bfactors

def parse_list(value):
    try:
        return ast.literal_eval(value)
    except ValueError:
        return value
    except SyntaxError:
        return value

def parse_csv(file_name):
    data = []
    with open(file_name, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data.append({k: parse_list(v) for k, v in row.items()})
    return data

def write_to_csv(file_path, headers, data):
    with open(file_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(data)

def get_pdb_from_bmrb(bmrb_id):
    bmrb_url = "https://api.bmrb.io/v2/search/get_pdb_ids_from_bmrb_id/" + bmrb_id
    # Send a request to the BMRB API and retrieve the response
    with urllib.request.urlopen(bmrb_url) as url:
        responses = json.loads(url.read().decode())

    print(responses)
    try:
        k = responses[0]
    except:
        print(bmrb_id + " not in BMRB.")
        return
    # look for 'Author Provided' match_type
    for response in responses:
        if response["match_type"] == 'Author Provided':
            pdb_id = response["pdb_id"]
            print("PDB ID for BMRB ID", bmrb_id, "is", pdb_id)
            return pdb_id
    # if 'Author Provided' match_type not found, look for 'Exact' match_type
    for response in responses:
        if response["match_type"] == 'Exact':
            pdb_id = response["pdb_id"]
            print("PDB ID for BMRB ID", bmrb_id, "is", pdb_id)
            return pdb_id
        
    return None

def get_CSList(pdb_code):
    local_path = '/home/tiburon/Desktop/ROT4/complex3/CSLists/' + pdb_code + ".csv"
    if exists(local_path):
        print("Already have HSQC data for " + pdb_code + ", continuing...")

    # Define the URL to access the BMRB API for the given PDB code
    # example url = https://api.bmrb.io/v2/search/get_bmrb_ids_from_pdb_id/2JNS
    bmrb_url = "https://api.bmrb.io/v2/search/get_bmrb_ids_from_pdb_id/" + pdb_code

    # Send a request to the BMRB API and retrieve the response
    with urllib.request.urlopen(bmrb_url) as url:
        responses = json.loads(url.read().decode())

    try:
        k = responses[0]
    except:
        print(pdb_code + " not in BMRB.")
        return

    # Check if the response contains a BMRB ID for the given PDB code
    if "bmrb_id" in responses[0]:
        bmrb_id = None
        for response in responses:
            if response["match_types"][0] == 'Exact':
                bmrb_id = response["bmrb_id"]

        if bmrb_id is None:
            print("Couldn't find exact match.")
            return

        print("BMRB ID for PDB code", pdb_code, "is", bmrb_id)

        # get Chemical shift list from bmrb
        bmrb_csl_url = "https://api.bmrb.io/v2/entry/"+str(bmrb_id)

        # Send a request to the BMRB API and retrieve the response
        with urllib.request.urlopen(bmrb_csl_url) as url:
            responses = json.loads(url.read().decode())
            saveframe_data = responses[bmrb_id]['saveframes']
            cs_frame_index = -1
            for i, data in enumerate(saveframe_data):
                #if data['name'] in ["assigned_chemical_shifts_1", "assigned_chem_shift_list_1", "chemical_shift_1", "ChemShifts"]:
                if data['category'] == "assigned_chemical_shifts":
                    cs_frame_index = i
                    break
            if cs_frame_index == -1:
                print("Couldn't find chemical shift lists... continuing...")
                return
            CSL_data = saveframe_data[cs_frame_index]['loops']
            loop_i = -1
            for i in range(0, len(CSL_data)):
                if CSL_data[i]['category'] == '_Atom_chem_shift':
                    loop_i = i
            if loop_i == -1:
                return

            CSL_data = CSL_data[loop_i]
            tags = CSL_data['tags']
            data = CSL_data['data']
            dat = [ d for d in data]

            write_to_csv(local_path, tags, dat)
            print("wrote to file " + local_path)

    else:
        print("No BMRB ID found for PDB code", pdb_code)
        return

def update_row(csv_filename, apo, bound, new_values, new_columns):
    try:
        # Load the DataFrame if the CSV file exists
        df = pd.read_csv(csv_filename, low_memory=False)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        # Create an empty DataFrame if the CSV file is empty or doesn't exist
        df = pd.DataFrame()

    data_dict = {col: val for col, val in zip(new_columns, new_values)}
    #print(data_dict)
    #print(data_dict.items())

    # Check if 'apo' and 'bound' columns exist
    if 'apo_bmrb' not in df.columns or 'holo_pdb' not in df.columns:
        df = df._append(data_dict, ignore_index=True)
    else:
        # Update or create the row
        row_index = df[(df['holo_pdb'] == bound)].index
        print("UPDATING ROW INDEX : " + str(row_index))
        if not row_index.empty:
            for col, val in data_dict.items():
                try:
                    df.loc[row_index[0], col] = val
                except:
                    data_dict['apo_bmrb'] = apo
                    data_dict['holo_pdb'] = bound
                    df = df._append(data_dict, ignore_index=True)
        else:
            data_dict['apo_bmrb'] = apo
            data_dict['holo_pdb'] = bound
            df = df._append(data_dict, ignore_index=True)

    # Save the DataFrame back to the CSV file
    # Format float columns to 4 decimal places
    float_columns = df.select_dtypes(include=['float64', 'float32']).columns
    for col in float_columns:
        df[col] = df[col].round(4)
    df.to_csv(csv_filename, index=False)

def rmsd(offset, spectrum1, spectrum2):
    aligned_spectrum1 = spectrum1 + offset
    return np.sqrt(np.mean((aligned_spectrum1 - spectrum2) ** 2))

def find_optimal_offset(spectrum1, spectrum2):
    spectrum1 = np.array(spectrum1)
    spectrum2 = np.array(spectrum2)
    initial_offset = 0
    result = minimize(rmsd, initial_offset, args=(spectrum1, spectrum2), method='L-BFGS-B')
    optimal_offset = result.x[0]
    return optimal_offset

def compute_DockQ_score(pdb_file1, pdb_file2, well_defined_res=None):
    model_pdb = str(pdb_file1)
    native_pdb = str(pdb_file2)
    DockQ = None
    s = 'DockQ \'' + model_pdb + '\' \'' + native_pdb + '\''
    print(s)
    output = execute_command_and_get_output(s).split('\n')

    try:
        # Parse using the actual output format
        DockQ = float([f.split(': ')[1] for f in output if f.strip().startswith('DockQ:')][0])
        iRMS = float([f.split(': ')[1] for f in output if f.strip().startswith('iRMSD:')][0])
        LRMS = float([f.split(': ')[1] for f in output if f.strip().startswith('LRMSD:')][0])
        Fnat = float([f.split(': ')[1] for f in output if f.strip().startswith('fnat:')][0])
        Fnonnat = float([f.split(': ')[1] for f in output if f.strip().startswith('fnonnat:')][0])
        F1 = float([f.split(': ')[1] for f in output if f.strip().startswith('F1:')][0])
        clashes = int([f.split(': ')[1] for f in output if f.strip().startswith('clashes:')][0])
        
        if True:
            print("Fnat = " + str(Fnat))
            print("Fnonnat = " + str(Fnonnat))
            print("iRMS = " + str(iRMS))
            print("LRMS = " + str(LRMS))
            print("DockQ = " + str(DockQ))
            print("F1 = " + str(F1))
            print("Clashes = " + str(clashes))
    except Exception as e:
        print(e)
        raise

    return iRMS, LRMS, DockQ, Fnat, Fnonnat, F1, clashes

def get_es_files(holo: str) -> List[str]:
    """Get ensemble selection shift files for a holo structure"""
    es_dir = PDB_FILES + holo.lower() + '_max_RPF_NLDR_CSPRank_files'
    if not exists(es_dir):
        return []
        
    basenames = [f[:f.rfind('.')] for f in listdir(es_dir)]
    shift_files = []
    
    for basename in basenames:
        # Try AFS directory first
        shift_file = f"{CS_Predictions}{holo}_AFS_shift_predictions/{basename}.csv"
        if not exists(shift_file):
            # Try AFS2 directory if not found
            shift_file = f"{CS_Predictions}{holo}_AFS2_shift_predictions/{basename}.csv"
            if not exists(shift_file):
                continue
        shift_files.append(shift_file)
        
    return shift_files

def get_apo_medoid_files(apo: str, holo: str) -> List[str]:
    """
    Get medoid chemical shift files for apo structure
    
    Args:
        apo: Apo structure identifier
        holo: Holo structure identifier
        
    Returns:
        List of paths to medoid shift files
    """
    # Try to find medoid file specific to this holo structure
    medoid_files = [f for f in listdir(apo_NMR_shift_dir) 
                    if f.find(holo) != -1 and f.find('medoid') != -1 and f.endswith('.csv')]
    if medoid_files:
        return [apo_NMR_shift_dir + f for f in medoid_files]
    
    # Try general apo medoid files
    medoid_files = [f for f in listdir(apo_NMR_shift_dir) 
                    if f.find(apo) != -1 and f.find('medoid') != -1 and f.endswith('.csv')]
    if medoid_files:
        return [apo_NMR_shift_dir + f for f in medoid_files]
        
    return []

def get_holo_medoid_files(holo: str) -> List[str]:
    """
    Get medoid chemical shift files for holo structure
    
    Args:
        holo: Holo structure identifier
        
    Returns:
        List of paths to medoid shift files
    """
    medoid_files = [f for f in listdir(holo_NMR_shift_dir) 
                    if f.find(holo) != -1 and f.find('medoid') != -1 and f.endswith('.csv')]
    return [holo_NMR_shift_dir + f for f in medoid_files]

def get_apo_files(apo: str, holo: str, source: str = 'NMR') -> List[str]:
    """
    Get chemical shift files for apo structure
    
    Args:
        apo: Apo structure identifier
        holo: Holo structure identifier
        source: Source of structure ('NMR' or 'AF2')
        
    Returns:
        List of paths to shift files
    """
    if source == 'AF2':
        # Get AF2 apo files
        return [apo_AF2_shift_dir + f for f in listdir(apo_AF2_shift_dir) 
                if f.find(apo) != -1 and f.endswith('.csv')]
    else:  # NMR files
        # First try files specific to this holo structure
        apo_files = [f for f in listdir(apo_NMR_shift_dir) 
                     if f.find(holo) != -1 and f.endswith('.csv')]
        if apo_files:
            return [apo_NMR_shift_dir + f for f in apo_files]
        
        # Then try general apo files
        apo_files = [f for f in listdir(apo_NMR_shift_dir) 
                     if f.find(apo) != -1 and len(f.split('_')) == 2 and f.endswith('.csv')]
        if apo_files:
            return [apo_NMR_shift_dir + f for f in apo_files]
            
        return []

def get_holo_files(holo: str, structure_source: str) -> List[str]:
    """Get chemical shift files for holo structure based on structure source"""
    if structure_source == 'NMR':
        # Get all NMR shift files for this structure
        return [holo_NMR_shift_dir + f for f in listdir(holo_NMR_shift_dir) 
                if f.find(holo) != -1 and f.endswith('.csv')]
    elif structure_source == 'AF2':
        # Get all AF2 shift files
        return [holo_AF2_shift_dir + f for f in listdir(holo_AF2_shift_dir) 
                if f.find(holo) != -1 and f.endswith('.csv')]
    else:
        return []

def get_shift_files(apo, holo, structure_source, CSmethod, basename = None):

    apo_shift_file = ""
    holo_shift_file = ""
    def get_apo_pred_files(apo, holo, medoid = False):
        
        if medoid:
            apo_pdb_files = [ NMR_apo_structure_dir + f for f in listdir(NMR_apo_structure_dir) if f.find(apo) != -1]

            new_pdb_files = []
            for pdb_file in apo_pdb_files:
                basename = pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('.')]
                file = [ apo_NMR_shift_dir + f for f in listdir(apo_NMR_shift_dir) if f.find(basename) != -1 and f.endswith('.csv')]
                if len(file) > 0:
                    new_pdb_files.append(pdb_file)
            apo_pdb_files = new_pdb_files
            apo_pdb_files_from_holo = [ apo for apo in apo_pdb_files if apo.find(holo) != -1]
            if len(apo_pdb_files_from_holo) > 0:
                apo_pdb_files = apo_pdb_files_from_holo
            print("APO PDB FILES")
            print(apo_pdb_files)
            medoid_file = find_medoid_structure(apo_pdb_files)
            medoid_file = medoid_file[medoid_file.rfind('/')+1:medoid_file.rfind('.')]
            apo_shift_files = [ apo_NMR_shift_dir + f for f in listdir(apo_NMR_shift_dir) if f.find(medoid_file) != -1 and f.endswith('.csv')]
            print("APO_SHIFT FILES")
            print(apo_shift_files)
            return apo_shift_files
        else:
            apo_shift_files = [apo_NMR_shift_dir + f for f in listdir(apo_NMR_shift_dir) if f.find(holo) != -1 and f.endswith('.csv')]
            if len(apo_shift_files) > 0:
                return apo_shift_files
            
            apo_shift_files = [apo_NMR_shift_dir + f for f in listdir(apo_NMR_shift_dir) if f.find(apo) != -1 and len(f.split('_')) == 2 and f.endswith('.csv')]
            if len(apo_shift_files) > 0:
                return apo_shift_files

            print("trying to find apo files for pdb id " + holo)        
            apo_shift_pdb = [f for f in listdir(apo_NMR_shift_dir) if f.find(apo) != -1 and f.endswith('.csv')][0].split('_')[1].strip()
            print("found alternative holo to look for apo files " + apo_shift_pdb)
            return get_apo_pred_files(apo, apo_shift_pdb)
        
    if structure_source.find("NMR") != -1 and structure_source.find('real') != -1:
        if CSmethod != "REAL":
            print("malformed structure_source/CSmethod combination")
            raise
        apo_shift_file = real_CSList_dir + apo.upper() +'.csv'
        holo_shift_file = real_CSList_dir + holo.upper() +'.csv'
    elif structure_source.find("NMR") != -1 and structure_source.find('pred') != -1:
        if structure_source.find('medoid') != -1:
            if CSmethod == "UCBShift":
                apo_shift_file = get_apo_pred_files(apo, holo, medoid = True)
                try_apo_files = [ apo_file for apo_file in apo_shift_file if apo_file.find(holo) != -1]
                if len(try_apo_files) > 0:
                    apo_shift_file = try_apo_files
                holo_pdb_files = [ NMR_holo_structure_dir + f for f in listdir(NMR_holo_structure_dir) if f.find(holo) != -1]
                new_holo_pdb_files = []
                for pdb_file in holo_pdb_files:
                    basename = pdb_file[pdb_file.rfind('/')+1:pdb_file.rfind('.')]
                    file = [ holo_NMR_shift_dir + f for f in listdir(holo_NMR_shift_dir) if f.find(basename) != -1 and f.endswith('.csv')]
                    if len(file) > 0:
                        new_holo_pdb_files.append(pdb_file)
                holo_pdb_files = new_holo_pdb_files
                if len(holo_pdb_files) == 0:
                    raise
                medoid_file = find_medoid_structure(holo_pdb_files)
                medoid_basename = medoid_file[medoid_file.rfind('/')+1:medoid_file.rfind('.')]
                holo_shift_file = [ holo_NMR_shift_dir + f for f in listdir(holo_NMR_shift_dir) if f.find(medoid_basename) != -1 and f.endswith('.csv')]
                print("NMR SHIFT FILES")
                print("MEDOID FILE = " + str(medoid_file))
                print("MEDOID APO SHIFT FILE = " + str(apo_shift_file))
                print("MEDOID HOLO SHIFT FILE = " + str(holo_shift_file))
        else:
            if CSmethod == "UCBShift":
                apo_shift_file = get_apo_pred_files(apo, holo)
                try_apo_files = [ apo_file for apo_file in apo_shift_file if apo_file.find(holo) != -1]
                if len(try_apo_files) > 0:
                    apo_shift_file = try_apo_files
                holo_shift_file =[ holo_NMR_shift_dir + f for f in listdir(holo_NMR_shift_dir) if f.find(holo) != -1 and f.endswith('.csv')]
            elif CSmethod == "SPARTA":
                apo_shift_file = apo_AF2_SPARTA_dir + apo +'.tab'
                holo_shift_file = holo_AF2_SPARTA_dir + holo +'.tab'
            elif CSmethod == "ShiftX":
                apo_shift_file = apo_AF2_ShiftX_dir + apo +'.pdb.cs'
                holo_shift_file = holo_AF2_ShiftX_dir + holo +'.pdb.cs'
            else:
                print("malformed structure_source/CSmethod combination")
                raise
    elif structure_source.find("AF2") != -1:
        if structure_source.find('top_rank') != -1:
            apo_shift_file = get_apo_pred_files(apo, holo, medoid = True)
            holo_shift_file = [ holo_AF2_shift_dir + f for f in listdir(holo_AF2_shift_dir) if f.find(holo + "_1") != -1 and f.endswith('.csv')]
            print("AF2 SHIFT FILES")
            print("MEDOID APO SHIFT FILE = " + str(apo_shift_file))
            print("MEDOID HOLO SHIFT FILE = " + str(holo_shift_file))
        else:
            if CSmethod == "UCBShift":
                apo_shift_file = get_apo_pred_files(apo, holo)
                try_apo_files = [ apo_file for apo_file in apo_shift_file if apo_file.find(holo) != -1]
                if len(try_apo_files) > 0:
                    apo_shift_file = try_apo_files
                holo_shift_file = [ holo_AF2_shift_dir + f for f in listdir(holo_AF2_shift_dir) if f.find(holo) != -1 and f.endswith('.csv')]
                #holo_shift_file = holo_AF2_shift_dir + holo +'.csv'
            elif CSmethod == "SPARTA":
                apo_shift_file = apo_AF2_SPARTA_dir + apo +'.tab'
                holo_shift_file = holo_AF2_SPARTA_dir + holo +'.tab'
            elif CSmethod == "ShiftX":
                apo_shift_file = apo_AF2_ShiftX_dir + apo +'.pdb.cs'
                holo_shift_file = holo_AF2_ShiftX_dir + holo +'.pdb.cs'
            else:
                print("malformed structure_source/CSmethod combination")
                raise
    elif structure_source.find("AF3") != -1:
        if CSmethod == "UCBShift":
            apo_shift_file = apo_AF3_shift_dir + apo +'.csv'
            holo_shift_file = holo_AF3_shift_dir + holo +'.csv'
        elif CSmethod == "SPARTA":
            apo_shift_file = apo_AF3_SPARTA_dir + apo +'.tab'
            holo_shift_file = holo_AF3_SPARTA_dir + holo +'.tab'
        elif CSmethod == "ShiftX":
            apo_shift_file = apo_AF3_ShiftX_dir + apo +'.pdb.cs'
            holo_shift_file = holo_AF3_ShiftX_dir + holo +'.pdb.cs'
        else:
            print("malformed structure_source/CSmethod combination")
            raise
    elif structure_source.find('ES') != -1:
        if CSmethod == 'UCBShift':
            apo_shift_file = get_apo_pred_files(apo, holo)
            if basename == "" or basename is None:

                ES_dir = PDB_FILES + holo.lower() + '_max_RPF_NLDR_CSPRank_files'
                basenames = [f[:f.rfind('.')] for f in listdir(ES_dir) if f.find(apo) == -1 ]
                holo_shift_file = []
                for basename in basenames:
                    t_holo_shift_file = CS_Predictions + holo + '_AFS_shift_predictions/'+basename+'.csv'
                    if not(exists(t_holo_shift_file)):
                        t_holo_shift_file = CS_Predictions + holo + '_AFS2_shift_predictions/'+basename+'.csv'
                        if not(exists(t_holo_shift_file)):
                            print("Could not locate " + basename + " ES shift file in AFS_ or AFS2_shift_predictions directories.")
                            raise
                        else:
                            holo_shift_file.append(t_holo_shift_file)
                    else:
                        holo_shift_file.append(t_holo_shift_file)
            else:             
                if type(basename) == str:   
                    holo_shift_file = CS_Predictions + holo + '_AFS_shift_predictions/'+basename+'.csv'
            if not(exists(holo_shift_file)):
                holo_shift_file = CS_Predictions + holo + '_AFS2_shift_predictions/'+basename+'.csv'
                if not(exists(holo_shift_file)):
                    print("Could not locate " + basename + " ES shift file in AFS_ or AFS2_shift_predictions directories.")
                    raise
                elif type(basename) == list:
                    holo_shift_file = []
                    for base in basename:
                        t_holo_shift_file = CS_Predictions + holo + '_AFS_shift_predictions/'+base+'.csv'
                        if not(exists(t_holo_shift_file)):
                            t_holo_shift_file = CS_Predictions + holo + '_AFS2_shift_predictions/'+base+'.csv'
                            if not(exists(t_holo_shift_file)):
                                print("Could not locate " + base + " ES shift file in AFS_ or AFS2_shift_predictions directories.")
                                raise
                            else:
                                holo_shift_file.append(t_holo_shift_file)
                        else:
                            holo_shift_file.append(t_holo_shift_file)
                else:
                    print("malformed basename")
                    raise
        else:
            print("Unimplemented CS Prediction method and Structure Source.")
            raise
    else:
        print("malformed structure_source.")
        raise

    return apo_shift_file, holo_shift_file

def calc_CSP_wrapper(apo, holo, well_defined_res, method = "MONTE", CSmethod = "REAL", structure_source = "NMR_real", match_seq = "", basename=""):

    apo_shift_files = ""
    holo_shift_files = ""
    apo_shift_files, holo_shift_files = get_shift_files(apo, holo, structure_source, CSmethod, basename=basename)

    if type(apo_shift_files) == str:
        apo_shift_files = [apo_shift_files]
    if type(holo_shift_files) == str:
        holo_shift_files = [holo_shift_files]

    print("CALCULATING CSPs FOR " + apo + " AND " + holo)
    print("METHOD = " + method)
    print("CS METHOD = " + CSmethod)
    print("STRUCTURE SOURCE = " + structure_source)

    # print("APO SHIFT FILES")
    # print(apo_shift_files)
    # print("HOLO SHIFT FILES")
    # print(holo_shift_files)

    if len(apo_shift_files) == 0 or len(holo_shift_files) == 0:
        print("No shift files found for " + apo + " or " + holo)
        return None
    
    return calc_CSP_from_ensemble(apo_shift_files, holo_shift_files, method=method, UCBShift_pred_holo=(CSmethod == "UCBShift"), UCBShift_pred_apo=(CSmethod == "UCBShift"))

def get_value_from_csv(csv_file, bound, apo, column):
    try:
        # Load the DataFrame if the CSV file exists
        df = pd.read_csv(csv_file)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        # Create an empty DataFrame if the CSV file is empty or doesn't exist
        df = pd.DataFrame()

    # Filter rows where 'bound' and 'apo' match the input
    matching_rows = df[(df['bound'] == bound) & (df['apo'] == apo)]

    # If there's no match, return None
    if matching_rows.empty:
        return None

    # If there's more than one match, this will return the value from the first match
    return matching_rows.iloc[0][column]


def get_longest_chain(pdb_filepath):
    with open(pdb_filepath, 'r') as pdb_file:
        lines = pdb_file.readlines()

    chain_residue_count = defaultdict(set)
    atom_lines = []
    seqs = {}
    last_res_ind = -1
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_lines.append(line)
            chain_id = line[21]
            if chain_id not in list(seqs):
                seqs[chain_id] = ""
            residue_index = int(line[22:26].strip())
            if residue_index != last_res_ind:
                last_res_ind = residue_index
                #print(line[17:20].strip())
                seqs[chain_id] += convert_aa_name(line[17:20].strip())
            chain_residue_count[chain_id].add(residue_index)

    # Identify the longest chain
    longest_chain = max(chain_residue_count, key=lambda k: len(chain_residue_count[k]))
    return longest_chain

def get_shortest_chain(pdb_filepath):
    with open(pdb_filepath, 'r') as pdb_file:
        lines = pdb_file.readlines()

    chain_residue_count = defaultdict(set)
    atom_lines = []
    seqs = {}
    last_res_ind = -1
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_lines.append(line)
            chain_id = line[21]
            if chain_id not in list(seqs):
                seqs[chain_id] = ""
            residue_index = int(line[22:26].strip())
            if residue_index != last_res_ind:
                last_res_ind = residue_index
                #print(line[17:20].strip())
                seqs[chain_id] += convert_aa_name(line[17:20].strip())
            chain_residue_count[chain_id].add(residue_index)

    # Identify the longest chain
    longest_chain = min(chain_residue_count, key=lambda k: len(chain_residue_count[k]))
    return longest_chain


def update_b_factors_longest_chain(pdb_filepath, b_factors, bound_seq, new_pdb_filepath):
    with open(pdb_filepath, 'r') as pdb_file:
        lines = pdb_file.readlines()

    chain_residue_count = defaultdict(set)
    atom_lines = []
    seqs = {}
    last_res_ind = -1
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_lines.append(line)
            chain_id = line[21]
            if chain_id not in list(seqs):
                seqs[chain_id] = ""
            residue_index = int(line[22:26].strip())
            if residue_index != last_res_ind:
                last_res_ind = residue_index
                #print(line[17:20].strip())
                seqs[chain_id] += convert_aa_name(line[17:20].strip())
            chain_residue_count[chain_id].add(residue_index)

    # Identify the longest chain
    longest_chain = max(chain_residue_count, key=lambda k: len(chain_residue_count[k]))
    new_bfactors = []

    #print(seqs[longest_chain])
    if len(chain_residue_count[longest_chain]) != len(b_factors):
        # align
        try:
            bound_aligned1, bound_aligned2 = align(bound_seq, seqs[longest_chain])
            print(bound_aligned1)
            print(bound_aligned2)
            print("HERE")

        except:
            return
        ind = 0
        for i,c in enumerate(bound_aligned1):
            
            if bound_aligned2[i] in ['_', '-']:
                ind += 1
                continue
            if bound_aligned1[i] in ['_', '-']:
                ind += 1
                new_bfactors.append(-1)
                continue
            #if ind >= len(bound_seq) or ind >= len(b_factors) or c != bound_seq[ind]:
            #    new_bfactors.append(0)
            #    x = 0
            #else:
            new_bfactors.append(b_factors[ind])
            ind += 1
        if len(new_bfactors) != len(chain_residue_count[longest_chain]):
            return
    else:
        new_bfactors = [l for l in b_factors]

    b_factor_dict = {index: b_factor for index, b_factor in zip(sorted(chain_residue_count[longest_chain]), new_bfactors)}
    #print(b_factor_dict)

    updated_lines = []
    for line in lines:
        if (line.startswith('ATOM') or line.startswith('HETATM')) and line[21] == longest_chain:
            residue_index = int(line[22:26].strip())
            new_b_factor = b_factor_dict[residue_index]
            updated_line = line[:60] + f'{new_b_factor:6.2f}' + line[66:]
            #print(new_b_factor)
            #print(updated_line)
            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)

    #print(updated_lines)
    print(new_pdb_filepath)
    with open(new_pdb_filepath, 'w') as pdb_file:
        pdb_file.writelines(updated_lines)

def get_interface_patch(contacts):
    return set(i for contact in contacts for i in contact)

def calculate_ics_ips(pdb_file1, pdb_file2, threshold = 5.0):

    contacts1 = get_interface_contacts(pdb_file1, threshold)
    contacts2 = get_interface_contacts(pdb_file2, threshold)

    # Calculate ICS
    precision = len(contacts1.intersection(contacts2)) / len(contacts1)
    recall = len(contacts2.intersection(contacts1)) / len(contacts2)
    ics = 0
    try:
        ics = 2 * (precision * recall) / (precision + recall)
    except:
        ics = 0

    # Calculate IPS
    patch1 = get_interface_patch(contacts1)
    patch2 = get_interface_patch(contacts2)
    intersection_size = len(patch1.intersection(patch2))
    union_size = len(patch1.union(patch2))
    ips = 0
    try:
        ips = intersection_size / union_size
    except:
        ips = 0

    return ics, ips
def get_atom_distance(atom1, atom2):
    vec = atom1-atom2
    distance = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
    return distance

def get_closest_distance(residue1, residue2):
    min_distance = 1000
    for atom1 in residue1:
        for atom2 in residue2:
            distance = get_atom_distance(atom1.coord, atom2.coord)
#            print(distance)
            if distance < min_distance:
                min_distance = distance
    return min_distance

def get_interface_contacts(pdb_file, threshold):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    residues = [residue for model in structure for chain in model for residue in chain]
    num_residues = len(residues)

    pairs = []
    for i in range(num_residues):
        residue1 = residues[i]
        chain1 = residue1.parent.id
        for j in range(i+1, num_residues):
            residue2 = residues[j]
            chain2 = residue2.parent.id

            if chain1 != chain2:
                distance = get_closest_distance(residue1, residue2)
                if distance <= threshold:
                    pairs.append((i, j))

    ret = set(pairs)
    return ret
def compute_lddt_score(pdb_file1, pdb_file2):
    # need to merge chains because our version of lddt only accepts single chain files
    merged_file1 = pdb_file1[:len(pdb_file1)-4]+'_merged.pdb'
    merged_file2 = pdb_file2[:len(pdb_file2)-4]+'_merged.pdb'
    merge_chains(pdb_file1, merged_file1)
    merge_chains(pdb_file2, merged_file2)

    clistring = './lddt ' + merged_file1 + ' ' + merged_file2
    print(clistring)

    result = subprocess.run(['./lddt', merged_file1, merged_file2],
                            capture_output=True, text=True)


    # Parse the global LDDT score from the output
    global_lddt_score = None
    for line in result.stdout.splitlines():
        if "Global LDDT score:" in line:
            global_lddt_score = float(line.split()[3])
            break

    if global_lddt_score is None:
        raise ValueError("Failed to parse global LDDT score from lddt output")

    print("lddt score = " + str(global_lddt_score))
    os.system('rm ' + merged_file1)
    os.system('rm ' + merged_file2)

    return global_lddt_score


def compute_TM_score(pdb_file1, pdb_file2):
    # need to merge chains because our version of lddt only accepts single chain files
    merged_file1 = pdb_file1[:len(pdb_file1)-4]+'_merged.pdb'
    merged_file2 = pdb_file2[:len(pdb_file2)-4]+'_merged.pdb'
    merge_chains(pdb_file1, merged_file1)
    merge_chains(pdb_file2, merged_file2)

    clistring = './TMalign ' + merged_file1 + ' ' + merged_file2
    print(clistring)

    # Call TMalign with the model and reference PDB files
    result = subprocess.run(['./TMalign', merged_file1, merged_file2],
                            capture_output=True, text=True)

    # Parse the TM-score from the output
    tm_score = None
    for line in result.stdout.splitlines():
        if "TM-score=" in line and "LN=" in line:
            tm_score = float(line.split()[1])
            break

    if tm_score is None:
        raise ValueError("Failed to parse TM-score from TMalign output")

    print("TM score = " + str(tm_score))

    os.system('rm ' + merged_file1)
    os.system('rm ' + merged_file2)

    return tm_score


def compute_structure_similarity(pdb_file1, pdb_file2, multimer=True):
    """
    Computes structural similarity between two PDB files using USalign.
    Optimized for protein-peptide complexes with multimeric support.
    
    Args:
        pdb_file1 (str): Path to first PDB file
        pdb_file2 (str): Path to second PDB file  
        multimer (bool): Whether to treat structures as multimeric complexes
        
    Returns:
        float: Similarity score (TM-score)
    """
    try:
        # Build USalign command with appropriate options
        # -mol prot: only align proteins
        # -mm 1: alignment of multi-chain oligomeric structures
        # -ter 1: align all chains of the first model (for asymmetric units)
        # -het 1: align both ATOM and HETATM residues (important for peptides)
        base_cmd = 'USalign'
        options = ['-mol prot']
        
        if multimer:
            options.extend(['-mm 1', '-ter 1'])
        options.append('-het 1')
        
        clistring = f'{base_cmd} {pdb_file1} {pdb_file2} {" ".join(options)}'
        print(clistring)
        result = subprocess.run(clistring, shell=True, capture_output=True, text=True)
        
        # Parse TM-score from output
        tm_score = None
        for line in result.stdout.splitlines():
            if "TM-score=" in line:
                # USalign outputs multiple TM-scores, we want the first one
                # which is normalized by the length of the first structure
                tm_score = float(line.split()[1])
                break
                
        if tm_score is None:
            raise ValueError("Failed to parse TM-score from USalign output")
            
        return tm_score
        
    except Exception as e:
        print(f"Error computing structure similarity: {e}")
        print("USalign output:")
        print(result.stdout)
        print("USalign error:")
        print(result.stderr)
        raise