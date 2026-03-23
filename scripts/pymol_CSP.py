# open_comp_exp_csp.py
from Bio.PDB import PDBParser, PDBIO, Select
import pymol
from pymol import cmd, stored
from os.path import exists
import os
from os import listdir
from os.path import isdir
import csv
from tqdm import tqdm
import os
from os import listdir
from os.path import exists, isdir
import mdtraj as md
import numpy as np
from util import *
import glob

def split_pdb_chains(pdb_filepath):
    # Parse the input PDB file
    parser = PDBParser()
    structure = parser.get_structure('input_structure', pdb_filepath)
    new_files = []

    # Iterate through chains and write each chain to a separate PDB file
    io = PDBIO()
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            output_filename = f'{os.path.splitext(pdb_filepath)[0]}_{chain_id}.pdb'
            new_files.append(output_filename)
            io.set_structure(chain)
            io.save(output_filename, select=Select())
    return new_files

def open_structures(bound_file, apo, label, CSPs, bound_seq):
    print(CSPs)
    bound = bound_file[bound_file.rfind('/')+1:bound_file.rfind('/')+5].lower()
    new_exp_file = bound_file[:bound_file.rfind('.')]+'_'+apo+'.pdb'
    print(new_exp_file)
    update_b_factors_longest_chain(bound_file, CSPs, bound_seq, new_exp_file)
    pdb = bound
    bound_path = new_exp_file
    bound = pdb
    pymol.finish_launching()
    new_files = split_pdb_chains(bound_path)
    file_sizes = [os.path.getsize(file_path) for file_path in new_files]

    # load in the larger of the two files:
    protein_file = new_files[file_sizes.index(max(file_sizes))]
    peptide_file = new_files[file_sizes.index(min(file_sizes))]

    #pymol.finish_launching()
    cmd.load(bound_path, label)

    if file_sizes.index(max(file_sizes)) == 0:
        pymol.cmd.alter(label + " and chain A", "chain='prot"+label+"'")
        pymol.cmd.alter(label + " and chain B", "chain='peptide"+label+"'")
    else:
        pymol.cmd.alter(label + " and chain B", "chain='prot"+label+"'")
        pymol.cmd.alter(label + " and chain A", "chain='peptide"+label+"'")

    # Zoom to chain A
    cmd.orient()

    cmd.viewport(800, 800)
    return new_files

def hide_structures(label, chain):
    selection_string = f"{label} and chain {chain}"
    pymol.cmd.hide('everything', selection_string)

def show_structures(label, chain):
    selection_string = f"{label} and chain {chain}"
    pymol.cmd.show('cartoon', selection_string)

def hide_CSP_surface(label, chain):
    selection_string = f"{label} and chain {chain}"
    pymol.cmd.hide('surface', selection_string)
    cmd.set('transparency', 0, selection_string)

def show_CSP_surface(label, chain, cutoffs):
    selection_string = f"{label} and chain {chain}"
    pymol.cmd.show('surface', selection_string)
    cmd.set('transparency', 0.5, selection_string)

    pymol.cmd.color('blue', selection_string)

    # Define the colors for each range
    colors = ['blue', 'cyan', 'white', 'yellow', 'orange', 'red']

    # Assuming cutoffs are sorted in increasing order and have length 4
    if len(cutoffs) != 5:
        raise ValueError("Cutoffs list must contain exactly four values.")

    for i, cutoff in enumerate(cutoffs):
        lower_bound = cutoffs[i - 1] if i > 0 else 0
        upper_bound = cutoff
        color_range_selection = f"({selection_string}) and b > {lower_bound} and b < {upper_bound}"
        pymol.cmd.color(colors[i], color_range_selection)

    # Color residues above the last cutoff
    last_cutoff_selection = f"({selection_string}) and b > {cutoffs[-1]}"
    pymol.cmd.color(colors[-1], last_cutoff_selection)

def show_CSP_ribbon(label, chain, cutoffs):
    selection_string = f"{label} and chain {chain}"

    # Set initial color for the ribbon
    pymol.cmd.color('white', selection_string)

    # Hide the surface and show the ribbon
    pymol.cmd.hide('surface', selection_string)
    #pymol.cmd.show('ribbon', selection_string)

    # Define the colors for each range
    colors = ['gray50', 'darksalmon', 'deepsalmon', 'firebrick']

    print(cutoffs)
    # Check if the cutoffs list contains exactly five values
    if len(cutoffs) != 3:
        raise ValueError("Cutoffs list must contain exactly five values.")

    color_range_selection = f"({selection_string}) and b < 0"
    pymol.cmd.color('white', color_range_selection)

    for i, cutoff in enumerate(cutoffs):
        lower_bound = cutoffs[i - 1] if i > 0 else 0
        upper_bound = cutoff
        color_range_selection = f"({selection_string}) and b > {lower_bound} and b < {upper_bound}"
        pymol.cmd.color(colors[i], color_range_selection)

    # Color residues above the last cutoff and show them as sticks
    last_cutoff_selection = f"({selection_string}) and b > {cutoffs[-1]}"
    pymol.cmd.color(colors[-1], last_cutoff_selection)
    pymol.cmd.show('sticks', last_cutoff_selection)

def show_atoms_hide_ribbon(label, chain):
    selection_string = f"{label} and chain {chain}"
    pymol.cmd.hide('everything', selection_string)
    #pymol.cmd.show('spheres', selection_string)
    pymol.cmd.show('sticks', selection_string)

def save_view_as_png(png_file_path, ray=1, dpi=300):
    # Enable ray tracing for higher quality
    if ray:
        pymol.cmd.ray()

    # Set transparent background
    pymol.cmd.set('ray_opaque_background', 0)
    
    # Save PNG with specified dpi

def process_png_files(bound, images_directory="images"):
    png_files = glob.glob(f"{images_directory}/{bound}*.png")
    for png in png_files:
        command = [
            "convert", png,
            "-background", "white",
            "-alpha", "remove",
            "-alpha", "off",
            png
        ]
        subprocess.run(command)

def process_structure(file_pref):
    pymol.cmd.color("red", file_pref + " and chain prot" + file_pref)
    if file_pref.find('NMR') != -1:
        pymol.cmd.color("yellow", file_pref + " and chain peptide" + file_pref)
    elif file_pref == "AF2":
        pymol.cmd.color("magenta", file_pref + " and chain peptide" + file_pref)
    elif file_pref == 'AF3':
        pymol.cmd.color("cyan", file_pref + " and chain peptide" + file_pref)
    elif file_pref == 'ES':
        pymol.cmd.color("blue", file_pref + " and chain peptide" + file_pref)
    else:
        print("received malformed file prefix: " + file_pref)
    pymol.cmd.color("green", file_pref + " and chain prot" + file_pref)

    pymol.cmd.orient()

    pymol.cmd.show('sticks', file_pref + ' and chain peptide' + file_pref)
    pymol.cmd.hide('cartoon', file_pref + ' and chain peptide' + file_pref)

def hide_residues(file_pref, well_defined_residues):
    pymol.cmd.hide('everything', file_pref + ' and chain prot' + file_pref + ' and resi -'+str(well_defined_res[0]))
    pymol.cmd.hide('everything', file_pref + ' and chain prot' + file_pref + ' and resi '+str(well_defined_res[1])+'-')

method = "MONTE"

CSmethod = "UCBShift"
while CSmethod not in ['SPARTA', 'UCBShift', 'ShiftX', 'consensus']:
    CSmethod = input("What CS prediciton method to use? [SPARTA, UCBShift, ShiftX, consensus]")

data_source_file = './data/CSP_UBQ.csv'

images_directory = './images/'

parsed_data = pd.read_csv(data_source_file)

apos = [str(data) for data in parsed_data['apo_bmrb']]
apo_pdbs = [str(data) for data in parsed_data['apo_pdb']]
holos = [data.lower() for data in parsed_data['holo_pdb']]
well_defined_residues = [data for data in parsed_data['Well_Defined_Residues']]
match_sequences = [data for data in parsed_data['match_seq']]

pdbs = input('Provide bound pdb ( e.g. "7jq8" ) ')
pdb = pdbs.strip()
holo = pdb.lower()

apo = str(apos[holos.index(holo.lower())]).lower()

# structures = input('what structures do you want to display [ AF2, AF3, NMR or ALL] ')
structures = 'ALL'
structures = [ structure.strip().replace(',','') for structure in structures.split(' ')]
if 'ALL' in structures:
    structures = ['AF2', 'AF3', 'NMR', 'ES']
z_scores = [0, 1, 3]

match_seq = match_sequences[apos.index(apo)]
well_defined_res = well_defined_residues[apos.index(apo)]
well_defined_res = [int(resi.strip())-1 for resi in well_defined_res.split(':')[1].split('..')]

new_files = []
cutoffs = []
structure = 'NMR_real'
CSPs, CSP_cutoff, bound_seq = calc_CSP_wrapper(apo, holo, well_defined_res, method=method, CSmethod='REAL', structure_source=structure, match_seq=match_seq)
CSP_below_thresh = [ C for C in CSPs if C < CSP_cutoff and C > 0 ]
for z_score in z_scores:
    cutoffs.append(calculate_z_score_threshold(CSP_below_thresh, z_score))

for structure in structures:
    bound_file = ""
    file_pref = structure 
    if structure == 'NMR':
        # load NMR w/ real shifts
        structure = 'NMR_real'
        file_pref = structure 
        bound_file = experimental_structures + 'exp_' + holo + '.pdb'
        genfiles = open_structures(bound_file, apo, file_pref, CSPs, bound_seq)
        for f in genfiles:
            new_files.append(f)
        process_structure(file_pref)
        show_CSP_ribbon(file_pref, 'prot' + file_pref, cutoffs)
        # hide_residues(file_pref, well_defined_res)
        # load NMR w/ pred shifts
        # structure = 'NMR_pred'
        # file_pref = structure 
        # bound_file = experimental_structures + 'exp_' + holo + '.pdb'
        # genfiles = open_structures(bound_file, apo, file_pref, CSPs, bound_seq)
        # for f in genfiles:
        #     new_files.append(f)
        # process_structure(file_pref)
        # show_CSP_ribbon(file_pref, 'prot' + file_pref, cutoffs)
        # hide_residues(file_pref, well_defined_res)
    elif structure == 'AF2':
        bound_file = computational_structures + 'comp_' + holo + '.pdb'
        genfiles = open_structures(bound_file, apo, file_pref, CSPs, bound_seq)
        for f in genfiles:
            new_files.append(f)
        process_structure(file_pref)
        show_CSP_ribbon(file_pref, 'prot' + file_pref, cutoffs)
    elif structure == 'ES':
        bound_dir = f'./PDB_FILES/{holo.lower()}_max_RPF_NLDR_CSPRank_files/'
        files = [ bound_dir + f for f in listdir(bound_dir) if f.endswith('.pdb') ]
        medoid_file = find_medoid_structure(files)
        genfiles = open_structures(medoid_file, apo, file_pref, CSPs, bound_seq)
        for f in genfiles:
            new_files.append(f)
        process_structure(file_pref)
        show_CSP_ribbon(file_pref, 'prot' + file_pref, cutoffs)
    else:
        print("received malformed structure selection: " + structure)
        continue

pymol.cmd.align('AF2 and chain protAF2', 'chain protNMR_real')
pymol.cmd.align('ES and chain protES', 'chain protNMR_real')
# pymol.cmd.align('AF3 and chain protAF3', 'chain protNMR_real')
pymol.cmd.hide('everything', 'hydrogen')
pymol.cmd.orient()

for file in new_files:
    os.system('rm ' + file)

