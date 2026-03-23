reinitialize
load ./outputs/2hug/2hug_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 7
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 10
color red, structure and chain A and resi 21
color red, structure and chain A and resi 27
color red, structure and chain A and resi 29
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 42
color red, structure and chain A and resi 43
color red, structure and chain A and resi 44
color red, structure and chain A and resi 46
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 7, 8, 9, 10, 21, 27, 29, 35, 36, 37, 38, 42, 43, 44, 46]
