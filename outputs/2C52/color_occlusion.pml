reinitialize
load ./outputs/2C52/2C52_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 5
color red, structure and chain A and resi 6
color red, structure and chain A and resi 8
color red, structure and chain A and resi 11
color red, structure and chain A and resi 14
color red, structure and chain A and resi 15
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 22
color red, structure and chain A and resi 43
color red, structure and chain A and resi 47
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [5, 6, 8, 11, 14, 15, 18, 19, 22, 43, 47, 51, 52]
