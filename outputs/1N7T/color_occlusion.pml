reinitialize
load ./outputs/1N7T/1N7T_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 23
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 37
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 79
color red, structure and chain A and resi 83
color red, structure and chain A and resi 87
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [23, 25, 26, 27, 28, 34, 35, 37, 48, 49, 50, 51, 79, 83, 87]
