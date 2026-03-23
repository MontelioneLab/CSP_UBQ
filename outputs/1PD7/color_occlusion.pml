reinitialize
load ./outputs/1PD7/1PD7_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 7
color red, structure and chain A and resi 10
color red, structure and chain A and resi 11
color red, structure and chain A and resi 13
color red, structure and chain A and resi 14
color red, structure and chain A and resi 15
color red, structure and chain A and resi 29
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 58
color red, structure and chain A and resi 76
color red, structure and chain A and resi 78
color red, structure and chain A and resi 79
color red, structure and chain A and resi 80
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 19
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 7, 10, 11, 13, 14, 15, 29, 31, 32, 35, 36, 38, 39, 58, 76, 78, 79, 80]
