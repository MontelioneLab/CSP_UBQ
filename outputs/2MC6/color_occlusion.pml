reinitialize
load ./outputs/2MC6/2MC6_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 2
color red, structure and chain A and resi 5
color red, structure and chain A and resi 6
color red, structure and chain A and resi 7
color red, structure and chain A and resi 8
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 19
color red, structure and chain A and resi 41
color red, structure and chain A and resi 47
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 54
color red, structure and chain A and resi 67
color red, structure and chain A and resi 68
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [2, 5, 6, 7, 8, 15, 16, 17, 19, 41, 47, 50, 51, 54, 67, 68]
