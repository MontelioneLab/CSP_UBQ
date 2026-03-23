reinitialize
load ./outputs/2MCN_2/2MCN_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 10
color red, structure and chain A and resi 11
color red, structure and chain A and resi 13
color red, structure and chain A and resi 14
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 32
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 53
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [8, 9, 10, 11, 13, 14, 16, 17, 32, 33, 34, 35, 36, 37, 51, 52, 53]
