reinitialize
load ./outputs/2jw1/2jw1_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 31
color red, structure and chain A and resi 34
color red, structure and chain A and resi 36
color red, structure and chain A and resi 56
color red, structure and chain A and resi 58
color red, structure and chain A and resi 61
color red, structure and chain A and resi 78
color red, structure and chain A and resi 81
color red, structure and chain A and resi 117
color red, structure and chain A and resi 118
color red, structure and chain A and resi 119
color red, structure and chain A and resi 122
color red, structure and chain A and resi 123
color red, structure and chain A and resi 125
color red, structure and chain A and resi 126
color red, structure and chain A and resi 127
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [28, 29, 31, 34, 36, 56, 58, 61, 78, 81, 117, 118, 119, 122, 123, 125, 126, 127]
