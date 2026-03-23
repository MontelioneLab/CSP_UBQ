reinitialize
load ./outputs/2YKA/2YKA_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 83
color red, structure and chain A and resi 85
color red, structure and chain A and resi 86
color red, structure and chain A and resi 87
color red, structure and chain A and resi 90
color red, structure and chain A and resi 91
color red, structure and chain A and resi 94
color red, structure and chain A and resi 97
color red, structure and chain A and resi 134
color red, structure and chain A and resi 135
color red, structure and chain A and resi 137
color red, structure and chain A and resi 138
color red, structure and chain A and resi 139
color red, structure and chain A and resi 140
color red, structure and chain A and resi 141
color red, structure and chain A and resi 145
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [83, 85, 86, 87, 90, 91, 94, 97, 134, 135, 137, 138, 139, 140, 141, 145]
