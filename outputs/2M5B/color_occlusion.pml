reinitialize
load ./outputs/2M5B/2M5B_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 81
color red, structure and chain A and resi 85
color red, structure and chain A and resi 89
color red, structure and chain A and resi 92
color red, structure and chain A and resi 93
color red, structure and chain A and resi 96
color red, structure and chain A and resi 99
color red, structure and chain A and resi 100
color red, structure and chain A and resi 113
color red, structure and chain A and resi 114
color red, structure and chain A and resi 117
color red, structure and chain A and resi 118
color red, structure and chain A and resi 126
color red, structure and chain A and resi 127
color red, structure and chain A and resi 129
color red, structure and chain A and resi 130
color red, structure and chain A and resi 183
color red, structure and chain A and resi 184
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [81, 85, 89, 92, 93, 96, 99, 100, 113, 114, 117, 118, 126, 127, 129, 130, 183, 184]
