reinitialize
load ./outputs/2nmb/2nmb_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 149
color red, structure and chain A and resi 150
color red, structure and chain A and resi 151
color red, structure and chain A and resi 160
color red, structure and chain A and resi 199
color red, structure and chain A and resi 202
color red, structure and chain A and resi 204
color red, structure and chain A and resi 206
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [149, 150, 151, 160, 199, 202, 204, 206]
