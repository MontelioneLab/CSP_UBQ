reinitialize
load ./outputs/2FFK/2FFK_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 52
color red, structure and chain A and resi 53
color red, structure and chain A and resi 64
color red, structure and chain A and resi 65
color red, structure and chain A and resi 78
color red, structure and chain A and resi 80
color red, structure and chain A and resi 93
color red, structure and chain A and resi 96
color red, structure and chain A and resi 141
color red, structure and chain A and resi 143
color red, structure and chain A and resi 145
color red, structure and chain A and resi 147
color red, structure and chain A and resi 181
color red, structure and chain A and resi 182
color red, structure and chain A and resi 183
color red, structure and chain A and resi 184
color red, structure and chain A and resi 185
color red, structure and chain A and resi 187
color red, structure and chain A and resi 217
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 19
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [52, 53, 64, 65, 78, 80, 93, 96, 141, 143, 145, 147, 181, 182, 183, 184, 185, 187, 217]
