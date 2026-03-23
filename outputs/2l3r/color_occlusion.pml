reinitialize
load ./outputs/2l3r/2l3r_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 142
color red, structure and chain A and resi 145
color red, structure and chain A and resi 147
color red, structure and chain A and resi 148
color red, structure and chain A and resi 152
color red, structure and chain A and resi 153
color red, structure and chain A and resi 188
color red, structure and chain A and resi 190
color red, structure and chain A and resi 191
color red, structure and chain A and resi 208
color red, structure and chain A and resi 211
color red, structure and chain A and resi 235
color red, structure and chain A and resi 236
color red, structure and chain A and resi 237
color red, structure and chain A and resi 238
color red, structure and chain A and resi 275
color red, structure and chain A and resi 276
color red, structure and chain A and resi 278
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [142, 145, 147, 148, 152, 153, 188, 190, 191, 208, 211, 235, 236, 237, 238, 275, 276, 278]
