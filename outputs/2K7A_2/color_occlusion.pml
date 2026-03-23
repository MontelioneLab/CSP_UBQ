reinitialize
load ./outputs/2K7A_2/2K7A_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 179
color red, structure and chain A and resi 180
color red, structure and chain A and resi 182
color red, structure and chain A and resi 185
color red, structure and chain A and resi 189
color red, structure and chain A and resi 205
color red, structure and chain A and resi 206
color red, structure and chain A and resi 207
color red, structure and chain A and resi 208
color red, structure and chain A and resi 222
color red, structure and chain A and resi 223
color red, structure and chain A and resi 224
color red, structure and chain A and resi 225
color red, structure and chain A and resi 227
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [179, 180, 182, 185, 189, 205, 206, 207, 208, 222, 223, 224, 225, 227]
