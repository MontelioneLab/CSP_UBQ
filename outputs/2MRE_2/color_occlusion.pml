reinitialize
load ./outputs/2MRE_2/2MRE_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 202
color red, structure and chain B and resi 203
color red, structure and chain B and resi 204
color red, structure and chain B and resi 205
color red, structure and chain B and resi 206
color red, structure and chain B and resi 207
color red, structure and chain B and resi 213
color red, structure and chain B and resi 216
color red, structure and chain B and resi 217
color red, structure and chain B and resi 219
color red, structure and chain B and resi 220
color red, structure and chain B and resi 221
color red, structure and chain B and resi 223
color red, structure and chain B and resi 224
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [202, 203, 204, 205, 206, 207, 213, 216, 217, 219, 220, 221, 223, 224]
