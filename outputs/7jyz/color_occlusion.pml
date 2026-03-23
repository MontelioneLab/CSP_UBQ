reinitialize
load ./outputs/7jyz/7jyz_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 230
color red, structure and chain B and resi 233
color red, structure and chain B and resi 234
color red, structure and chain B and resi 237
color red, structure and chain B and resi 247
color red, structure and chain B and resi 248
color red, structure and chain B and resi 251
color red, structure and chain B and resi 267
color red, structure and chain B and resi 268
color red, structure and chain B and resi 269
color red, structure and chain B and resi 270
color red, structure and chain B and resi 271
color red, structure and chain B and resi 272
color red, structure and chain B and resi 273
color red, structure and chain B and resi 274
color red, structure and chain B and resi 275
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [230, 233, 234, 237, 247, 248, 251, 267, 268, 269, 270, 271, 272, 273, 274, 275]
