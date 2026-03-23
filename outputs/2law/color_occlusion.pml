reinitialize
load ./outputs/2law/2law_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 237
color red, structure and chain A and resi 239
color red, structure and chain A and resi 241
color red, structure and chain A and resi 247
color red, structure and chain A and resi 249
color red, structure and chain A and resi 251
color red, structure and chain A and resi 254
color red, structure and chain A and resi 255
color red, structure and chain A and resi 256
color red, structure and chain A and resi 258
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [237, 239, 241, 247, 249, 251, 254, 255, 256, 258]
