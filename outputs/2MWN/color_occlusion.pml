reinitialize
load ./outputs/2MWN/2MWN_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 316
color red, structure and chain B and resi 323
color red, structure and chain B and resi 324
color red, structure and chain B and resi 325
color red, structure and chain B and resi 358
color red, structure and chain B and resi 360
color red, structure and chain B and resi 362
color red, structure and chain B and resi 365
color red, structure and chain B and resi 367
color red, structure and chain B and resi 369
color red, structure and chain B and resi 379
color red, structure and chain B and resi 380
color red, structure and chain B and resi 381
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [316, 323, 324, 325, 358, 360, 362, 365, 367, 369, 379, 380, 381]
