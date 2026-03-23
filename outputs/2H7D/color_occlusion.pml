reinitialize
load ./outputs/2H7D/2H7D_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 318
color red, structure and chain A and resi 324
color red, structure and chain A and resi 325
color red, structure and chain A and resi 354
color red, structure and chain A and resi 357
color red, structure and chain A and resi 358
color red, structure and chain A and resi 359
color red, structure and chain A and resi 360
color red, structure and chain A and resi 362
color red, structure and chain A and resi 363
color red, structure and chain A and resi 364
color red, structure and chain A and resi 368
color red, structure and chain A and resi 369
color red, structure and chain A and resi 379
color red, structure and chain A and resi 380
color red, structure and chain A and resi 381
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [318, 324, 325, 354, 357, 358, 359, 360, 362, 363, 364, 368, 369, 379, 380, 381]
