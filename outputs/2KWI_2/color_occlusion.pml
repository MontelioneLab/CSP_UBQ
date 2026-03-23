reinitialize
load ./outputs/2KWI_2/2KWI_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 409
color red, structure and chain B and resi 413
color red, structure and chain B and resi 416
color red, structure and chain B and resi 417
color red, structure and chain B and resi 421
color red, structure and chain B and resi 422
color red, structure and chain B and resi 426
color red, structure and chain B and resi 427
color red, structure and chain B and resi 429
color red, structure and chain B and resi 430
color red, structure and chain B and resi 433
color red, structure and chain B and resi 434
color red, structure and chain B and resi 437
color red, structure and chain B and resi 440
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [409, 413, 416, 417, 421, 422, 426, 427, 429, 430, 433, 434, 437, 440]
