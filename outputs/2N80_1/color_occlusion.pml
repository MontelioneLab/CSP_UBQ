reinitialize
load ./outputs/2N80_1/2N80_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 342
color red, structure and chain A and resi 343
color red, structure and chain A and resi 346
color red, structure and chain A and resi 349
color red, structure and chain A and resi 350
color red, structure and chain A and resi 352
color red, structure and chain A and resi 353
color red, structure and chain A and resi 354
color red, structure and chain A and resi 355
color red, structure and chain A and resi 356
color red, structure and chain A and resi 409
color red, structure and chain A and resi 410
color red, structure and chain A and resi 412
color red, structure and chain A and resi 413
color red, structure and chain A and resi 416
color red, structure and chain A and resi 419
color red, structure and chain A and resi 420
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [342, 343, 346, 349, 350, 352, 353, 354, 355, 356, 409, 410, 412, 413, 416, 419, 420]
