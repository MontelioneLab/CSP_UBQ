reinitialize
load ./outputs/6FDP/6FDP_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 286
color red, structure and chain A and resi 290
color red, structure and chain A and resi 293
color red, structure and chain A and resi 294
color red, structure and chain A and resi 317
color red, structure and chain A and resi 321
color red, structure and chain A and resi 324
color red, structure and chain A and resi 327
color red, structure and chain A and resi 328
color red, structure and chain A and resi 350
color red, structure and chain A and resi 351
color red, structure and chain A and resi 354
color red, structure and chain A and resi 355
color red, structure and chain A and resi 358
color red, structure and chain A and resi 380
color red, structure and chain A and resi 383
color red, structure and chain A and resi 384
color red, structure and chain A and resi 385
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [286, 290, 293, 294, 317, 321, 324, 327, 328, 350, 351, 354, 355, 358, 380, 383, 384, 385]
