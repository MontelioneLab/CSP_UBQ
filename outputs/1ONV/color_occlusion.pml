reinitialize
load ./outputs/1ONV/1ONV_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 470
color red, structure and chain A and resi 471
color red, structure and chain A and resi 474
color red, structure and chain A and resi 478
color red, structure and chain A and resi 486
color red, structure and chain A and resi 490
color red, structure and chain A and resi 493
color red, structure and chain A and resi 494
color red, structure and chain A and resi 497
color red, structure and chain A and resi 498
color red, structure and chain A and resi 501
color red, structure and chain A and resi 502
color red, structure and chain A and resi 504
color red, structure and chain A and resi 513
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [470, 471, 474, 478, 486, 490, 493, 494, 497, 498, 501, 502, 504, 513]
