reinitialize
load ./outputs/5I22/5I22_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 529
color red, structure and chain A and resi 535
color red, structure and chain A and resi 537
color red, structure and chain A and resi 538
color red, structure and chain A and resi 556
color red, structure and chain A and resi 557
color red, structure and chain A and resi 559
color red, structure and chain A and resi 560
color red, structure and chain A and resi 562
color red, structure and chain A and resi 583
color red, structure and chain A and resi 585
color red, structure and chain A and resi 587
color red, structure and chain A and resi 588
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [529, 535, 537, 538, 556, 557, 559, 560, 562, 583, 585, 587, 588]
