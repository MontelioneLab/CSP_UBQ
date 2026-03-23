reinitialize
load ./outputs/2RQG/2RQG_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 559
color red, structure and chain B and resi 560
color red, structure and chain B and resi 564
color red, structure and chain B and resi 583
color red, structure and chain B and resi 584
color red, structure and chain B and resi 586
color red, structure and chain B and resi 587
color red, structure and chain B and resi 588
color red, structure and chain B and resi 592
color red, structure and chain B and resi 606
color red, structure and chain B and resi 609
color red, structure and chain B and resi 613
color red, structure and chain B and resi 614
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [559, 560, 564, 583, 584, 586, 587, 588, 592, 606, 609, 613, 614]
