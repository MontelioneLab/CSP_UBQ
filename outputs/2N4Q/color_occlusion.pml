reinitialize
load ./outputs/2N4Q/2N4Q_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 504
color red, structure and chain B and resi 508
color red, structure and chain B and resi 511
color red, structure and chain B and resi 524
color red, structure and chain B and resi 527
color red, structure and chain B and resi 542
color red, structure and chain B and resi 543
color red, structure and chain B and resi 544
color red, structure and chain B and resi 545
color red, structure and chain B and resi 546
color red, structure and chain B and resi 547
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [504, 508, 511, 524, 527, 542, 543, 544, 545, 546, 547]
