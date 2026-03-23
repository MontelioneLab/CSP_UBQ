reinitialize
load ./outputs/2MV7/2MV7_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 507
color red, structure and chain A and resi 511
color red, structure and chain A and resi 520
color red, structure and chain A and resi 523
color red, structure and chain A and resi 524
color red, structure and chain A and resi 527
color red, structure and chain A and resi 528
color red, structure and chain A and resi 531
color red, structure and chain A and resi 538
color red, structure and chain A and resi 541
color red, structure and chain A and resi 542
color red, structure and chain A and resi 543
color red, structure and chain A and resi 544
color red, structure and chain A and resi 545
color red, structure and chain A and resi 546
color red, structure and chain A and resi 547
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [507, 511, 520, 523, 524, 527, 528, 531, 538, 541, 542, 543, 544, 545, 546, 547]
