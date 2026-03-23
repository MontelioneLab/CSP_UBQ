reinitialize
load ./outputs/2L29_1/2L29_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 1543
color red, structure and chain A and resi 1544
color red, structure and chain A and resi 1567
color red, structure and chain A and resi 1569
color red, structure and chain A and resi 1570
color red, structure and chain A and resi 1572
color red, structure and chain A and resi 1597
color red, structure and chain A and resi 1600
color red, structure and chain A and resi 1601
color red, structure and chain A and resi 1602
color red, structure and chain A and resi 1629
color red, structure and chain A and resi 1631
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1543, 1544, 1567, 1569, 1570, 1572, 1597, 1600, 1601, 1602, 1629, 1631]
