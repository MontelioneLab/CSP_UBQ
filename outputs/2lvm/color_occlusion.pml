reinitialize
load ./outputs/2lvm/2lvm_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 1495
color red, structure and chain A and resi 1498
color red, structure and chain A and resi 1500
color red, structure and chain A and resi 1502
color red, structure and chain A and resi 1519
color red, structure and chain A and resi 1521
color red, structure and chain A and resi 1523
color red, structure and chain A and resi 1545
color red, structure and chain A and resi 1547
color red, structure and chain A and resi 1550
color red, structure and chain A and resi 1551
color red, structure and chain A and resi 1552
color red, structure and chain A and resi 1553
color red, structure and chain A and resi 1587
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1495, 1498, 1500, 1502, 1519, 1521, 1523, 1545, 1547, 1550, 1551, 1552, 1553, 1587]
