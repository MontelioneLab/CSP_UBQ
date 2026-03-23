reinitialize
load ./outputs/2kfh/2kfh_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 69
color red, structure and chain A and resi 73
color red, structure and chain A and resi 76
color red, structure and chain A and resi 83
color red, structure and chain A and resi 86
color red, structure and chain A and resi 87
color red, structure and chain A and resi 90
color red, structure and chain A and resi 94
color red, structure and chain A and resi 96
color red, structure and chain A and resi 98
color red, structure and chain A and resi 100
color red, structure and chain A and resi 105
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [69, 73, 76, 83, 86, 87, 90, 94, 96, 98, 100, 105]
