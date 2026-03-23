reinitialize
load ./outputs/2mnu/2mnu_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 11
color red, structure and chain A and resi 12
color red, structure and chain A and resi 14
color red, structure and chain A and resi 17
color red, structure and chain A and resi 19
color red, structure and chain A and resi 22
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 72
color red, structure and chain A and resi 73
color red, structure and chain A and resi 74
color red, structure and chain A and resi 75
color red, structure and chain A and resi 76
color red, structure and chain A and resi 77
color red, structure and chain A and resi 78
color red, structure and chain A and resi 79
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [11, 12, 14, 17, 19, 22, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79]
