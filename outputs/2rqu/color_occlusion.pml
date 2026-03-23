reinitialize
load ./outputs/2rqu/2rqu_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 1075
color red, structure and chain A and resi 1076
color red, structure and chain A and resi 1077
color red, structure and chain A and resi 1079
color red, structure and chain A and resi 1081
color red, structure and chain A and resi 1082
color red, structure and chain A and resi 1084
color red, structure and chain A and resi 1085
color red, structure and chain A and resi 1090
color red, structure and chain A and resi 1100
color red, structure and chain A and resi 1101
color red, structure and chain A and resi 1103
color red, structure and chain A and resi 1104
color red, structure and chain A and resi 1118
color red, structure and chain A and resi 1122
color red, structure and chain A and resi 1123
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1075, 1076, 1077, 1079, 1081, 1082, 1084, 1085, 1090, 1100, 1101, 1103, 1104, 1118, 1122, 1123]
