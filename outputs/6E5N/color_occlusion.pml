reinitialize
load ./outputs/6E5N/6E5N_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 1050
color red, structure and chain B and resi 1051
color red, structure and chain B and resi 1053
color red, structure and chain B and resi 1054
color red, structure and chain B and resi 1055
color red, structure and chain B and resi 1058
color red, structure and chain B and resi 1062
color red, structure and chain B and resi 1087
color red, structure and chain B and resi 1088
color red, structure and chain B and resi 1090
color red, structure and chain B and resi 1091
color red, structure and chain B and resi 1117
color red, structure and chain B and resi 1120
color red, structure and chain B and resi 1121
color red, structure and chain B and resi 1124
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1050, 1051, 1053, 1054, 1055, 1058, 1062, 1087, 1088, 1090, 1091, 1117, 1120, 1121, 1124]
