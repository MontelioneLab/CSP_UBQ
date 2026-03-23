reinitialize
load ./outputs/2n1g/2n1g_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 1160
color red, structure and chain A and resi 1161
color red, structure and chain A and resi 1162
color red, structure and chain A and resi 1167
color red, structure and chain A and resi 1170
color red, structure and chain A and resi 1171
color red, structure and chain A and resi 1172
color red, structure and chain A and resi 1174
color red, structure and chain A and resi 1175
color red, structure and chain A and resi 1179
color red, structure and chain A and resi 1183
color red, structure and chain A and resi 1185
color red, structure and chain A and resi 1186
color red, structure and chain A and resi 1189
color red, structure and chain A and resi 1190
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1160, 1161, 1162, 1167, 1170, 1171, 1172, 1174, 1175, 1179, 1183, 1185, 1186, 1189, 1190]
