reinitialize
load ./outputs/2FIN_1/2FIN_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 53
color red, structure and chain A and resi 75
color red, structure and chain A and resi 78
color red, structure and chain A and resi 80
color red, structure and chain A and resi 141
color red, structure and chain A and resi 143
color red, structure and chain A and resi 147
color red, structure and chain A and resi 181
color red, structure and chain A and resi 182
color red, structure and chain A and resi 183
color red, structure and chain A and resi 184
color red, structure and chain A and resi 185
color red, structure and chain A and resi 186
color red, structure and chain A and resi 187
color red, structure and chain A and resi 188
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [53, 75, 78, 80, 141, 143, 147, 181, 182, 183, 184, 185, 186, 187, 188]
