reinitialize
load ./outputs/1DDM/1DDM_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 144
color red, structure and chain A and resi 145
color red, structure and chain A and resi 146
color red, structure and chain A and resi 147
color red, structure and chain A and resi 149
color red, structure and chain A and resi 150
color red, structure and chain A and resi 151
color red, structure and chain A and resi 152
color red, structure and chain A and resi 153
color red, structure and chain A and resi 188
color red, structure and chain A and resi 192
color red, structure and chain A and resi 195
color red, structure and chain A and resi 198
color red, structure and chain A and resi 199
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [144, 145, 146, 147, 149, 150, 151, 152, 153, 188, 192, 195, 198, 199]
