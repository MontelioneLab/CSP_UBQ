reinitialize
load ./outputs/2B0F/2B0F_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 22
color red, structure and chain A and resi 40
color red, structure and chain A and resi 124
color red, structure and chain A and resi 125
color red, structure and chain A and resi 126
color red, structure and chain A and resi 127
color red, structure and chain A and resi 128
color red, structure and chain A and resi 141
color red, structure and chain A and resi 142
color red, structure and chain A and resi 143
color red, structure and chain A and resi 144
color red, structure and chain A and resi 146
color red, structure and chain A and resi 160
color red, structure and chain A and resi 161
color red, structure and chain A and resi 162
color red, structure and chain A and resi 163
color red, structure and chain A and resi 164
color red, structure and chain A and resi 169
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [22, 40, 124, 125, 126, 127, 128, 141, 142, 143, 144, 146, 160, 161, 162, 163, 164, 169]
