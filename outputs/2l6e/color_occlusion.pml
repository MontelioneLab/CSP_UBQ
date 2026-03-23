reinitialize
load ./outputs/2l6e/2l6e_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 165
color red, structure and chain A and resi 166
color red, structure and chain A and resi 169
color red, structure and chain A and resi 172
color red, structure and chain A and resi 173
color red, structure and chain A and resi 179
color red, structure and chain A and resi 183
color red, structure and chain A and resi 186
color red, structure and chain A and resi 187
color red, structure and chain A and resi 190
color red, structure and chain A and resi 209
color red, structure and chain A and resi 210
color red, structure and chain A and resi 211
color red, structure and chain A and resi 212
color red, structure and chain A and resi 214
color red, structure and chain A and resi 215
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [165, 166, 169, 172, 173, 179, 183, 186, 187, 190, 209, 210, 211, 212, 214, 215]
