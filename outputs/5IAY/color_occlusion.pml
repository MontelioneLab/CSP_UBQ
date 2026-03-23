reinitialize
load ./outputs/5IAY/5IAY_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 142
color red, structure and chain A and resi 153
color red, structure and chain A and resi 189
color red, structure and chain A and resi 190
color red, structure and chain A and resi 207
color red, structure and chain A and resi 228
color red, structure and chain A and resi 234
color red, structure and chain A and resi 235
color red, structure and chain A and resi 236
color red, structure and chain A and resi 237
color red, structure and chain A and resi 238
color red, structure and chain A and resi 278
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [142, 153, 189, 190, 207, 228, 234, 235, 236, 237, 238, 278]
