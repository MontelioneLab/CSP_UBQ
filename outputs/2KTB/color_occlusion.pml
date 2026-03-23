reinitialize
load ./outputs/2KTB/2KTB_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 198
color red, structure and chain B and resi 202
color red, structure and chain B and resi 203
color red, structure and chain B and resi 206
color red, structure and chain B and resi 207
color red, structure and chain B and resi 262
color red, structure and chain B and resi 263
color red, structure and chain B and resi 266
color red, structure and chain B and resi 267
color red, structure and chain B and resi 268
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [198, 202, 203, 206, 207, 262, 263, 266, 267, 268]
