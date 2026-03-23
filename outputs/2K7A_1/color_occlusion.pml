reinitialize
load ./outputs/2K7A_1/2K7A_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 280
color red, structure and chain B and resi 281
color red, structure and chain B and resi 282
color red, structure and chain B and resi 283
color red, structure and chain B and resi 290
color red, structure and chain B and resi 292
color red, structure and chain B and resi 295
color red, structure and chain B and resi 307
color red, structure and chain B and resi 309
color red, structure and chain B and resi 328
color red, structure and chain B and resi 329
color red, structure and chain B and resi 330
color red, structure and chain B and resi 331
color red, structure and chain B and resi 332
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [280, 281, 282, 283, 290, 292, 295, 307, 309, 328, 329, 330, 331, 332]
