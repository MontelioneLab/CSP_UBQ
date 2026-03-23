reinitialize
load ./outputs/2KNH/2KNH_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 273
color red, structure and chain A and resi 276
color red, structure and chain A and resi 277
color red, structure and chain A and resi 280
color red, structure and chain A and resi 284
color red, structure and chain A and resi 325
color red, structure and chain A and resi 326
color red, structure and chain A and resi 329
color red, structure and chain A and resi 332
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [273, 276, 277, 280, 284, 325, 326, 329, 332]
