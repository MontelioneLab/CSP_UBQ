reinitialize
load ./outputs/2L29_2/2L29_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 6
color red, structure and chain B and resi 11
color red, structure and chain B and resi 15
color red, structure and chain B and resi 18
color red, structure and chain B and resi 19
color red, structure and chain B and resi 22
color red, structure and chain B and resi 23
color red, structure and chain B and resi 50
color red, structure and chain B and resi 52
color red, structure and chain B and resi 53
color red, structure and chain B and resi 57
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 11, 15, 18, 19, 22, 23, 50, 52, 53, 57]
