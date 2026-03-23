reinitialize
load ./outputs/6ijq/6ijq_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 52
color red, structure and chain B and resi 53
color red, structure and chain B and resi 54
color red, structure and chain B and resi 56
color red, structure and chain B and resi 57
color red, structure and chain B and resi 61
color red, structure and chain B and resi 98
color red, structure and chain B and resi 99
color red, structure and chain B and resi 151
color red, structure and chain B and resi 154
color red, structure and chain B and resi 155
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [52, 53, 54, 56, 57, 61, 98, 99, 151, 154, 155]
