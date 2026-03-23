reinitialize
load ./outputs/6OQJ_2/6OQJ_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 102
color red, structure and chain B and resi 117
color red, structure and chain B and resi 121
color red, structure and chain B and resi 126
color red, structure and chain B and resi 128
color red, structure and chain B and resi 129
color red, structure and chain B and resi 130
color red, structure and chain B and resi 132
color red, structure and chain B and resi 133
color red, structure and chain B and resi 136
color red, structure and chain B and resi 139
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [102, 117, 121, 126, 128, 129, 130, 132, 133, 136, 139]
