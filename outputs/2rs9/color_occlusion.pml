reinitialize
load ./outputs/2rs9/2rs9_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 77
color red, structure and chain B and resi 119
color red, structure and chain B and resi 122
color red, structure and chain B and resi 123
color red, structure and chain B and resi 124
color red, structure and chain B and resi 126
color red, structure and chain B and resi 127
color red, structure and chain B and resi 129
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [77, 119, 122, 123, 124, 126, 127, 129]
