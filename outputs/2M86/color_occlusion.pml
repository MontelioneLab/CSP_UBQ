reinitialize
load ./outputs/2M86/2M86_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 115
color red, structure and chain B and resi 116
color red, structure and chain B and resi 119
color red, structure and chain B and resi 120
color red, structure and chain B and resi 123
color red, structure and chain B and resi 127
color red, structure and chain B and resi 137
color red, structure and chain B and resi 138
color red, structure and chain B and resi 159
color red, structure and chain B and resi 161
color red, structure and chain B and resi 164
color red, structure and chain B and resi 165
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [115, 116, 119, 120, 123, 127, 137, 138, 159, 161, 164, 165]
