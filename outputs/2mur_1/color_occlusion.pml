reinitialize
load ./outputs/2mur_1/2mur_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 108
color red, structure and chain B and resi 142
color red, structure and chain B and resi 144
color red, structure and chain B and resi 145
color red, structure and chain B and resi 146
color red, structure and chain B and resi 147
color red, structure and chain B and resi 148
color red, structure and chain B and resi 149
color red, structure and chain B and resi 166
color red, structure and chain B and resi 168
color red, structure and chain B and resi 170
color red, structure and chain B and resi 171
color red, structure and chain B and resi 172
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [108, 142, 144, 145, 146, 147, 148, 149, 166, 168, 170, 171, 172]
