reinitialize
load ./outputs/2RUI/2RUI_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 103
color red, structure and chain A and resi 110
color red, structure and chain A and resi 111
color red, structure and chain A and resi 124
color red, structure and chain A and resi 125
color red, structure and chain A and resi 126
color red, structure and chain A and resi 168
color red, structure and chain A and resi 169
color red, structure and chain A and resi 170
color red, structure and chain A and resi 171
color red, structure and chain A and resi 185
color red, structure and chain A and resi 186
color red, structure and chain A and resi 187
color red, structure and chain A and resi 196
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [59, 60, 103, 110, 111, 124, 125, 126, 168, 169, 170, 171, 185, 186, 187, 196]
