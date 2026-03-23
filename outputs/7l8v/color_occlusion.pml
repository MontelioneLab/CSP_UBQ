reinitialize
load ./outputs/7l8v/7l8v_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 88
color red, structure and chain A and resi 89
color red, structure and chain A and resi 94
color red, structure and chain A and resi 96
color red, structure and chain A and resi 98
color red, structure and chain A and resi 100
color red, structure and chain A and resi 105
color red, structure and chain A and resi 110
color red, structure and chain A and resi 125
color red, structure and chain A and resi 130
color red, structure and chain A and resi 132
color red, structure and chain A and resi 134
color red, structure and chain A and resi 136
color red, structure and chain A and resi 141
color red, structure and chain A and resi 145
color red, structure and chain A and resi 146
color cyan, structure and chain C
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: C
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [88, 89, 94, 96, 98, 100, 105, 110, 125, 130, 132, 134, 136, 141, 145, 146]
