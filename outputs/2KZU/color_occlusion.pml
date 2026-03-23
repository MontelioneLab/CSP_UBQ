reinitialize
load ./outputs/2KZU/2KZU_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 83
color red, structure and chain A and resi 84
color red, structure and chain A and resi 87
color red, structure and chain A and resi 91
color red, structure and chain A and resi 121
color red, structure and chain A and resi 122
color red, structure and chain A and resi 124
color red, structure and chain A and resi 125
color red, structure and chain A and resi 127
color red, structure and chain A and resi 128
color red, structure and chain A and resi 132
color red, structure and chain A and resi 135
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [83, 84, 87, 91, 121, 122, 124, 125, 127, 128, 132, 135]
