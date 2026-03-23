reinitialize
load ./outputs/2MOW/2MOW_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 21
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 74
color red, structure and chain A and resi 127
color red, structure and chain A and resi 130
color red, structure and chain A and resi 133
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [21, 25, 26, 27, 28, 29, 30, 70, 71, 74, 127, 130, 133]
