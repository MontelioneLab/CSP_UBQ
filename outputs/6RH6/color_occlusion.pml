reinitialize
load ./outputs/6RH6/6RH6_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 38
color red, structure and chain A and resi 85
color red, structure and chain A and resi 86
color red, structure and chain A and resi 87
color red, structure and chain A and resi 88
color red, structure and chain A and resi 91
color red, structure and chain A and resi 113
color red, structure and chain A and resi 114
color red, structure and chain A and resi 117
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [33, 34, 38, 85, 86, 87, 88, 91, 113, 114, 117]
