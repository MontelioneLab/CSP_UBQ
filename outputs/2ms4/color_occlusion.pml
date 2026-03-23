reinitialize
load ./outputs/2ms4/2ms4_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 55
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 63
color red, structure and chain A and resi 101
color red, structure and chain A and resi 102
color red, structure and chain A and resi 103
color red, structure and chain A and resi 110
color red, structure and chain A and resi 111
color red, structure and chain A and resi 121
color red, structure and chain A and resi 122
color red, structure and chain A and resi 148
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [55, 57, 58, 59, 60, 61, 63, 101, 102, 103, 110, 111, 121, 122, 148]
