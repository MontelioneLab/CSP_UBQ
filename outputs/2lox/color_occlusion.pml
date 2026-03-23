reinitialize
load ./outputs/2lox/2lox_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 11
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 55
color red, structure and chain A and resi 57
color red, structure and chain A and resi 59
color red, structure and chain A and resi 61
color red, structure and chain A and resi 88
color red, structure and chain A and resi 101
color red, structure and chain A and resi 105
color red, structure and chain A and resi 108
color red, structure and chain A and resi 109
color red, structure and chain A and resi 112
color red, structure and chain A and resi 113
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [11, 47, 48, 49, 50, 51, 52, 55, 57, 59, 61, 88, 101, 105, 108, 109, 112, 113]
