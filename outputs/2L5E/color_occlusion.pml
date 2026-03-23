reinitialize
load ./outputs/2L5E/2L5E_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 59
color red, structure and chain A and resi 68
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 73
color red, structure and chain A and resi 111
color red, structure and chain A and resi 112
color red, structure and chain A and resi 115
color red, structure and chain A and resi 116
color red, structure and chain A and resi 120
color red, structure and chain A and resi 121
color red, structure and chain A and resi 122
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [57, 58, 59, 68, 69, 70, 73, 111, 112, 115, 116, 120, 121, 122]
