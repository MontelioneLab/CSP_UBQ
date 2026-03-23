reinitialize
load ./outputs/2L1R/2L1R_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 23
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 41
color red, structure and chain A and resi 45
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 56
color red, structure and chain A and resi 60
color red, structure and chain A and resi 65
color red, structure and chain A and resi 67
color red, structure and chain A and resi 71
color red, structure and chain A and resi 80
color red, structure and chain A and resi 81
color red, structure and chain A and resi 84
color red, structure and chain A and resi 85
color red, structure and chain A and resi 91
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [23, 26, 27, 41, 45, 47, 48, 56, 60, 65, 67, 71, 80, 81, 84, 85, 91]
