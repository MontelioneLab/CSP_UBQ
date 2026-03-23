reinitialize
load ./outputs/2KYL/2KYL_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 41
color red, structure and chain A and resi 44
color red, structure and chain A and resi 73
color red, structure and chain A and resi 80
color red, structure and chain A and resi 81
color red, structure and chain A and resi 94
color red, structure and chain A and resi 96
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 41, 44, 73, 80, 81, 94, 96]
