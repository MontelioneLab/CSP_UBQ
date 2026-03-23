reinitialize
load ./outputs/2M0U/2M0U_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 22
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 40
color red, structure and chain A and resi 72
color red, structure and chain A and resi 76
color red, structure and chain A and resi 79
color red, structure and chain A and resi 80
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [22, 23, 24, 25, 26, 27, 28, 29, 40, 72, 76, 79, 80]
