reinitialize
load ./outputs/2M5A_1/2M5A_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 5
color red, structure and chain A and resi 9
color red, structure and chain A and resi 10
color red, structure and chain A and resi 13
color red, structure and chain A and resi 14
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 20
color red, structure and chain A and resi 24
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [5, 9, 10, 13, 14, 17, 18, 20, 24, 27, 28, 31, 32, 34, 35, 36]
