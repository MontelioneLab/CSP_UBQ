reinitialize
load ./outputs/2M0G_1/2M0G_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 30
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 40
color red, structure and chain A and resi 42
color red, structure and chain A and resi 49
color red, structure and chain A and resi 53
color red, structure and chain A and resi 57
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [15, 16, 21, 22, 23, 24, 27, 28, 30, 33, 34, 38, 39, 40, 42, 49, 53, 57]
