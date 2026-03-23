reinitialize
load ./outputs/6B1G/6B1G_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 1
color red, structure and chain A and resi 4
color red, structure and chain A and resi 14
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 27
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 36
color red, structure and chain A and resi 80
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1, 4, 14, 15, 16, 17, 18, 19, 20, 21, 22, 27, 29, 30, 33, 34, 36, 80]
