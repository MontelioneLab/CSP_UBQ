reinitialize
load ./outputs/6IVU/6IVU_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 11
color red, structure and chain A and resi 16
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 49
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 8, 9, 11, 16, 20, 21, 22, 23, 24, 25, 26, 35, 36, 37, 38, 49]
