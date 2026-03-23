reinitialize
load ./outputs/2L4T/2L4T_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 20
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 39
color red, structure and chain A and resi 40
color red, structure and chain A and resi 42
color red, structure and chain A and resi 43
color red, structure and chain A and resi 59
color red, structure and chain A and resi 90
color red, structure and chain A and resi 97
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [20, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42, 43, 59, 90, 97]
