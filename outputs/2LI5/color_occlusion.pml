reinitialize
load ./outputs/2LI5/2LI5_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 4
color red, structure and chain A and resi 5
color red, structure and chain A and resi 7
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 10
color red, structure and chain A and resi 20
color red, structure and chain A and resi 24
color red, structure and chain A and resi 28
color red, structure and chain A and resi 45
color red, structure and chain A and resi 46
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 55
color red, structure and chain A and resi 63
color red, structure and chain A and resi 66
color red, structure and chain A and resi 67
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 20
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [4, 5, 7, 8, 9, 10, 20, 24, 28, 45, 46, 48, 49, 50, 51, 52, 55, 63, 66, 67]
