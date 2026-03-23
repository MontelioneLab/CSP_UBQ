reinitialize
load ./outputs/6h8c/6h8c_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 4
color red, structure and chain A and resi 5
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 13
color red, structure and chain A and resi 17
color red, structure and chain A and resi 25
color red, structure and chain A and resi 28
color red, structure and chain A and resi 32
color red, structure and chain A and resi 34
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 63
color red, structure and chain A and resi 66
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 19
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [4, 5, 8, 9, 13, 17, 25, 28, 32, 34, 46, 47, 48, 49, 50, 51, 52, 63, 66]
