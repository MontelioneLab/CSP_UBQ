reinitialize
load ./outputs/2LZ6_1/2LZ6_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 8
color red, structure and chain A and resi 42
color red, structure and chain A and resi 44
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 68
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 72
color red, structure and chain A and resi 74
color red, structure and chain A and resi 75
color red, structure and chain A and resi 76
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 8, 42, 44, 47, 48, 49, 68, 70, 71, 72, 74, 75, 76]
