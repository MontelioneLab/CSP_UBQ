reinitialize
load ./outputs/2L0Y/2L0Y_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 9
color red, structure and chain A and resi 10
color red, structure and chain A and resi 11
color red, structure and chain A and resi 23
color red, structure and chain A and resi 27
color red, structure and chain A and resi 38
color red, structure and chain A and resi 43
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 49
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [9, 10, 11, 23, 27, 38, 43, 46, 47, 49]
