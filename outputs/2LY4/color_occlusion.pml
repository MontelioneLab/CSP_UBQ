reinitialize
load ./outputs/2LY4/2LY4_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 12
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 23
color red, structure and chain A and resi 27
color red, structure and chain A and resi 34
color red, structure and chain A and resi 37
color red, structure and chain A and resi 41
color red, structure and chain A and resi 45
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [12, 15, 16, 19, 20, 23, 27, 34, 37, 41, 45]
