reinitialize
load ./outputs/2mbh/2mbh_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 8
color red, structure and chain A and resi 44
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 68
color red, structure and chain A and resi 70
color red, structure and chain A and resi 72
color red, structure and chain A and resi 75
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 8, 44, 46, 47, 68, 70, 72, 75]
