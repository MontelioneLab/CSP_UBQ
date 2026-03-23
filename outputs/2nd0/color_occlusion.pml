reinitialize
load ./outputs/2nd0/2nd0_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 12
color red, structure and chain A and resi 19
color red, structure and chain A and resi 30
color red, structure and chain A and resi 31
color red, structure and chain A and resi 34
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 53
color red, structure and chain A and resi 54
color red, structure and chain A and resi 55
color red, structure and chain A and resi 56
color red, structure and chain A and resi 57
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [12, 19, 30, 31, 34, 49, 50, 51, 52, 53, 54, 55, 56, 57]
