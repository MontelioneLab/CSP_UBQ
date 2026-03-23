reinitialize
load ./outputs/2mkr/2mkr_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 57
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 86
color red, structure and chain A and resi 88
color red, structure and chain A and resi 89
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [49, 50, 51, 52, 57, 59, 60, 61, 86, 88, 89]
