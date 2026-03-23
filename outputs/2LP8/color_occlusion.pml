reinitialize
load ./outputs/2LP8/2LP8_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 57
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 71
color red, structure and chain A and resi 72
color red, structure and chain A and resi 85
color red, structure and chain A and resi 86
color red, structure and chain A and resi 89
color red, structure and chain A and resi 90
color red, structure and chain A and resi 98
color red, structure and chain A and resi 99
color red, structure and chain A and resi 106
color red, structure and chain A and resi 155
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [57, 60, 61, 71, 72, 85, 86, 89, 90, 98, 99, 106, 155]
