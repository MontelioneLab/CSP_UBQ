reinitialize
load ./outputs/2mps/2mps_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 51
color red, structure and chain A and resi 54
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 61
color red, structure and chain A and resi 62
color red, structure and chain A and resi 72
color red, structure and chain A and resi 75
color red, structure and chain A and resi 91
color red, structure and chain A and resi 93
color red, structure and chain A and resi 94
color red, structure and chain A and resi 96
color red, structure and chain A and resi 99
color red, structure and chain A and resi 100
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [51, 54, 57, 58, 61, 62, 72, 75, 91, 93, 94, 96, 99, 100]
