reinitialize
load ./outputs/2K6Q/2K6Q_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 23
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 32
color red, structure and chain A and resi 33
color red, structure and chain A and resi 49
color red, structure and chain A and resi 51
color red, structure and chain A and resi 52
color red, structure and chain A and resi 53
color red, structure and chain A and resi 54
color red, structure and chain A and resi 55
color red, structure and chain A and resi 62
color red, structure and chain A and resi 66
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 108
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [19, 20, 23, 26, 27, 32, 33, 49, 51, 52, 53, 54, 55, 62, 66, 69, 70, 108]
