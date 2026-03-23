reinitialize
load ./outputs/2l0i/2l0i_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 62
color red, structure and chain A and resi 65
color red, structure and chain A and resi 66
color red, structure and chain A and resi 69
color red, structure and chain A and resi 108
color red, structure and chain A and resi 109
color red, structure and chain A and resi 112
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [16, 17, 18, 19, 62, 65, 66, 69, 108, 109, 112]
