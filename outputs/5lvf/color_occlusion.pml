reinitialize
load ./outputs/5lvf/5lvf_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 65
color red, structure and chain A and resi 66
color red, structure and chain A and resi 69
color red, structure and chain A and resi 72
color red, structure and chain A and resi 73
color red, structure and chain A and resi 75
color red, structure and chain A and resi 108
color red, structure and chain A and resi 109
color red, structure and chain A and resi 112
color red, structure and chain A and resi 116
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [17, 18, 19, 20, 65, 66, 69, 72, 73, 75, 108, 109, 112, 116]
