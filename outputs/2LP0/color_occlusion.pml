reinitialize
load ./outputs/2LP0/2LP0_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 24
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 63
color red, structure and chain A and resi 67
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 74
color red, structure and chain A and resi 75
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [17, 18, 19, 20, 21, 22, 24, 59, 60, 63, 67, 70, 71, 74, 75]
