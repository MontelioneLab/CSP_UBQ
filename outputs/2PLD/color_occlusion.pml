reinitialize
load ./outputs/2PLD/2PLD_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 18
color red, structure and chain A and resi 37
color red, structure and chain A and resi 39
color red, structure and chain A and resi 46
color red, structure and chain A and resi 49
color red, structure and chain A and resi 56
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 59
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 88
color red, structure and chain A and resi 89
color red, structure and chain A and resi 90
color red, structure and chain A and resi 91
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [18, 37, 39, 46, 49, 56, 57, 58, 59, 69, 70, 71, 88, 89, 90, 91]
