reinitialize
load ./outputs/2mkp/2mkp_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain C
color red, structure and chain C and resi 15
color red, structure and chain C and resi 18
color red, structure and chain C and resi 19
color red, structure and chain C and resi 22
color red, structure and chain C and resi 44
color red, structure and chain C and resi 45
color red, structure and chain C and resi 48
color red, structure and chain C and resi 63
color red, structure and chain C and resi 65
color red, structure and chain C and resi 67
color red, structure and chain C and resi 69
color red, structure and chain C and resi 71
color red, structure and chain C and resi 76
color red, structure and chain C and resi 80
color red, structure and chain C and resi 81
color red, structure and chain C and resi 84
color red, structure and chain C and resi 86
color red, structure and chain C and resi 87
color cyan, structure and chain I
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: C
# Ligand chain: I
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [15, 18, 19, 22, 44, 45, 48, 63, 65, 67, 69, 71, 76, 80, 81, 84, 86, 87]
