reinitialize
load ./outputs/2MZP/2MZP_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain C
color red, structure and chain C and resi 26
color red, structure and chain C and resi 48
color red, structure and chain C and resi 65
color red, structure and chain C and resi 67
color red, structure and chain C and resi 69
color red, structure and chain C and resi 71
color red, structure and chain C and resi 76
color red, structure and chain C and resi 85
color cyan, structure and chain I
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: C
# Ligand chain: I
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [26, 48, 65, 67, 69, 71, 76, 85]
