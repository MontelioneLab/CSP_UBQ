reinitialize
load ./outputs/2oj2/2oj2_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 32
color red, structure and chain A and resi 37
color red, structure and chain A and resi 40
color red, structure and chain A and resi 56
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 72
color red, structure and chain A and resi 74
color red, structure and chain A and resi 76
color red, structure and chain A and resi 77
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [32, 37, 40, 56, 57, 58, 59, 60, 72, 74, 76, 77]
