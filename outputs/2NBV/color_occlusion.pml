reinitialize
load ./outputs/2NBV/2NBV_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 54
color red, structure and chain A and resi 55
color red, structure and chain A and resi 56
color red, structure and chain A and resi 73
color red, structure and chain A and resi 74
color red, structure and chain A and resi 75
color red, structure and chain A and resi 76
color red, structure and chain A and resi 77
color red, structure and chain A and resi 78
color red, structure and chain A and resi 79
color red, structure and chain A and resi 98
color red, structure and chain A and resi 100
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [54, 55, 56, 73, 74, 75, 76, 77, 78, 79, 98, 100]
