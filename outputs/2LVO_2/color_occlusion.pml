reinitialize
load ./outputs/2LVO_2/2LVO_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 6
color red, structure and chain A and resi 8
color red, structure and chain A and resi 42
color red, structure and chain A and resi 44
color red, structure and chain A and resi 47
color red, structure and chain A and resi 68
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 72
color red, structure and chain A and resi 73
color red, structure and chain A and resi 74
color cyan, structure and chain C
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: C
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 8, 42, 44, 47, 68, 70, 71, 72, 73, 74]
