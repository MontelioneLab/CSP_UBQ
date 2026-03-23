reinitialize
load ./outputs/2N9P/2N9P_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 11
color red, structure and chain A and resi 13
color red, structure and chain A and resi 14
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 27
color red, structure and chain A and resi 29
color red, structure and chain A and resi 32
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color cyan, structure and chain C
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: C
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [11, 13, 14, 15, 16, 27, 29, 32, 35, 36, 37, 38, 39]
