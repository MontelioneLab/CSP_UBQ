reinitialize
load ./outputs/2N8T/2N8T_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 21
color red, structure and chain A and resi 23
color red, structure and chain A and resi 25
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [15, 16, 21, 23, 25, 27, 28, 29, 30, 31, 32]
