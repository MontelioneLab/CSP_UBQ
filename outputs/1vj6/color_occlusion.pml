reinitialize
load ./outputs/1VJ6/1VJ6_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 20
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 45
color red, structure and chain A and resi 78
color red, structure and chain A and resi 82
color red, structure and chain A and resi 86
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [20, 24, 25, 26, 27, 28, 29, 30, 45, 78, 82, 86]
