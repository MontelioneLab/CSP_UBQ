reinitialize
load ./outputs/2l7u/2l7u_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 2
color red, structure and chain A and resi 3
color red, structure and chain A and resi 34
color red, structure and chain A and resi 41
color red, structure and chain A and resi 76
color red, structure and chain A and resi 78
color red, structure and chain A and resi 90
color red, structure and chain A and resi 92
color red, structure and chain A and resi 94
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [2, 3, 34, 41, 76, 78, 90, 92, 94]
