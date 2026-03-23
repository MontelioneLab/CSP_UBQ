reinitialize
load ./outputs/2lsj/2lsj_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 39
color red, structure and chain A and resi 40
color red, structure and chain A and resi 42
color red, structure and chain A and resi 43
color red, structure and chain A and resi 47
color red, structure and chain A and resi 51
color red, structure and chain A and resi 53
color red, structure and chain A and resi 54
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [28, 29, 39, 40, 42, 43, 47, 51, 53, 54]
