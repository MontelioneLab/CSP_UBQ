reinitialize
load ./outputs/2kq0/2kq0_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 20
color red, structure and chain A and resi 24
color red, structure and chain A and resi 30
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 39
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [20, 24, 30, 31, 32, 35, 36, 37, 39]
