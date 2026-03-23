reinitialize
load ./outputs/6K5R/6K5R_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 17
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 31
color red, structure and chain A and resi 32
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 42
color red, structure and chain A and resi 47
color red, structure and chain A and resi 50
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [17, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 42, 47, 50]
