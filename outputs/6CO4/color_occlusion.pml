reinitialize
load ./outputs/6CO4/6CO4_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 27
color red, structure and chain A and resi 31
color red, structure and chain A and resi 33
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 40
color red, structure and chain A and resi 42
color red, structure and chain A and resi 44
color red, structure and chain A and resi 64
color red, structure and chain A and resi 87
color red, structure and chain A and resi 88
color red, structure and chain A and resi 92
color red, structure and chain A and resi 106
color red, structure and chain A and resi 110
color red, structure and chain A and resi 111
color red, structure and chain A and resi 112
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [27, 31, 33, 36, 37, 38, 39, 40, 42, 44, 64, 87, 88, 92, 106, 110, 111, 112]
