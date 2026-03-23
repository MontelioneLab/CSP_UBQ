reinitialize
load ./outputs/6GC3/6GC3_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 30
color red, structure and chain A and resi 67
color red, structure and chain A and resi 70
color red, structure and chain A and resi 123
color red, structure and chain A and resi 127
color red, structure and chain A and resi 130
color red, structure and chain A and resi 133
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [20, 21, 22, 24, 25, 26, 27, 28, 30, 67, 70, 123, 127, 130, 133]
