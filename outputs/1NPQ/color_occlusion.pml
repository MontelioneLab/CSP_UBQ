reinitialize
load ./outputs/1NPQ/1NPQ_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 18
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 25
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 30
color red, structure and chain A and resi 41
color red, structure and chain A and resi 45
color red, structure and chain A and resi 46
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 66
color red, structure and chain A and resi 82
color red, structure and chain A and resi 85
color red, structure and chain A and resi 86
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 19
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [18, 20, 21, 22, 23, 24, 25, 28, 29, 30, 41, 45, 46, 48, 49, 66, 82, 85, 86]
