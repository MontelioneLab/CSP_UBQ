reinitialize
load ./outputs/7jq8/7jq8_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 5
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 12
color red, structure and chain A and resi 16
color red, structure and chain A and resi 20
color red, structure and chain A and resi 21
color red, structure and chain A and resi 23
color red, structure and chain A and resi 24
color red, structure and chain A and resi 27
color red, structure and chain A and resi 43
color red, structure and chain A and resi 44
color red, structure and chain A and resi 45
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 50
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [5, 8, 9, 12, 16, 20, 21, 23, 24, 27, 43, 44, 45, 46, 47, 48, 49, 50]
