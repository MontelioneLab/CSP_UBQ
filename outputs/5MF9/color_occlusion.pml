reinitialize
load ./outputs/5MF9/5MF9_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 18
color red, structure and chain A and resi 21
color red, structure and chain A and resi 23
color red, structure and chain A and resi 25
color red, structure and chain A and resi 32
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 39
color red, structure and chain A and resi 41
color red, structure and chain A and resi 43
color red, structure and chain A and resi 46
color red, structure and chain A and resi 48
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [18, 21, 23, 25, 32, 33, 34, 36, 37, 39, 41, 43, 46, 48]
