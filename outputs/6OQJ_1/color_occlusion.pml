reinitialize
load ./outputs/6OQJ_1/6OQJ_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 37
color red, structure and chain A and resi 39
color red, structure and chain A and resi 41
color red, structure and chain A and resi 42
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 62
color red, structure and chain A and resi 63
color red, structure and chain A and resi 67
color red, structure and chain A and resi 76
color red, structure and chain A and resi 77
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [37, 39, 41, 42, 46, 47, 60, 61, 62, 63, 67, 76, 77]
