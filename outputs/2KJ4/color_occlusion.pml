reinitialize
load ./outputs/2KJ4/2KJ4_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 40
color red, structure and chain A and resi 42
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 50
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 62
color red, structure and chain A and resi 63
color red, structure and chain A and resi 67
color red, structure and chain A and resi 75
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
# Binding site residue numbers: [40, 42, 46, 47, 50, 60, 61, 62, 63, 67, 75, 76, 77]
