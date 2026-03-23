reinitialize
load ./outputs/2LE8/2LE8_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 46
color red, structure and chain A and resi 47
color red, structure and chain A and resi 51
color red, structure and chain A and resi 54
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 60
color red, structure and chain A and resi 63
color red, structure and chain A and resi 64
color red, structure and chain A and resi 67
color red, structure and chain A and resi 71
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [46, 47, 51, 54, 57, 58, 60, 63, 64, 67, 71]
