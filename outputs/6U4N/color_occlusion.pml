reinitialize
load ./outputs/6U4N/6U4N_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 42
color red, structure and chain A and resi 43
color red, structure and chain A and resi 45
color red, structure and chain A and resi 46
color red, structure and chain A and resi 50
color red, structure and chain A and resi 55
color red, structure and chain A and resi 57
color red, structure and chain A and resi 58
color red, structure and chain A and resi 59
color red, structure and chain A and resi 73
color red, structure and chain A and resi 75
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [42, 43, 45, 46, 50, 55, 57, 58, 59, 73, 75]
