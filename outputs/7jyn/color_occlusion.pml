reinitialize
load ./outputs/7jyn/7jyn_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 30
color red, structure and chain A and resi 37
color red, structure and chain A and resi 41
color red, structure and chain A and resi 45
color red, structure and chain A and resi 48
color red, structure and chain A and resi 49
color red, structure and chain A and resi 52
color red, structure and chain A and resi 68
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 72
color red, structure and chain A and resi 73
color red, structure and chain A and resi 74
color red, structure and chain A and resi 75
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [30, 37, 41, 45, 48, 49, 52, 68, 69, 70, 71, 72, 73, 74, 75]
