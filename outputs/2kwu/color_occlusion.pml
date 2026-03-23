reinitialize
load ./outputs/2kwu/2kwu_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 7
color red, structure and chain B and resi 8
color red, structure and chain B and resi 9
color red, structure and chain B and resi 10
color red, structure and chain B and resi 42
color red, structure and chain B and resi 44
color red, structure and chain B and resi 49
color red, structure and chain B and resi 68
color red, structure and chain B and resi 72
color red, structure and chain B and resi 73
color red, structure and chain B and resi 75
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [7, 8, 9, 10, 42, 44, 49, 68, 72, 73, 75]
