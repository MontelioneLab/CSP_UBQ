reinitialize
load ./outputs/2K6D/2K6D_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 7
color red, structure and chain B and resi 8
color red, structure and chain B and resi 9
color red, structure and chain B and resi 42
color red, structure and chain B and resi 44
color red, structure and chain B and resi 45
color red, structure and chain B and resi 46
color red, structure and chain B and resi 47
color red, structure and chain B and resi 48
color red, structure and chain B and resi 49
color red, structure and chain B and resi 68
color red, structure and chain B and resi 69
color red, structure and chain B and resi 70
color red, structure and chain B and resi 72
color red, structure and chain B and resi 73
color red, structure and chain B and resi 74
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [7, 8, 9, 42, 44, 45, 46, 47, 48, 49, 68, 69, 70, 72, 73, 74]
