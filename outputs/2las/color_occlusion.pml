reinitialize
load ./outputs/2las/2las_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 21
color red, structure and chain A and resi 33
color red, structure and chain A and resi 34
color red, structure and chain A and resi 35
color red, structure and chain A and resi 36
color red, structure and chain A and resi 37
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 41
color red, structure and chain A and resi 42
color red, structure and chain A and resi 46
color red, structure and chain A and resi 50
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [21, 33, 34, 35, 36, 37, 38, 39, 41, 42, 46, 50]
