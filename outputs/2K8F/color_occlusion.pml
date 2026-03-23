reinitialize
load ./outputs/2K8F/2K8F_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 8
color red, structure and chain A and resi 9
color red, structure and chain A and resi 12
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 41
color red, structure and chain A and resi 42
color red, structure and chain A and resi 44
color red, structure and chain A and resi 45
color red, structure and chain A and resi 47
color red, structure and chain A and resi 48
color red, structure and chain A and resi 52
color red, structure and chain A and resi 57
color red, structure and chain A and resi 59
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [8, 9, 12, 15, 16, 38, 39, 41, 42, 44, 45, 47, 48, 52, 57, 59]
