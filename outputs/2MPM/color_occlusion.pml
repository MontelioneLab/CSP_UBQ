reinitialize
load ./outputs/2MPM/2MPM_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 13
color red, structure and chain A and resi 14
color red, structure and chain A and resi 15
color red, structure and chain A and resi 16
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 45
color red, structure and chain A and resi 47
color red, structure and chain A and resi 49
color red, structure and chain A and resi 55
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [13, 14, 15, 16, 17, 18, 45, 47, 49, 55]
