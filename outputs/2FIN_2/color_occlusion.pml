reinitialize
load ./outputs/2FIN_2/2FIN_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 8
color red, structure and chain B and resi 10
color red, structure and chain B and resi 12
color red, structure and chain B and resi 13
color red, structure and chain B and resi 14
color red, structure and chain B and resi 15
color red, structure and chain B and resi 16
color red, structure and chain B and resi 17
color red, structure and chain B and resi 18
color red, structure and chain B and resi 24
color red, structure and chain B and resi 45
color red, structure and chain B and resi 46
color red, structure and chain B and resi 47
color red, structure and chain B and resi 48
color red, structure and chain B and resi 49
color red, structure and chain B and resi 50
color red, structure and chain B and resi 51
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 17
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [8, 10, 12, 13, 14, 15, 16, 17, 18, 24, 45, 46, 47, 48, 49, 50, 51]
