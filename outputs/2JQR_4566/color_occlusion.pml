reinitialize
load ./outputs/2JQR_4566/2JQR_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 44
color red, structure and chain B and resi 46
color red, structure and chain B and resi 48
color red, structure and chain B and resi 50
color red, structure and chain B and resi 51
color red, structure and chain B and resi 52
color red, structure and chain B and resi 55
color red, structure and chain B and resi 76
color red, structure and chain B and resi 77
color red, structure and chain B and resi 79
color red, structure and chain B and resi 80
color red, structure and chain B and resi 92
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [44, 46, 48, 50, 51, 52, 55, 76, 77, 79, 80, 92]
