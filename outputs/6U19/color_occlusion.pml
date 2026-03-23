reinitialize
load ./outputs/6U19/6U19_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 25
color red, structure and chain B and resi 26
color red, structure and chain B and resi 27
color red, structure and chain B and resi 29
color red, structure and chain B and resi 30
color red, structure and chain B and resi 33
color red, structure and chain B and resi 44
color red, structure and chain B and resi 49
color red, structure and chain B and resi 54
color red, structure and chain B and resi 65
color red, structure and chain B and resi 66
color red, structure and chain B and resi 69
color red, structure and chain B and resi 70
color red, structure and chain B and resi 76
color red, structure and chain B and resi 77
color red, structure and chain B and resi 83
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 16
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [25, 26, 27, 29, 30, 33, 44, 49, 54, 65, 66, 69, 70, 76, 77, 83]
