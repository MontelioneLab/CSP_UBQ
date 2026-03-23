reinitialize
load ./outputs/2MMA/2MMA_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 46
color red, structure and chain B and resi 47
color red, structure and chain B and resi 48
color red, structure and chain B and resi 50
color red, structure and chain B and resi 52
color red, structure and chain B and resi 62
color red, structure and chain B and resi 66
color red, structure and chain B and resi 67
color red, structure and chain B and resi 68
color red, structure and chain B and resi 69
color red, structure and chain B and resi 70
color red, structure and chain B and resi 71
color red, structure and chain B and resi 72
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [46, 47, 48, 50, 52, 62, 66, 67, 68, 69, 70, 71, 72]
