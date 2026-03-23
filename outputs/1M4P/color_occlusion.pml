reinitialize
load ./outputs/1M4P/1M4P_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 58
color red, structure and chain A and resi 61
color red, structure and chain A and resi 68
color red, structure and chain A and resi 69
color red, structure and chain A and resi 70
color red, structure and chain A and resi 71
color red, structure and chain A and resi 92
color red, structure and chain A and resi 95
color red, structure and chain A and resi 139
color red, structure and chain A and resi 141
color red, structure and chain A and resi 142
color red, structure and chain A and resi 143
color red, structure and chain A and resi 144
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 13
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [58, 61, 68, 69, 70, 71, 92, 95, 139, 141, 142, 143, 144]
