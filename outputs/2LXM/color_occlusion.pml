reinitialize
load ./outputs/2LXM/2LXM_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 15
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 22
color red, structure and chain A and resi 26
color red, structure and chain A and resi 27
color red, structure and chain A and resi 30
color red, structure and chain A and resi 96
color red, structure and chain A and resi 115
color red, structure and chain A and resi 116
color red, structure and chain A and resi 119
color red, structure and chain A and resi 144
color red, structure and chain A and resi 146
color red, structure and chain A and resi 147
color red, structure and chain A and resi 150
color red, structure and chain A and resi 153
color red, structure and chain A and resi 154
color red, structure and chain A and resi 157
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [15, 18, 19, 22, 26, 27, 30, 96, 115, 116, 119, 144, 146, 147, 150, 153, 154, 157]
