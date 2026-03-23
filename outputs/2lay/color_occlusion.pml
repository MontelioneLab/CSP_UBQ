reinitialize
load ./outputs/2lay/2lay_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 178
color red, structure and chain A and resi 186
color red, structure and chain A and resi 188
color red, structure and chain A and resi 190
color red, structure and chain A and resi 192
color red, structure and chain A and resi 197
color red, structure and chain A and resi 199
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 7
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [178, 186, 188, 190, 192, 197, 199]
