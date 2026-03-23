reinitialize
load ./outputs/2m0v/2m0v_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 163
color red, structure and chain A and resi 164
color red, structure and chain A and resi 165
color red, structure and chain A and resi 166
color red, structure and chain A and resi 167
color red, structure and chain A and resi 168
color red, structure and chain A and resi 169
color red, structure and chain A and resi 170
color red, structure and chain A and resi 180
color red, structure and chain A and resi 212
color red, structure and chain A and resi 219
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [163, 164, 165, 166, 167, 168, 169, 170, 180, 212, 219]
