reinitialize
load ./outputs/2KID/2KID_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 92
color red, structure and chain A and resi 93
color red, structure and chain A and resi 97
color red, structure and chain A and resi 104
color red, structure and chain A and resi 105
color red, structure and chain A and resi 118
color red, structure and chain A and resi 120
color red, structure and chain A and resi 162
color red, structure and chain A and resi 163
color red, structure and chain A and resi 164
color red, structure and chain A and resi 165
color red, structure and chain A and resi 166
color red, structure and chain A and resi 168
color red, structure and chain A and resi 169
color red, structure and chain A and resi 182
color red, structure and chain A and resi 184
color red, structure and chain A and resi 194
color red, structure and chain A and resi 197
color cyan, structure and chain C
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: C
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [92, 93, 97, 104, 105, 118, 120, 162, 163, 164, 165, 166, 168, 169, 182, 184, 194, 197]
