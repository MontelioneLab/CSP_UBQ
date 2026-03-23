reinitialize
load ./outputs/2KVQ/2KVQ_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain G
color red, structure and chain G and resi 140
color red, structure and chain G and resi 141
color red, structure and chain G and resi 164
color red, structure and chain G and resi 165
color red, structure and chain G and resi 167
color red, structure and chain G and resi 169
color red, structure and chain G and resi 170
color red, structure and chain G and resi 171
color red, structure and chain G and resi 172
color cyan, structure and chain E
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: G
# Ligand chain: E
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [140, 141, 164, 165, 167, 169, 170, 171, 172]
