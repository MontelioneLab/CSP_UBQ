reinitialize
load ./outputs/2JMX/2JMX_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 14
color red, structure and chain A and resi 17
color red, structure and chain A and resi 21
color red, structure and chain A and resi 22
color red, structure and chain A and resi 25
color red, structure and chain A and resi 28
color red, structure and chain A and resi 29
color red, structure and chain A and resi 88
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [14, 17, 21, 22, 25, 28, 29, 88]
