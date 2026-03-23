reinitialize
load ./outputs/2n1a/2n1a_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 17
color red, structure and chain A and resi 20
color red, structure and chain A and resi 40
color red, structure and chain A and resi 77
color red, structure and chain A and resi 81
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 5
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [17, 20, 40, 77, 81]
