reinitialize
load ./outputs/2RSE_2/2RSE_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 1085
color red, structure and chain B and resi 1095
color red, structure and chain B and resi 1099
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 3
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1085, 1095, 1099]
