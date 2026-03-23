reinitialize
load ./outputs/2BN5/2BN5_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 660
color red, structure and chain A and resi 661
color red, structure and chain A and resi 665
color red, structure and chain A and resi 674
color red, structure and chain A and resi 677
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 5
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [660, 661, 665, 674, 677]
