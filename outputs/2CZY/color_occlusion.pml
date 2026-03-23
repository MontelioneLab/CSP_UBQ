reinitialize
load ./outputs/2CZY/2CZY_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 32
color red, structure and chain A and resi 34
color red, structure and chain A and resi 38
color red, structure and chain A and resi 41
color red, structure and chain A and resi 59
color red, structure and chain A and resi 62
color red, structure and chain A and resi 63
color red, structure and chain A and resi 65
color red, structure and chain A and resi 66
color red, structure and chain A and resi 75
color red, structure and chain A and resi 93
color red, structure and chain A and resi 96
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [32, 34, 38, 41, 59, 62, 63, 65, 66, 75, 93, 96]
