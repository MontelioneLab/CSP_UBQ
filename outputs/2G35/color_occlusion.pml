reinitialize
load ./outputs/2G35/2G35_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 50
color red, structure and chain A and resi 54
color red, structure and chain A and resi 55
color red, structure and chain A and resi 56
color red, structure and chain A and resi 63
color red, structure and chain A and resi 92
color red, structure and chain A and resi 95
color red, structure and chain A and resi 96
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [50, 54, 55, 56, 63, 92, 95, 96]
