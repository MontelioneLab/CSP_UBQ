reinitialize
load ./outputs/2IPA/2IPA_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 67
color red, structure and chain B and resi 69
color red, structure and chain B and resi 72
color red, structure and chain B and resi 88
color red, structure and chain B and resi 89
color red, structure and chain B and resi 90
color red, structure and chain B and resi 91
color red, structure and chain B and resi 92
color red, structure and chain B and resi 93
color red, structure and chain B and resi 94
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [67, 69, 72, 88, 89, 90, 91, 92, 93, 94]
