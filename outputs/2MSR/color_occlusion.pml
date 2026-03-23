reinitialize
load ./outputs/2MSR/2MSR_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 359
color red, structure and chain B and resi 360
color red, structure and chain B and resi 363
color red, structure and chain B and resi 399
color red, structure and chain B and resi 402
color red, structure and chain B and resi 406
color red, structure and chain B and resi 407
color red, structure and chain B and resi 408
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 8
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [359, 360, 363, 399, 402, 406, 407, 408]
