reinitialize
load ./outputs/2ROL/2ROL_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 494
color red, structure and chain A and resi 495
color red, structure and chain A and resi 512
color red, structure and chain A and resi 513
color red, structure and chain A and resi 514
color red, structure and chain A and resi 528
color red, structure and chain A and resi 529
color red, structure and chain A and resi 530
color red, structure and chain A and resi 531
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 9
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [494, 495, 512, 513, 514, 528, 529, 530, 531]
