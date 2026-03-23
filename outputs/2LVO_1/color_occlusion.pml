reinitialize
load ./outputs/2LVO_1/2LVO_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain C
color red, structure and chain C and resi 466
color red, structure and chain C and resi 467
color red, structure and chain C and resi 468
color red, structure and chain C and resi 470
color red, structure and chain C and resi 486
color red, structure and chain C and resi 487
color red, structure and chain C and resi 488
color red, structure and chain C and resi 490
color red, structure and chain C and resi 491
color red, structure and chain C and resi 494
color red, structure and chain C and resi 495
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: C
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [466, 467, 468, 470, 486, 487, 488, 490, 491, 494, 495]
