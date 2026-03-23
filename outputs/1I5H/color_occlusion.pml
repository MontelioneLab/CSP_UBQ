reinitialize
load ./outputs/1I5H/1I5H_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain W
color red, structure and chain W and resi 470
color red, structure and chain W and resi 476
color red, structure and chain W and resi 478
color red, structure and chain W and resi 483
color red, structure and chain W and resi 485
color red, structure and chain W and resi 486
color red, structure and chain W and resi 487
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: W
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 7
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [470, 476, 478, 483, 485, 486, 487]
