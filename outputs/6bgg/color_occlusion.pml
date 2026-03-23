reinitialize
load ./outputs/6bgg/6bgg_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 577
color red, structure and chain B and resi 578
color red, structure and chain B and resi 581
color red, structure and chain B and resi 592
color red, structure and chain B and resi 613
color red, structure and chain B and resi 614
color red, structure and chain B and resi 615
color red, structure and chain B and resi 616
color red, structure and chain B and resi 617
color red, structure and chain B and resi 618
color red, structure and chain B and resi 619
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 11
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [577, 578, 581, 592, 613, 614, 615, 616, 617, 618, 619]
