reinitialize
load ./outputs/2N3K/2N3K_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 615
color red, structure and chain A and resi 616
color red, structure and chain A and resi 627
color red, structure and chain A and resi 630
color red, structure and chain A and resi 631
color red, structure and chain A and resi 634
color red, structure and chain A and resi 650
color red, structure and chain A and resi 651
color red, structure and chain A and resi 652
color red, structure and chain A and resi 653
color red, structure and chain A and resi 654
color red, structure and chain A and resi 655
color red, structure and chain A and resi 656
color red, structure and chain A and resi 657
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [615, 616, 627, 630, 631, 634, 650, 651, 652, 653, 654, 655, 656, 657]
