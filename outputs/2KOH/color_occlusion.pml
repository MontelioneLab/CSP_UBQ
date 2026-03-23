reinitialize
load ./outputs/2KOH/2KOH_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 600
color red, structure and chain A and resi 601
color red, structure and chain A and resi 602
color red, structure and chain A and resi 603
color red, structure and chain A and resi 604
color red, structure and chain A and resi 605
color red, structure and chain A and resi 606
color red, structure and chain A and resi 607
color red, structure and chain A and resi 608
color red, structure and chain A and resi 609
color red, structure and chain A and resi 611
color red, structure and chain A and resi 612
color red, structure and chain A and resi 622
color red, structure and chain A and resi 625
color red, structure and chain A and resi 659
color red, structure and chain A and resi 662
color red, structure and chain A and resi 663
color red, structure and chain A and resi 666
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 18
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 611, 612, 622, 625, 659, 662, 663, 666]
