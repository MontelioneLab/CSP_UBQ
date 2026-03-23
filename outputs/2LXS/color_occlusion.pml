reinitialize
load ./outputs/2LXS/2LXS_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 607
color red, structure and chain A and resi 611
color red, structure and chain A and resi 612
color red, structure and chain A and resi 614
color red, structure and chain A and resi 624
color red, structure and chain A and resi 631
color red, structure and chain A and resi 638
color red, structure and chain A and resi 656
color red, structure and chain A and resi 660
color red, structure and chain A and resi 663
color red, structure and chain A and resi 667
color red, structure and chain A and resi 668
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 12
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [607, 611, 612, 614, 624, 631, 638, 656, 660, 663, 667, 668]
