reinitialize
load ./outputs/2JMF/2JMF_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 534
color red, structure and chain A and resi 536
color red, structure and chain A and resi 538
color red, structure and chain A and resi 540
color red, structure and chain A and resi 542
color red, structure and chain A and resi 545
color red, structure and chain A and resi 547
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 7
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [534, 536, 538, 540, 542, 545, 547]
