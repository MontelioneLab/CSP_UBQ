reinitialize
load ./outputs/2MH0/2MH0_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 1731
color red, structure and chain B and resi 1734
color red, structure and chain B and resi 1735
color red, structure and chain B and resi 1738
color red, structure and chain B and resi 1760
color red, structure and chain B and resi 1761
color red, structure and chain B and resi 1763
color red, structure and chain B and resi 1764
color red, structure and chain B and resi 1767
color red, structure and chain B and resi 1781
color red, structure and chain B and resi 1783
color red, structure and chain B and resi 1784
color red, structure and chain B and resi 1787
color red, structure and chain B and resi 1788
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 14
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [1731, 1734, 1735, 1738, 1760, 1761, 1763, 1764, 1767, 1781, 1783, 1784, 1787, 1788]
