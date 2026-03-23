reinitialize
load ./outputs/2n9e/2n9e_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 30
color red, structure and chain B and resi 31
color red, structure and chain B and resi 32
color red, structure and chain B and resi 33
color red, structure and chain B and resi 34
color red, structure and chain B and resi 35
color red, structure and chain B and resi 37
color red, structure and chain B and resi 38
color red, structure and chain B and resi 42
color red, structure and chain B and resi 50
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [30, 31, 32, 33, 34, 35, 37, 38, 42, 50]
