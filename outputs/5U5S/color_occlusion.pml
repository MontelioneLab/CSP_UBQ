reinitialize
load ./outputs/5U5S/5U5S_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 27
color red, structure and chain A and resi 28
color red, structure and chain A and resi 33
color red, structure and chain A and resi 38
color red, structure and chain A and resi 39
color red, structure and chain A and resi 40
color red, structure and chain A and resi 41
color red, structure and chain A and resi 42
color red, structure and chain A and resi 46
color red, structure and chain A and resi 82
color red, structure and chain A and resi 85
color red, structure and chain A and resi 86
color red, structure and chain A and resi 90
color red, structure and chain A and resi 91
color red, structure and chain A and resi 92
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [27, 28, 33, 38, 39, 40, 41, 42, 46, 82, 85, 86, 90, 91, 92]
