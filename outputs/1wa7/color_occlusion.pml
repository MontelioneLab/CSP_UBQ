reinitialize
load ./outputs/1wa7/1wa7_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain A
color red, structure and chain A and resi 17
color red, structure and chain A and resi 18
color red, structure and chain A and resi 19
color red, structure and chain A and resi 23
color red, structure and chain A and resi 26
color red, structure and chain A and resi 31
color red, structure and chain A and resi 41
color red, structure and chain A and resi 42
color red, structure and chain A and resi 43
color red, structure and chain A and resi 44
color red, structure and chain A and resi 57
color red, structure and chain A and resi 59
color red, structure and chain A and resi 60
color red, structure and chain A and resi 61
color red, structure and chain A and resi 62
color cyan, structure and chain B
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Non-binding receptor residues: gray30
# Binding site residues (red): 15
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [17, 18, 19, 23, 26, 31, 41, 42, 43, 44, 57, 59, 60, 61, 62]
