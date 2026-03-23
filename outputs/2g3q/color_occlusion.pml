reinitialize
load ./outputs/2g3q/2g3q_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain B
color red, structure and chain B and resi 6
color red, structure and chain B and resi 42
color red, structure and chain B and resi 44
color red, structure and chain B and resi 46
color red, structure and chain B and resi 47
color red, structure and chain B and resi 49
color red, structure and chain B and resi 68
color red, structure and chain B and resi 70
color red, structure and chain B and resi 71
color red, structure and chain B and resi 73
color cyan, structure and chain A
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Non-binding receptor residues: gray30
# Binding site residues (red): 10
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [6, 42, 44, 46, 47, 49, 68, 70, 71, 73]
