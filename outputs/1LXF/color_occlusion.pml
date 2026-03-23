reinitialize
load ./outputs/1LXF/1LXF_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color gray30, structure and chain C
color red, structure and chain C and resi 20
color red, structure and chain C and resi 22
color red, structure and chain C and resi 23
color red, structure and chain C and resi 26
color red, structure and chain C and resi 41
color red, structure and chain C and resi 44
color red, structure and chain C and resi 45
color red, structure and chain C and resi 47
color red, structure and chain C and resi 48
color red, structure and chain C and resi 60
color red, structure and chain C and resi 64
color red, structure and chain C and resi 69
color red, structure and chain C and resi 72
color red, structure and chain C and resi 77
color red, structure and chain C and resi 80
color red, structure and chain C and resi 81
color red, structure and chain C and resi 84
color red, structure and chain C and resi 85
color red, structure and chain C and resi 88
color cyan, structure and chain I
set cartoon_transparency, 0.2, structure
# Binding Site Analysis Summary:
# Receptor chain: C
# Ligand chain: I
# Non-binding receptor residues: gray30
# Binding site residues (red): 19
# Binding site definition: is_occluded OR passes_ca_filter OR is_interacting OR min_any_atom_distance_sub_2A
# Binding site residue numbers: [20, 22, 23, 26, 41, 44, 45, 47, 48, 60, 64, 69, 72, 77, 80, 81, 84, 85, 88]
