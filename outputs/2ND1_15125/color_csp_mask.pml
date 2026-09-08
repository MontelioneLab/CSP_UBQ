reinitialize
load ./outputs/2ND1_15125/2nd1_csp.pdb, 2nd1_structure
hide everything, 2nd1_structure
show cartoon, 2nd1_structure
color cyan, 2nd1_structure and chain B
color gray30, 2nd1_structure and chain A
color red, 2nd1_structure and chain A and resi 3
color red, 2nd1_structure and chain A and resi 4
color red, 2nd1_structure and chain A and resi 6
color red, 2nd1_structure and chain A and resi 52
color red, 2nd1_structure and chain A and resi 82
set cartoon_transparency, 0.2, 2nd1_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 5
