reinitialize
load ./outputs/2LP0_25142/2LP0_csp.pdb, 2LP0_structure
hide everything, 2LP0_structure
show cartoon, 2LP0_structure
color cyan, 2LP0_structure and chain B
color gray30, 2LP0_structure and chain A
color red, 2LP0_structure and chain A and resi 18
color red, 2LP0_structure and chain A and resi 24
color red, 2LP0_structure and chain A and resi 69
color red, 2LP0_structure and chain A and resi 70
color red, 2LP0_structure and chain A and resi 71
color red, 2LP0_structure and chain A and resi 73
color red, 2LP0_structure and chain A and resi 74
color red, 2LP0_structure and chain A and resi 76
set cartoon_transparency, 0.2, 2LP0_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 8
