reinitialize
load ./outputs/2KFH_15279/2kfh_csp.pdb, 2kfh_structure
hide everything, 2kfh_structure
show cartoon, 2kfh_structure
color cyan, 2kfh_structure and chain B
color gray30, 2kfh_structure and chain A
color red, 2kfh_structure and chain A and resi 41
color red, 2kfh_structure and chain A and resi 69
color red, 2kfh_structure and chain A and resi 74
color red, 2kfh_structure and chain A and resi 76
color red, 2kfh_structure and chain A and resi 77
color red, 2kfh_structure and chain A and resi 90
color red, 2kfh_structure and chain A and resi 91
color red, 2kfh_structure and chain A and resi 94
color red, 2kfh_structure and chain A and resi 137
set cartoon_transparency, 0.2, 2kfh_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 9
