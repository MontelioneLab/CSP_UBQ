reinitialize
load ./outputs/2YS5_11095/2YS5_csp.pdb, 2YS5_structure
hide everything, 2YS5_structure
show cartoon, 2YS5_structure
color cyan, 2YS5_structure and chain B
color gray30, 2YS5_structure and chain A
color red, 2YS5_structure and chain A and resi 10
color red, 2YS5_structure and chain A and resi 25
color red, 2YS5_structure and chain A and resi 28
color red, 2YS5_structure and chain A and resi 33
color red, 2YS5_structure and chain A and resi 34
color red, 2YS5_structure and chain A and resi 35
color red, 2YS5_structure and chain A and resi 38
color red, 2YS5_structure and chain A and resi 49
color red, 2YS5_structure and chain A and resi 54
color red, 2YS5_structure and chain A and resi 55
color red, 2YS5_structure and chain A and resi 60
color red, 2YS5_structure and chain A and resi 110
set cartoon_transparency, 0.2, 2YS5_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 12
