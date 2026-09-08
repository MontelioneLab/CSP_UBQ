reinitialize
load ./outputs/2ND0_15125/2nd0_csp.pdb, 2nd0_structure
hide everything, 2nd0_structure
show cartoon, 2nd0_structure
color cyan, 2nd0_structure and chain B
color gray30, 2nd0_structure and chain A
color red, 2nd0_structure and chain A and resi 3
color red, 2nd0_structure and chain A and resi 4
color red, 2nd0_structure and chain A and resi 6
color red, 2nd0_structure and chain A and resi 52
color red, 2nd0_structure and chain A and resi 82
set cartoon_transparency, 0.2, 2nd0_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 5
