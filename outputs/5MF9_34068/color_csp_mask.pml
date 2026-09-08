reinitialize
load ./outputs/5MF9_34068/5MF9_csp.pdb, 5MF9_structure
hide everything, 5MF9_structure
show cartoon, 5MF9_structure
color cyan, 5MF9_structure and chain B
color gray30, 5MF9_structure and chain A
color red, 5MF9_structure and chain A and resi 34
color red, 5MF9_structure and chain A and resi 36
color red, 5MF9_structure and chain A and resi 40
color red, 5MF9_structure and chain A and resi 43
color red, 5MF9_structure and chain A and resi 47
color red, 5MF9_structure and chain A and resi 64
set cartoon_transparency, 0.2, 5MF9_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 6
