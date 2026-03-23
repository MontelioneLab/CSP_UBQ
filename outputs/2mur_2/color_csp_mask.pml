reinitialize
load ./outputs/2mur_2/2mur_csp.pdb, 2mur_structure
hide everything, 2mur_structure
show cartoon, 2mur_structure
color cyan, 2mur_structure and chain B
color gray30, 2mur_structure and chain A
color red, 2mur_structure and chain A and resi 10
color red, 2mur_structure and chain A and resi 13
color red, 2mur_structure and chain A and resi 28
color red, 2mur_structure and chain A and resi 29
color red, 2mur_structure and chain A and resi 32
color red, 2mur_structure and chain A and resi 34
color red, 2mur_structure and chain A and resi 35
color red, 2mur_structure and chain A and resi 36
color red, 2mur_structure and chain A and resi 37
color red, 2mur_structure and chain A and resi 38
color red, 2mur_structure and chain A and resi 39
color red, 2mur_structure and chain A and resi 40
color red, 2mur_structure and chain A and resi 41
color red, 2mur_structure and chain A and resi 42
color red, 2mur_structure and chain A and resi 43
color red, 2mur_structure and chain A and resi 44
set cartoon_transparency, 0.2, 2mur_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 16
