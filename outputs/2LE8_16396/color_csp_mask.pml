reinitialize
load ./outputs/2LE8_16396/2LE8_csp.pdb, 2LE8_structure
hide everything, 2LE8_structure
show cartoon, 2LE8_structure
color cyan, 2LE8_structure and chain B
color gray30, 2LE8_structure and chain A
color red, 2LE8_structure and chain A and resi 47
color red, 2LE8_structure and chain A and resi 51
color red, 2LE8_structure and chain A and resi 52
color red, 2LE8_structure and chain A and resi 53
color red, 2LE8_structure and chain A and resi 54
color red, 2LE8_structure and chain A and resi 55
color red, 2LE8_structure and chain A and resi 61
color red, 2LE8_structure and chain A and resi 68
set cartoon_transparency, 0.2, 2LE8_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 8
