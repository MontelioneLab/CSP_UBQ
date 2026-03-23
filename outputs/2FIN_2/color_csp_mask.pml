reinitialize
load ./outputs/2FIN_2/2FIN_csp.pdb, 2FIN_structure
hide everything, 2FIN_structure
show cartoon, 2FIN_structure
color cyan, 2FIN_structure and chain A
color gray30, 2FIN_structure and chain B
color red, 2FIN_structure and chain B and resi 9
color red, 2FIN_structure and chain B and resi 10
color red, 2FIN_structure and chain B and resi 13
color red, 2FIN_structure and chain B and resi 14
color red, 2FIN_structure and chain B and resi 15
color red, 2FIN_structure and chain B and resi 16
color red, 2FIN_structure and chain B and resi 17
color red, 2FIN_structure and chain B and resi 18
color red, 2FIN_structure and chain B and resi 20
color red, 2FIN_structure and chain B and resi 24
color red, 2FIN_structure and chain B and resi 27
color red, 2FIN_structure and chain B and resi 28
color red, 2FIN_structure and chain B and resi 44
color red, 2FIN_structure and chain B and resi 45
color red, 2FIN_structure and chain B and resi 46
color red, 2FIN_structure and chain B and resi 48
color red, 2FIN_structure and chain B and resi 49
color red, 2FIN_structure and chain B and resi 51
color red, 2FIN_structure and chain B and resi 55
color red, 2FIN_structure and chain B and resi 56
color red, 2FIN_structure and chain B and resi 65
set cartoon_transparency, 0.2, 2FIN_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 21
