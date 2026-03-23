reinitialize
load ./outputs/2N8T/2N8T_csp.pdb, 2N8T_structure
hide everything, 2N8T_structure
show cartoon, 2N8T_structure
color cyan, 2N8T_structure and chain B
color gray30, 2N8T_structure and chain A
color red, 2N8T_structure and chain A and resi 11
color red, 2N8T_structure and chain A and resi 15
color red, 2N8T_structure and chain A and resi 19
color red, 2N8T_structure and chain A and resi 20
color red, 2N8T_structure and chain A and resi 21
color red, 2N8T_structure and chain A and resi 22
color red, 2N8T_structure and chain A and resi 23
color red, 2N8T_structure and chain A and resi 28
color red, 2N8T_structure and chain A and resi 29
color red, 2N8T_structure and chain A and resi 30
color red, 2N8T_structure and chain A and resi 31
color red, 2N8T_structure and chain A and resi 32
color red, 2N8T_structure and chain A and resi 33
color red, 2N8T_structure and chain A and resi 34
set cartoon_transparency, 0.2, 2N8T_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 14
