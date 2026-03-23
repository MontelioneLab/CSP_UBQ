reinitialize
load ./outputs/2N9P/2N9P_csp.pdb, 2N9P_structure
hide everything, 2N9P_structure
show cartoon, 2N9P_structure
color cyan, 2N9P_structure and chain C
color gray30, 2N9P_structure and chain A
color red, 2N9P_structure and chain A and resi 15
color red, 2N9P_structure and chain A and resi 23
color red, 2N9P_structure and chain A and resi 27
color red, 2N9P_structure and chain A and resi 28
color red, 2N9P_structure and chain A and resi 29
color red, 2N9P_structure and chain A and resi 32
color red, 2N9P_structure and chain A and resi 34
color red, 2N9P_structure and chain A and resi 35
color red, 2N9P_structure and chain A and resi 36
color red, 2N9P_structure and chain A and resi 37
color red, 2N9P_structure and chain A and resi 38
color red, 2N9P_structure and chain A and resi 39
set cartoon_transparency, 0.2, 2N9P_structure
# Color scheme:
# Ligand chain (C): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 12
