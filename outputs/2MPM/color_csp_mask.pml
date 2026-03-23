reinitialize
load ./outputs/2MPM/2MPM_csp.pdb, 2MPM_structure
hide everything, 2MPM_structure
show cartoon, 2MPM_structure
color cyan, 2MPM_structure and chain B
color gray30, 2MPM_structure and chain A
color red, 2MPM_structure and chain A and resi 14
color red, 2MPM_structure and chain A and resi 16
color red, 2MPM_structure and chain A and resi 17
color red, 2MPM_structure and chain A and resi 18
color red, 2MPM_structure and chain A and resi 26
color red, 2MPM_structure and chain A and resi 27
color red, 2MPM_structure and chain A and resi 28
color red, 2MPM_structure and chain A and resi 29
color red, 2MPM_structure and chain A and resi 34
color red, 2MPM_structure and chain A and resi 38
color red, 2MPM_structure and chain A and resi 39
color red, 2MPM_structure and chain A and resi 45
color red, 2MPM_structure and chain A and resi 46
color red, 2MPM_structure and chain A and resi 47
color red, 2MPM_structure and chain A and resi 48
color red, 2MPM_structure and chain A and resi 50
color red, 2MPM_structure and chain A and resi 51
color red, 2MPM_structure and chain A and resi 69
color red, 2MPM_structure and chain A and resi 73
set cartoon_transparency, 0.2, 2MPM_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 19
