reinitialize
load ./outputs/2kq0/2kq0_csp.pdb, 2kq0_structure
hide everything, 2kq0_structure
show cartoon, 2kq0_structure
color cyan, 2kq0_structure and chain B
color gray30, 2kq0_structure and chain A
color red, 2kq0_structure and chain A and resi 20
color red, 2kq0_structure and chain A and resi 22
color red, 2kq0_structure and chain A and resi 28
color red, 2kq0_structure and chain A and resi 29
color red, 2kq0_structure and chain A and resi 30
color red, 2kq0_structure and chain A and resi 34
color red, 2kq0_structure and chain A and resi 35
color red, 2kq0_structure and chain A and resi 36
color red, 2kq0_structure and chain A and resi 37
color red, 2kq0_structure and chain A and resi 38
color red, 2kq0_structure and chain A and resi 39
color red, 2kq0_structure and chain A and resi 40
set cartoon_transparency, 0.2, 2kq0_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 12
