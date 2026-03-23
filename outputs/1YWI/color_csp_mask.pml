reinitialize
load ./outputs/1YWI/1YWI_csp.pdb, 1YWI_structure
hide everything, 1YWI_structure
show cartoon, 1YWI_structure
color cyan, 1YWI_structure and chain B
color gray30, 1YWI_structure and chain A
color red, 1YWI_structure and chain A and resi 18
color red, 1YWI_structure and chain A and resi 20
color red, 1YWI_structure and chain A and resi 22
color red, 1YWI_structure and chain A and resi 27
color red, 1YWI_structure and chain A and resi 28
color red, 1YWI_structure and chain A and resi 31
color red, 1YWI_structure and chain A and resi 37
color red, 1YWI_structure and chain A and resi 38
color red, 1YWI_structure and chain A and resi 39
color red, 1YWI_structure and chain A and resi 41
set cartoon_transparency, 0.2, 1YWI_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 10
