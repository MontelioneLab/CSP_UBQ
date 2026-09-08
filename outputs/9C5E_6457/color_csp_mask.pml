reinitialize
load ./outputs/9C5E_6457/9C5E_csp.pdb, 9C5E_structure
hide everything, 9C5E_structure
show cartoon, 9C5E_structure
color cyan, 9C5E_structure and chain A
color gray30, 9C5E_structure and chain B
color red, 9C5E_structure and chain B and resi 2
color red, 9C5E_structure and chain B and resi 6
color red, 9C5E_structure and chain B and resi 25
color red, 9C5E_structure and chain B and resi 28
color red, 9C5E_structure and chain B and resi 44
color red, 9C5E_structure and chain B and resi 55
color red, 9C5E_structure and chain B and resi 68
color red, 9C5E_structure and chain B and resi 69
set cartoon_transparency, 0.2, 9C5E_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 8
