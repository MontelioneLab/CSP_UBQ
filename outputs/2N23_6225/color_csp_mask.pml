reinitialize
load ./outputs/2N23_6225/2n23_csp.pdb, 2n23_structure
hide everything, 2n23_structure
show cartoon, 2n23_structure
color cyan, 2n23_structure and chain B
color gray30, 2n23_structure and chain A
color red, 2n23_structure and chain A and resi 48
color red, 2n23_structure and chain A and resi 50
color red, 2n23_structure and chain A and resi 51
color red, 2n23_structure and chain A and resi 54
color red, 2n23_structure and chain A and resi 55
color red, 2n23_structure and chain A and resi 59
color red, 2n23_structure and chain A and resi 60
color red, 2n23_structure and chain A and resi 61
color red, 2n23_structure and chain A and resi 63
color red, 2n23_structure and chain A and resi 85
color red, 2n23_structure and chain A and resi 87
color red, 2n23_structure and chain A and resi 88
color red, 2n23_structure and chain A and resi 89
color red, 2n23_structure and chain A and resi 90
color red, 2n23_structure and chain A and resi 91
set cartoon_transparency, 0.2, 2n23_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 15
