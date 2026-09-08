reinitialize
load ./outputs/5URN_6225/5urn_csp.pdb, 5urn_structure
hide everything, 5urn_structure
show cartoon, 5urn_structure
color cyan, 5urn_structure and chain B
color gray30, 5urn_structure and chain A
color red, 5urn_structure and chain A and resi 48
color red, 5urn_structure and chain A and resi 50
color red, 5urn_structure and chain A and resi 59
color red, 5urn_structure and chain A and resi 60
color red, 5urn_structure and chain A and resi 61
color red, 5urn_structure and chain A and resi 62
color red, 5urn_structure and chain A and resi 87
color red, 5urn_structure and chain A and resi 88
color red, 5urn_structure and chain A and resi 89
color red, 5urn_structure and chain A and resi 90
color red, 5urn_structure and chain A and resi 91
set cartoon_transparency, 0.2, 5urn_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 11
