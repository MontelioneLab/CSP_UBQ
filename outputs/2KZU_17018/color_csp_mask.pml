reinitialize
load ./outputs/2KZU_17018/2KZU_csp.pdb, 2KZU_structure
hide everything, 2KZU_structure
show cartoon, 2KZU_structure
color cyan, 2KZU_structure and chain B
color gray30, 2KZU_structure and chain A
color red, 2KZU_structure and chain A and resi 85
color red, 2KZU_structure and chain A and resi 98
color red, 2KZU_structure and chain A and resi 125
color red, 2KZU_structure and chain A and resi 126
color red, 2KZU_structure and chain A and resi 128
set cartoon_transparency, 0.2, 2KZU_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 5
