reinitialize
load ./outputs/2MMA_19849/2MMA_csp.pdb, 2MMA_structure
hide everything, 2MMA_structure
show cartoon, 2MMA_structure
color cyan, 2MMA_structure and chain A
color gray30, 2MMA_structure and chain B
color red, 2MMA_structure and chain B and resi 9
color red, 2MMA_structure and chain B and resi 10
color red, 2MMA_structure and chain B and resi 19
color red, 2MMA_structure and chain B and resi 20
color red, 2MMA_structure and chain B and resi 47
color red, 2MMA_structure and chain B and resi 49
color red, 2MMA_structure and chain B and resi 53
color red, 2MMA_structure and chain B and resi 54
color red, 2MMA_structure and chain B and resi 61
color red, 2MMA_structure and chain B and resi 66
color red, 2MMA_structure and chain B and resi 77
color red, 2MMA_structure and chain B and resi 78
set cartoon_transparency, 0.2, 2MMA_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 12
