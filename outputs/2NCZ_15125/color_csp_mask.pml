reinitialize
load ./outputs/2NCZ_15125/2ncz_csp.pdb, 2ncz_structure
hide everything, 2ncz_structure
show cartoon, 2ncz_structure
color cyan, 2ncz_structure and chain B
color gray30, 2ncz_structure and chain A
color red, 2ncz_structure and chain A and resi 2
color red, 2ncz_structure and chain A and resi 3
color red, 2ncz_structure and chain A and resi 4
color red, 2ncz_structure and chain A and resi 6
color red, 2ncz_structure and chain A and resi 12
color red, 2ncz_structure and chain A and resi 16
color red, 2ncz_structure and chain A and resi 30
color red, 2ncz_structure and chain A and resi 35
color red, 2ncz_structure and chain A and resi 46
color red, 2ncz_structure and chain A and resi 52
color red, 2ncz_structure and chain A and resi 54
color red, 2ncz_structure and chain A and resi 56
color red, 2ncz_structure and chain A and resi 82
set cartoon_transparency, 0.2, 2ncz_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 13
