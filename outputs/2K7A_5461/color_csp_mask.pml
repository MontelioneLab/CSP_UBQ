reinitialize
load ./outputs/2K7A_5461/2K7A_csp.pdb, 2K7A_structure
hide everything, 2K7A_structure
show cartoon, 2K7A_structure
color cyan, 2K7A_structure and chain A
color gray30, 2K7A_structure and chain B
color red, 2K7A_structure and chain B and resi 280
color red, 2K7A_structure and chain B and resi 281
color red, 2K7A_structure and chain B and resi 282
color red, 2K7A_structure and chain B and resi 283
color red, 2K7A_structure and chain B and resi 285
color red, 2K7A_structure and chain B and resi 311
color red, 2K7A_structure and chain B and resi 323
color red, 2K7A_structure and chain B and resi 324
color red, 2K7A_structure and chain B and resi 332
color red, 2K7A_structure and chain B and resi 334
set cartoon_transparency, 0.2, 2K7A_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 10
