reinitialize
load ./outputs/2KFT_6374/2KFT_csp.pdb, 2KFT_structure
hide everything, 2KFT_structure
show cartoon, 2KFT_structure
color cyan, 2KFT_structure and chain B
color gray30, 2KFT_structure and chain A
color red, 2KFT_structure and chain A and resi 297
color red, 2KFT_structure and chain A and resi 299
color red, 2KFT_structure and chain A and resi 312
color red, 2KFT_structure and chain A and resi 320
color red, 2KFT_structure and chain A and resi 346
color red, 2KFT_structure and chain A and resi 347
set cartoon_transparency, 0.2, 2KFT_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 6
