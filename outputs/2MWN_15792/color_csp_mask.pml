reinitialize
load ./outputs/2MWN_15792/2MWN_csp.pdb, 2MWN_structure
hide everything, 2MWN_structure
show cartoon, 2MWN_structure
color cyan, 2MWN_structure and chain A
color gray30, 2MWN_structure and chain B
color red, 2MWN_structure and chain B and resi 364
color red, 2MWN_structure and chain B and resi 370
color red, 2MWN_structure and chain B and resi 374
color red, 2MWN_structure and chain B and resi 382
color red, 2MWN_structure and chain B and resi 389
color red, 2MWN_structure and chain B and resi 393
color red, 2MWN_structure and chain B and resi 398
color red, 2MWN_structure and chain B and resi 400
set cartoon_transparency, 0.2, 2MWN_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 8
