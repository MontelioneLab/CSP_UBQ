reinitialize
load ./outputs/1I5H/1I5H_csp.pdb, 1I5H_structure
hide everything, 1I5H_structure
show cartoon, 1I5H_structure
color cyan, 1I5H_structure and chain B
color gray30, 1I5H_structure and chain W
color red, 1I5H_structure and chain W and resi 472
color red, 1I5H_structure and chain W and resi 481
color red, 1I5H_structure and chain W and resi 483
color red, 1I5H_structure and chain W and resi 484
color red, 1I5H_structure and chain W and resi 485
color red, 1I5H_structure and chain W and resi 486
color red, 1I5H_structure and chain W and resi 488
color red, 1I5H_structure and chain W and resi 494
set cartoon_transparency, 0.2, 1I5H_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (W): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 8
