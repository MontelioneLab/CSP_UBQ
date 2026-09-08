reinitialize
load ./outputs/2LTO_10276/2LTO_csp.pdb, 2LTO_structure
hide everything, 2LTO_structure
show cartoon, 2LTO_structure
color cyan, 2LTO_structure and chain B
color gray30, 2LTO_structure and chain A
color red, 2LTO_structure and chain A and resi 554
color red, 2LTO_structure and chain A and resi 556
color red, 2LTO_structure and chain A and resi 566
color red, 2LTO_structure and chain A and resi 568
color red, 2LTO_structure and chain A and resi 573
color red, 2LTO_structure and chain A and resi 574
color red, 2LTO_structure and chain A and resi 575
color red, 2LTO_structure and chain A and resi 576
color red, 2LTO_structure and chain A and resi 591
color red, 2LTO_structure and chain A and resi 593
color red, 2LTO_structure and chain A and resi 595
color red, 2LTO_structure and chain A and resi 596
color red, 2LTO_structure and chain A and resi 598
set cartoon_transparency, 0.2, 2LTO_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 13
