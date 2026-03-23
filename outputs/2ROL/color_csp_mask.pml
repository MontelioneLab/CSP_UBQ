reinitialize
load ./outputs/2ROL/2ROL_csp.pdb, 2ROL_structure
hide everything, 2ROL_structure
show cartoon, 2ROL_structure
color cyan, 2ROL_structure and chain B
color gray30, 2ROL_structure and chain A
color red, 2ROL_structure and chain A and resi 479
color red, 2ROL_structure and chain A and resi 480
color red, 2ROL_structure and chain A and resi 482
color red, 2ROL_structure and chain A and resi 484
color red, 2ROL_structure and chain A and resi 490
color red, 2ROL_structure and chain A and resi 492
color red, 2ROL_structure and chain A and resi 495
color red, 2ROL_structure and chain A and resi 496
color red, 2ROL_structure and chain A and resi 507
color red, 2ROL_structure and chain A and resi 510
color red, 2ROL_structure and chain A and resi 512
color red, 2ROL_structure and chain A and resi 514
color red, 2ROL_structure and chain A and resi 515
color red, 2ROL_structure and chain A and resi 526
color red, 2ROL_structure and chain A and resi 529
color red, 2ROL_structure and chain A and resi 530
color red, 2ROL_structure and chain A and resi 531
color red, 2ROL_structure and chain A and resi 533
color red, 2ROL_structure and chain A and resi 535
color red, 2ROL_structure and chain A and resi 537
set cartoon_transparency, 0.2, 2ROL_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 20
