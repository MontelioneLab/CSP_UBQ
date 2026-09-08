reinitialize
load ./outputs/2JMF_6262/2JMF_csp.pdb, 2JMF_structure
hide everything, 2JMF_structure
show cartoon, 2JMF_structure
color cyan, 2JMF_structure and chain B
color gray30, 2JMF_structure and chain A
color red, 2JMF_structure and chain A and resi 519
color red, 2JMF_structure and chain A and resi 520
color red, 2JMF_structure and chain A and resi 521
color red, 2JMF_structure and chain A and resi 529
color red, 2JMF_structure and chain A and resi 533
color red, 2JMF_structure and chain A and resi 553
set cartoon_transparency, 0.2, 2JMF_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 6
