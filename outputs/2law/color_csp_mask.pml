reinitialize
load ./outputs/2law/2law_csp.pdb, 2law_structure
hide everything, 2law_structure
show cartoon, 2law_structure
color cyan, 2law_structure and chain B
color gray30, 2law_structure and chain A
color red, 2law_structure and chain A and resi 239
color red, 2law_structure and chain A and resi 243
color red, 2law_structure and chain A and resi 247
color red, 2law_structure and chain A and resi 248
color red, 2law_structure and chain A and resi 249
color red, 2law_structure and chain A and resi 250
color red, 2law_structure and chain A and resi 251
color red, 2law_structure and chain A and resi 254
color red, 2law_structure and chain A and resi 255
color red, 2law_structure and chain A and resi 259
color red, 2law_structure and chain A and resi 262
set cartoon_transparency, 0.2, 2law_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 11
