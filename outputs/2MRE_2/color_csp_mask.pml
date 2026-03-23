reinitialize
load ./outputs/2MRE_2/2MRE_csp.pdb, 2MRE_structure
hide everything, 2MRE_structure
show cartoon, 2MRE_structure
color cyan, 2MRE_structure and chain A
color gray30, 2MRE_structure and chain B
color red, 2MRE_structure and chain B and resi 201
color red, 2MRE_structure and chain B and resi 202
color red, 2MRE_structure and chain B and resi 211
color red, 2MRE_structure and chain B and resi 213
color red, 2MRE_structure and chain B and resi 217
color red, 2MRE_structure and chain B and resi 218
color red, 2MRE_structure and chain B and resi 220
color red, 2MRE_structure and chain B and resi 221
color red, 2MRE_structure and chain B and resi 223
color red, 2MRE_structure and chain B and resi 224
set cartoon_transparency, 0.2, 2MRE_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (B): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 10
