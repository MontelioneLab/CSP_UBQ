reinitialize
load ./outputs/2lay/2lay_csp.pdb, 2lay_structure
hide everything, 2lay_structure
show cartoon, 2lay_structure
color cyan, 2lay_structure and chain B
color gray30, 2lay_structure and chain A
color red, 2lay_structure and chain A and resi 177
color red, 2lay_structure and chain A and resi 186
color red, 2lay_structure and chain A and resi 187
color red, 2lay_structure and chain A and resi 190
color red, 2lay_structure and chain A and resi 191
color red, 2lay_structure and chain A and resi 194
color red, 2lay_structure and chain A and resi 195
color red, 2lay_structure and chain A and resi 196
color red, 2lay_structure and chain A and resi 197
color red, 2lay_structure and chain A and resi 198
color red, 2lay_structure and chain A and resi 200
color red, 2lay_structure and chain A and resi 201
color red, 2lay_structure and chain A and resi 203
color red, 2lay_structure and chain A and resi 205
set cartoon_transparency, 0.2, 2lay_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 14
