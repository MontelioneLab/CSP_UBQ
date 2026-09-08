reinitialize
load ./outputs/2LAS_50128/2las_csp.pdb, 2las_structure
hide everything, 2las_structure
show cartoon, 2las_structure
color cyan, 2las_structure and chain B
color gray30, 2las_structure and chain A
color red, 2las_structure and chain A and resi 21
color red, 2las_structure and chain A and resi 27
color red, 2las_structure and chain A and resi 43
color red, 2las_structure and chain A and resi 46
color red, 2las_structure and chain A and resi 48
color red, 2las_structure and chain A and resi 49
color red, 2las_structure and chain A and resi 50
color red, 2las_structure and chain A and resi 51
color red, 2las_structure and chain A and resi 52
color red, 2las_structure and chain A and resi 53
color red, 2las_structure and chain A and resi 56
color red, 2las_structure and chain A and resi 57
color red, 2las_structure and chain A and resi 59
color red, 2las_structure and chain A and resi 60
color red, 2las_structure and chain A and resi 61
color red, 2las_structure and chain A and resi 62
color red, 2las_structure and chain A and resi 63
color red, 2las_structure and chain A and resi 71
color red, 2las_structure and chain A and resi 86
color red, 2las_structure and chain A and resi 87
color red, 2las_structure and chain A and resi 90
color red, 2las_structure and chain A and resi 93
color red, 2las_structure and chain A and resi 97
set cartoon_transparency, 0.2, 2las_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 23
