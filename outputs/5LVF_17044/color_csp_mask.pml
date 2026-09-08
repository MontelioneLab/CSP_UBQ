reinitialize
load ./outputs/5LVF_17044/5lvf_csp.pdb, 5lvf_structure
hide everything, 5lvf_structure
show cartoon, 5lvf_structure
color cyan, 5lvf_structure and chain B
color gray30, 5lvf_structure and chain A
color red, 5lvf_structure and chain A and resi 16
color red, 5lvf_structure and chain A and resi 22
color red, 5lvf_structure and chain A and resi 24
color red, 5lvf_structure and chain A and resi 26
color red, 5lvf_structure and chain A and resi 27
color red, 5lvf_structure and chain A and resi 59
color red, 5lvf_structure and chain A and resi 67
color red, 5lvf_structure and chain A and resi 69
color red, 5lvf_structure and chain A and resi 71
color red, 5lvf_structure and chain A and resi 75
color red, 5lvf_structure and chain A and resi 109
color red, 5lvf_structure and chain A and resi 113
color red, 5lvf_structure and chain A and resi 116
color red, 5lvf_structure and chain A and resi 117
color red, 5lvf_structure and chain A and resi 118
color red, 5lvf_structure and chain A and resi 135
color red, 5lvf_structure and chain A and resi 136
set cartoon_transparency, 0.2, 5lvf_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 17
