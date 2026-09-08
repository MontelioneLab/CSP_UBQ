reinitialize
load ./outputs/2PLD_5318/2PLD_csp.pdb, 2PLD_structure
hide everything, 2PLD_structure
show cartoon, 2PLD_structure
color cyan, 2PLD_structure and chain B
color gray30, 2PLD_structure and chain A
color red, 2PLD_structure and chain A and resi 5
color red, 2PLD_structure and chain A and resi 18
color red, 2PLD_structure and chain A and resi 44
color red, 2PLD_structure and chain A and resi 46
color red, 2PLD_structure and chain A and resi 47
color red, 2PLD_structure and chain A and resi 50
color red, 2PLD_structure and chain A and resi 51
color red, 2PLD_structure and chain A and resi 56
color red, 2PLD_structure and chain A and resi 57
color red, 2PLD_structure and chain A and resi 58
color red, 2PLD_structure and chain A and resi 59
color red, 2PLD_structure and chain A and resi 61
color red, 2PLD_structure and chain A and resi 69
color red, 2PLD_structure and chain A and resi 72
color red, 2PLD_structure and chain A and resi 89
color red, 2PLD_structure and chain A and resi 90
color red, 2PLD_structure and chain A and resi 91
color red, 2PLD_structure and chain A and resi 92
color red, 2PLD_structure and chain A and resi 93
color red, 2PLD_structure and chain A and resi 104
set cartoon_transparency, 0.2, 2PLD_structure
# Color scheme:
# Ligand chain (B): cyan
# Receptor chain (A): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 20
