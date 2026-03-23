reinitialize
load ./outputs/2LVO_1/2LVO_csp.pdb, 2LVO_structure
hide everything, 2LVO_structure
show cartoon, 2LVO_structure
color cyan, 2LVO_structure and chain A
color gray30, 2LVO_structure and chain C
color red, 2LVO_structure and chain C and resi 456
color red, 2LVO_structure and chain C and resi 459
color red, 2LVO_structure and chain C and resi 460
color red, 2LVO_structure and chain C and resi 461
color red, 2LVO_structure and chain C and resi 463
color red, 2LVO_structure and chain C and resi 466
color red, 2LVO_structure and chain C and resi 467
color red, 2LVO_structure and chain C and resi 468
color red, 2LVO_structure and chain C and resi 470
color red, 2LVO_structure and chain C and resi 471
color red, 2LVO_structure and chain C and resi 473
color red, 2LVO_structure and chain C and resi 487
color red, 2LVO_structure and chain C and resi 488
color red, 2LVO_structure and chain C and resi 489
color red, 2LVO_structure and chain C and resi 492
color red, 2LVO_structure and chain C and resi 494
color red, 2LVO_structure and chain C and resi 499
set cartoon_transparency, 0.2, 2LVO_structure
# Color scheme:
# Ligand chain (A): cyan
# Receptor chain (C): gray30 (non-significant/no CSP), red (significant CSP via significant)
# Significant residues: 17
