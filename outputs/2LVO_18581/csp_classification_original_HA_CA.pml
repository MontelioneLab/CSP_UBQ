reinitialize
load ./outputs/2LVO_18581/2LVO_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain C
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 5 residues
color tp_color, structure and chain C and resi 466
color tp_color, structure and chain C and resi 467
color tp_color, structure and chain C and resi 486
color tp_color, structure and chain C and resi 490
color tp_color, structure and chain C and resi 494
# FP (Sig. CSP -- Allosteric): 13 residues
color fp_color, structure and chain C and resi 461
color fp_color, structure and chain C and resi 464
color fp_color, structure and chain C and resi 465
color fp_color, structure and chain C and resi 473
color fp_color, structure and chain C and resi 475
color fp_color, structure and chain C and resi 476
color fp_color, structure and chain C and resi 477
color fp_color, structure and chain C and resi 478
color fp_color, structure and chain C and resi 481
color fp_color, structure and chain C and resi 483
color fp_color, structure and chain C and resi 484
color fp_color, structure and chain C and resi 493
color fp_color, structure and chain C and resi 496
# TN (low CSP -- Allosteric): 21 residues
color tn_color, structure and chain C and resi 456
color tn_color, structure and chain C and resi 457
color tn_color, structure and chain C and resi 458
color tn_color, structure and chain C and resi 459
color tn_color, structure and chain C and resi 460
color tn_color, structure and chain C and resi 462
color tn_color, structure and chain C and resi 463
color tn_color, structure and chain C and resi 471
color tn_color, structure and chain C and resi 474
color tn_color, structure and chain C and resi 479
color tn_color, structure and chain C and resi 480
color tn_color, structure and chain C and resi 482
color tn_color, structure and chain C and resi 485
color tn_color, structure and chain C and resi 489
color tn_color, structure and chain C and resi 492
color tn_color, structure and chain C and resi 497
color tn_color, structure and chain C and resi 498
color tn_color, structure and chain C and resi 499
color tn_color, structure and chain C and resi 500
color tn_color, structure and chain C and resi 502
color tn_color, structure and chain C and resi 504
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain C and resi 470
color fn_color, structure and chain C and resi 488
color fn_color, structure and chain C and resi 491
color fn_color, structure and chain C and resi 495
set cartoon_transparency, 0.2, structure
set cartoon_fancy_helices, 1
set cartoon_ring_mode, 1
# CSP Classification Color Scheme:
# TP (#2ecc71): Sig. CSP in Union Site
# FP (#9b59b6): Sig. CSP -- Allosteric
# TN (#3498db): low CSP -- Allosteric
# FN (#f39c12): low CSP in Union Site
# Cyan: Peptide chain
# Gray: Protein residues without CSP data
# CSP Classification Analysis Summary:
# Receptor chain: C
# Ligand chain: A
# TP (Sig. CSP in Union Site): 5
# FP (Sig. CSP -- Allosteric): 13
# TN (low CSP -- Allosteric): 21
# FN (low CSP in Union Site): 4
# Residues without CSP data: 7
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/GTM_Research/01_Active_Projects/new_CSP_UBQ/CSP_UBQ/outputs/2LVO_18581/2LVO_csp.pdb
