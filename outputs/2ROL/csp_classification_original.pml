reinitialize
load ./outputs/2ROL/2ROL_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 6 residues
color tp_color, structure and chain A and resi 495
color tp_color, structure and chain A and resi 512
color tp_color, structure and chain A and resi 514
color tp_color, structure and chain A and resi 529
color tp_color, structure and chain A and resi 530
color tp_color, structure and chain A and resi 531
# FP (Sig. CSP -- Allosteric): 14 residues
color fp_color, structure and chain A and resi 479
color fp_color, structure and chain A and resi 480
color fp_color, structure and chain A and resi 482
color fp_color, structure and chain A and resi 484
color fp_color, structure and chain A and resi 490
color fp_color, structure and chain A and resi 492
color fp_color, structure and chain A and resi 496
color fp_color, structure and chain A and resi 507
color fp_color, structure and chain A and resi 510
color fp_color, structure and chain A and resi 515
color fp_color, structure and chain A and resi 526
color fp_color, structure and chain A and resi 533
color fp_color, structure and chain A and resi 535
color fp_color, structure and chain A and resi 537
# TN (low CSP -- Allosteric): 34 residues
color tn_color, structure and chain A and resi 477
color tn_color, structure and chain A and resi 481
color tn_color, structure and chain A and resi 483
color tn_color, structure and chain A and resi 485
color tn_color, structure and chain A and resi 486
color tn_color, structure and chain A and resi 487
color tn_color, structure and chain A and resi 488
color tn_color, structure and chain A and resi 489
color tn_color, structure and chain A and resi 491
color tn_color, structure and chain A and resi 493
color tn_color, structure and chain A and resi 497
color tn_color, structure and chain A and resi 498
color tn_color, structure and chain A and resi 499
color tn_color, structure and chain A and resi 500
color tn_color, structure and chain A and resi 501
color tn_color, structure and chain A and resi 502
color tn_color, structure and chain A and resi 503
color tn_color, structure and chain A and resi 504
color tn_color, structure and chain A and resi 505
color tn_color, structure and chain A and resi 506
color tn_color, structure and chain A and resi 508
color tn_color, structure and chain A and resi 509
color tn_color, structure and chain A and resi 511
color tn_color, structure and chain A and resi 516
color tn_color, structure and chain A and resi 517
color tn_color, structure and chain A and resi 518
color tn_color, structure and chain A and resi 519
color tn_color, structure and chain A and resi 521
color tn_color, structure and chain A and resi 522
color tn_color, structure and chain A and resi 523
color tn_color, structure and chain A and resi 524
color tn_color, structure and chain A and resi 525
color tn_color, structure and chain A and resi 527
color tn_color, structure and chain A and resi 532
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
# Receptor chain: A
# Ligand chain: B
# TP (Sig. CSP in Union Site): 6
# FP (Sig. CSP -- Allosteric): 14
# TN (low CSP -- Allosteric): 34
# FN (low CSP in Union Site): 0
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2ROL/2ROL_csp.pdb
