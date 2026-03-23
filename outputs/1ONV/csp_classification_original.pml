reinitialize
load ./outputs/1ONV/1ONV_csp.pdb, structure
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
color tp_color, structure and chain A and resi 470
color tp_color, structure and chain A and resi 474
color tp_color, structure and chain A and resi 490
color tp_color, structure and chain A and resi 493
color tp_color, structure and chain A and resi 494
color tp_color, structure and chain A and resi 498
# FP (Sig. CSP -- Allosteric): 16 residues
color fp_color, structure and chain A and resi 452
color fp_color, structure and chain A and resi 454
color fp_color, structure and chain A and resi 461
color fp_color, structure and chain A and resi 472
color fp_color, structure and chain A and resi 475
color fp_color, structure and chain A and resi 476
color fp_color, structure and chain A and resi 477
color fp_color, structure and chain A and resi 481
color fp_color, structure and chain A and resi 482
color fp_color, structure and chain A and resi 483
color fp_color, structure and chain A and resi 485
color fp_color, structure and chain A and resi 491
color fp_color, structure and chain A and resi 492
color fp_color, structure and chain A and resi 495
color fp_color, structure and chain A and resi 496
color fp_color, structure and chain A and resi 512
# TN (low CSP -- Allosteric): 33 residues
color tn_color, structure and chain A and resi 453
color tn_color, structure and chain A and resi 455
color tn_color, structure and chain A and resi 456
color tn_color, structure and chain A and resi 457
color tn_color, structure and chain A and resi 458
color tn_color, structure and chain A and resi 459
color tn_color, structure and chain A and resi 460
color tn_color, structure and chain A and resi 462
color tn_color, structure and chain A and resi 463
color tn_color, structure and chain A and resi 464
color tn_color, structure and chain A and resi 465
color tn_color, structure and chain A and resi 466
color tn_color, structure and chain A and resi 468
color tn_color, structure and chain A and resi 469
color tn_color, structure and chain A and resi 473
color tn_color, structure and chain A and resi 484
color tn_color, structure and chain A and resi 487
color tn_color, structure and chain A and resi 488
color tn_color, structure and chain A and resi 489
color tn_color, structure and chain A and resi 499
color tn_color, structure and chain A and resi 500
color tn_color, structure and chain A and resi 503
color tn_color, structure and chain A and resi 505
color tn_color, structure and chain A and resi 506
color tn_color, structure and chain A and resi 507
color tn_color, structure and chain A and resi 508
color tn_color, structure and chain A and resi 509
color tn_color, structure and chain A and resi 510
color tn_color, structure and chain A and resi 511
color tn_color, structure and chain A and resi 514
color tn_color, structure and chain A and resi 515
color tn_color, structure and chain A and resi 516
color tn_color, structure and chain A and resi 517
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain A and resi 471
color fn_color, structure and chain A and resi 478
color fn_color, structure and chain A and resi 497
color fn_color, structure and chain A and resi 501
color fn_color, structure and chain A and resi 504
color fn_color, structure and chain A and resi 513
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
# FP (Sig. CSP -- Allosteric): 16
# TN (low CSP -- Allosteric): 33
# FN (low CSP in Union Site): 6
# Residues without CSP data: 6
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/1ONV/1ONV_csp.pdb
