reinitialize
load ./outputs/5I22/5I22_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 12 residues
color tp_color, structure and chain A and resi 529
color tp_color, structure and chain A and resi 535
color tp_color, structure and chain A and resi 537
color tp_color, structure and chain A and resi 538
color tp_color, structure and chain A and resi 556
color tp_color, structure and chain A and resi 557
color tp_color, structure and chain A and resi 559
color tp_color, structure and chain A and resi 560
color tp_color, structure and chain A and resi 562
color tp_color, structure and chain A and resi 583
color tp_color, structure and chain A and resi 587
color tp_color, structure and chain A and resi 588
# FP (Sig. CSP -- Allosteric): 17 residues
color fp_color, structure and chain A and resi 515
color fp_color, structure and chain A and resi 516
color fp_color, structure and chain A and resi 530
color fp_color, structure and chain A and resi 531
color fp_color, structure and chain A and resi 532
color fp_color, structure and chain A and resi 533
color fp_color, structure and chain A and resi 534
color fp_color, structure and chain A and resi 536
color fp_color, structure and chain A and resi 539
color fp_color, structure and chain A and resi 540
color fp_color, structure and chain A and resi 546
color fp_color, structure and chain A and resi 548
color fp_color, structure and chain A and resi 553
color fp_color, structure and chain A and resi 558
color fp_color, structure and chain A and resi 578
color fp_color, structure and chain A and resi 579
color fp_color, structure and chain A and resi 584
# TN (low CSP -- Allosteric): 42 residues
color tn_color, structure and chain A and resi 517
color tn_color, structure and chain A and resi 520
color tn_color, structure and chain A and resi 521
color tn_color, structure and chain A and resi 522
color tn_color, structure and chain A and resi 523
color tn_color, structure and chain A and resi 524
color tn_color, structure and chain A and resi 525
color tn_color, structure and chain A and resi 526
color tn_color, structure and chain A and resi 527
color tn_color, structure and chain A and resi 528
color tn_color, structure and chain A and resi 541
color tn_color, structure and chain A and resi 542
color tn_color, structure and chain A and resi 543
color tn_color, structure and chain A and resi 544
color tn_color, structure and chain A and resi 545
color tn_color, structure and chain A and resi 547
color tn_color, structure and chain A and resi 549
color tn_color, structure and chain A and resi 550
color tn_color, structure and chain A and resi 552
color tn_color, structure and chain A and resi 554
color tn_color, structure and chain A and resi 561
color tn_color, structure and chain A and resi 563
color tn_color, structure and chain A and resi 564
color tn_color, structure and chain A and resi 565
color tn_color, structure and chain A and resi 566
color tn_color, structure and chain A and resi 567
color tn_color, structure and chain A and resi 568
color tn_color, structure and chain A and resi 569
color tn_color, structure and chain A and resi 570
color tn_color, structure and chain A and resi 571
color tn_color, structure and chain A and resi 572
color tn_color, structure and chain A and resi 573
color tn_color, structure and chain A and resi 574
color tn_color, structure and chain A and resi 577
color tn_color, structure and chain A and resi 580
color tn_color, structure and chain A and resi 581
color tn_color, structure and chain A and resi 582
color tn_color, structure and chain A and resi 586
color tn_color, structure and chain A and resi 589
color tn_color, structure and chain A and resi 590
color tn_color, structure and chain A and resi 591
color tn_color, structure and chain A and resi 592
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
# TP (Sig. CSP in Union Site): 12
# FP (Sig. CSP -- Allosteric): 17
# TN (low CSP -- Allosteric): 42
# FN (low CSP in Union Site): 0
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/5I22/5I22_csp.pdb
