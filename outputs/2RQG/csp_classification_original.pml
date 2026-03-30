reinitialize
load ./outputs/2RQG/2RQG_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 10 residues
color tp_color, structure and chain B and resi 564
color tp_color, structure and chain B and resi 583
color tp_color, structure and chain B and resi 584
color tp_color, structure and chain B and resi 586
color tp_color, structure and chain B and resi 587
color tp_color, structure and chain B and resi 588
color tp_color, structure and chain B and resi 606
color tp_color, structure and chain B and resi 609
color tp_color, structure and chain B and resi 613
color tp_color, structure and chain B and resi 614
# FP (Sig. CSP -- Allosteric): 24 residues
color fp_color, structure and chain B and resi 541
color fp_color, structure and chain B and resi 561
color fp_color, structure and chain B and resi 563
color fp_color, structure and chain B and resi 566
color fp_color, structure and chain B and resi 567
color fp_color, structure and chain B and resi 571
color fp_color, structure and chain B and resi 577
color fp_color, structure and chain B and resi 580
color fp_color, structure and chain B and resi 581
color fp_color, structure and chain B and resi 582
color fp_color, structure and chain B and resi 585
color fp_color, structure and chain B and resi 589
color fp_color, structure and chain B and resi 593
color fp_color, structure and chain B and resi 596
color fp_color, structure and chain B and resi 607
color fp_color, structure and chain B and resi 610
color fp_color, structure and chain B and resi 612
color fp_color, structure and chain B and resi 616
color fp_color, structure and chain B and resi 617
color fp_color, structure and chain B and resi 618
color fp_color, structure and chain B and resi 619
color fp_color, structure and chain B and resi 620
color fp_color, structure and chain B and resi 622
color fp_color, structure and chain B and resi 623
# TN (low CSP -- Allosteric): 38 residues
color tn_color, structure and chain B and resi 542
color tn_color, structure and chain B and resi 543
color tn_color, structure and chain B and resi 545
color tn_color, structure and chain B and resi 546
color tn_color, structure and chain B and resi 547
color tn_color, structure and chain B and resi 549
color tn_color, structure and chain B and resi 550
color tn_color, structure and chain B and resi 551
color tn_color, structure and chain B and resi 552
color tn_color, structure and chain B and resi 553
color tn_color, structure and chain B and resi 556
color tn_color, structure and chain B and resi 557
color tn_color, structure and chain B and resi 558
color tn_color, structure and chain B and resi 562
color tn_color, structure and chain B and resi 565
color tn_color, structure and chain B and resi 569
color tn_color, structure and chain B and resi 570
color tn_color, structure and chain B and resi 572
color tn_color, structure and chain B and resi 573
color tn_color, structure and chain B and resi 574
color tn_color, structure and chain B and resi 576
color tn_color, structure and chain B and resi 578
color tn_color, structure and chain B and resi 590
color tn_color, structure and chain B and resi 591
color tn_color, structure and chain B and resi 594
color tn_color, structure and chain B and resi 595
color tn_color, structure and chain B and resi 597
color tn_color, structure and chain B and resi 598
color tn_color, structure and chain B and resi 599
color tn_color, structure and chain B and resi 601
color tn_color, structure and chain B and resi 602
color tn_color, structure and chain B and resi 603
color tn_color, structure and chain B and resi 604
color tn_color, structure and chain B and resi 605
color tn_color, structure and chain B and resi 608
color tn_color, structure and chain B and resi 611
color tn_color, structure and chain B and resi 615
color tn_color, structure and chain B and resi 621
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain B and resi 559
color fn_color, structure and chain B and resi 560
color fn_color, structure and chain B and resi 592
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
# Receptor chain: B
# Ligand chain: A
# TP (Sig. CSP in Union Site): 10
# FP (Sig. CSP -- Allosteric): 24
# TN (low CSP -- Allosteric): 38
# FN (low CSP in Union Site): 3
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2RQG/2RQG_csp.pdb
