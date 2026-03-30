reinitialize
load ./outputs/6bgg/6bgg_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 5 residues
color tp_color, structure and chain B and resi 578
color tp_color, structure and chain B and resi 592
color tp_color, structure and chain B and resi 615
color tp_color, structure and chain B and resi 616
color tp_color, structure and chain B and resi 617
# FP (Sig. CSP -- Allosteric): 34 residues
color fp_color, structure and chain B and resi 558
color fp_color, structure and chain B and resi 560
color fp_color, structure and chain B and resi 561
color fp_color, structure and chain B and resi 562
color fp_color, structure and chain B and resi 563
color fp_color, structure and chain B and resi 568
color fp_color, structure and chain B and resi 569
color fp_color, structure and chain B and resi 570
color fp_color, structure and chain B and resi 575
color fp_color, structure and chain B and resi 583
color fp_color, structure and chain B and resi 584
color fp_color, structure and chain B and resi 585
color fp_color, structure and chain B and resi 586
color fp_color, structure and chain B and resi 587
color fp_color, structure and chain B and resi 591
color fp_color, structure and chain B and resi 593
color fp_color, structure and chain B and resi 594
color fp_color, structure and chain B and resi 595
color fp_color, structure and chain B and resi 596
color fp_color, structure and chain B and resi 597
color fp_color, structure and chain B and resi 601
color fp_color, structure and chain B and resi 607
color fp_color, structure and chain B and resi 608
color fp_color, structure and chain B and resi 610
color fp_color, structure and chain B and resi 621
color fp_color, structure and chain B and resi 624
color fp_color, structure and chain B and resi 625
color fp_color, structure and chain B and resi 628
color fp_color, structure and chain B and resi 631
color fp_color, structure and chain B and resi 634
color fp_color, structure and chain B and resi 635
color fp_color, structure and chain B and resi 639
color fp_color, structure and chain B and resi 642
color fp_color, structure and chain B and resi 643
# TN (low CSP -- Allosteric): 31 residues
color tn_color, structure and chain B and resi 565
color tn_color, structure and chain B and resi 566
color tn_color, structure and chain B and resi 572
color tn_color, structure and chain B and resi 573
color tn_color, structure and chain B and resi 574
color tn_color, structure and chain B and resi 576
color tn_color, structure and chain B and resi 579
color tn_color, structure and chain B and resi 580
color tn_color, structure and chain B and resi 582
color tn_color, structure and chain B and resi 589
color tn_color, structure and chain B and resi 590
color tn_color, structure and chain B and resi 598
color tn_color, structure and chain B and resi 599
color tn_color, structure and chain B and resi 600
color tn_color, structure and chain B and resi 602
color tn_color, structure and chain B and resi 603
color tn_color, structure and chain B and resi 605
color tn_color, structure and chain B and resi 606
color tn_color, structure and chain B and resi 609
color tn_color, structure and chain B and resi 612
color tn_color, structure and chain B and resi 620
color tn_color, structure and chain B and resi 622
color tn_color, structure and chain B and resi 626
color tn_color, structure and chain B and resi 627
color tn_color, structure and chain B and resi 629
color tn_color, structure and chain B and resi 630
color tn_color, structure and chain B and resi 632
color tn_color, structure and chain B and resi 633
color tn_color, structure and chain B and resi 636
color tn_color, structure and chain B and resi 637
color tn_color, structure and chain B and resi 638
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain B and resi 577
color fn_color, structure and chain B and resi 581
color fn_color, structure and chain B and resi 613
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
# TP (Sig. CSP in Union Site): 5
# FP (Sig. CSP -- Allosteric): 34
# TN (low CSP -- Allosteric): 31
# FN (low CSP in Union Site): 3
# Residues without CSP data: 15
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/6bgg/6bgg_csp.pdb
