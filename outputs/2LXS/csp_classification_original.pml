reinitialize
load ./outputs/2LXS/2LXS_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 7 residues
color tp_color, structure and chain A and resi 607
color tp_color, structure and chain A and resi 611
color tp_color, structure and chain A and resi 612
color tp_color, structure and chain A and resi 614
color tp_color, structure and chain A and resi 624
color tp_color, structure and chain A and resi 631
color tp_color, structure and chain A and resi 660
# FP (Sig. CSP -- Allosteric): 25 residues
color fp_color, structure and chain A and resi 587
color fp_color, structure and chain A and resi 603
color fp_color, structure and chain A and resi 605
color fp_color, structure and chain A and resi 609
color fp_color, structure and chain A and resi 616
color fp_color, structure and chain A and resi 618
color fp_color, structure and chain A and resi 619
color fp_color, structure and chain A and resi 620
color fp_color, structure and chain A and resi 621
color fp_color, structure and chain A and resi 622
color fp_color, structure and chain A and resi 623
color fp_color, structure and chain A and resi 625
color fp_color, structure and chain A and resi 626
color fp_color, structure and chain A and resi 627
color fp_color, structure and chain A and resi 630
color fp_color, structure and chain A and resi 632
color fp_color, structure and chain A and resi 634
color fp_color, structure and chain A and resi 641
color fp_color, structure and chain A and resi 650
color fp_color, structure and chain A and resi 651
color fp_color, structure and chain A and resi 652
color fp_color, structure and chain A and resi 653
color fp_color, structure and chain A and resi 655
color fp_color, structure and chain A and resi 659
color fp_color, structure and chain A and resi 671
# TN (low CSP -- Allosteric): 46 residues
color tn_color, structure and chain A and resi 588
color tn_color, structure and chain A and resi 589
color tn_color, structure and chain A and resi 590
color tn_color, structure and chain A and resi 591
color tn_color, structure and chain A and resi 592
color tn_color, structure and chain A and resi 593
color tn_color, structure and chain A and resi 594
color tn_color, structure and chain A and resi 595
color tn_color, structure and chain A and resi 596
color tn_color, structure and chain A and resi 597
color tn_color, structure and chain A and resi 598
color tn_color, structure and chain A and resi 599
color tn_color, structure and chain A and resi 600
color tn_color, structure and chain A and resi 601
color tn_color, structure and chain A and resi 602
color tn_color, structure and chain A and resi 604
color tn_color, structure and chain A and resi 606
color tn_color, structure and chain A and resi 608
color tn_color, structure and chain A and resi 610
color tn_color, structure and chain A and resi 628
color tn_color, structure and chain A and resi 629
color tn_color, structure and chain A and resi 633
color tn_color, structure and chain A and resi 635
color tn_color, structure and chain A and resi 636
color tn_color, structure and chain A and resi 637
color tn_color, structure and chain A and resi 639
color tn_color, structure and chain A and resi 640
color tn_color, structure and chain A and resi 642
color tn_color, structure and chain A and resi 643
color tn_color, structure and chain A and resi 644
color tn_color, structure and chain A and resi 645
color tn_color, structure and chain A and resi 646
color tn_color, structure and chain A and resi 647
color tn_color, structure and chain A and resi 648
color tn_color, structure and chain A and resi 649
color tn_color, structure and chain A and resi 654
color tn_color, structure and chain A and resi 657
color tn_color, structure and chain A and resi 658
color tn_color, structure and chain A and resi 661
color tn_color, structure and chain A and resi 662
color tn_color, structure and chain A and resi 664
color tn_color, structure and chain A and resi 665
color tn_color, structure and chain A and resi 666
color tn_color, structure and chain A and resi 669
color tn_color, structure and chain A and resi 670
color tn_color, structure and chain A and resi 672
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain A and resi 638
color fn_color, structure and chain A and resi 656
color fn_color, structure and chain A and resi 663
color fn_color, structure and chain A and resi 667
color fn_color, structure and chain A and resi 668
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
# TP (Sig. CSP in Union Site): 7
# FP (Sig. CSP -- Allosteric): 25
# TN (low CSP -- Allosteric): 46
# FN (low CSP in Union Site): 5
# Residues without CSP data: 4
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2LXS/2LXS_csp.pdb
