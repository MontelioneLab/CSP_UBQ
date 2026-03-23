reinitialize
load ./outputs/2KOH/2KOH_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 10 residues
color tp_color, structure and chain A and resi 600
color tp_color, structure and chain A and resi 601
color tp_color, structure and chain A and resi 602
color tp_color, structure and chain A and resi 603
color tp_color, structure and chain A and resi 604
color tp_color, structure and chain A and resi 605
color tp_color, structure and chain A and resi 606
color tp_color, structure and chain A and resi 607
color tp_color, structure and chain A and resi 611
color tp_color, structure and chain A and resi 622
# FP (Sig. CSP -- Allosteric): 41 residues
color fp_color, structure and chain A and resi 582
color fp_color, structure and chain A and resi 583
color fp_color, structure and chain A and resi 589
color fp_color, structure and chain A and resi 590
color fp_color, structure and chain A and resi 593
color fp_color, structure and chain A and resi 594
color fp_color, structure and chain A and resi 595
color fp_color, structure and chain A and resi 599
color fp_color, structure and chain A and resi 614
color fp_color, structure and chain A and resi 618
color fp_color, structure and chain A and resi 619
color fp_color, structure and chain A and resi 620
color fp_color, structure and chain A and resi 623
color fp_color, structure and chain A and resi 628
color fp_color, structure and chain A and resi 629
color fp_color, structure and chain A and resi 630
color fp_color, structure and chain A and resi 634
color fp_color, structure and chain A and resi 639
color fp_color, structure and chain A and resi 642
color fp_color, structure and chain A and resi 645
color fp_color, structure and chain A and resi 646
color fp_color, structure and chain A and resi 647
color fp_color, structure and chain A and resi 648
color fp_color, structure and chain A and resi 654
color fp_color, structure and chain A and resi 655
color fp_color, structure and chain A and resi 656
color fp_color, structure and chain A and resi 657
color fp_color, structure and chain A and resi 658
color fp_color, structure and chain A and resi 660
color fp_color, structure and chain A and resi 661
color fp_color, structure and chain A and resi 664
color fp_color, structure and chain A and resi 665
color fp_color, structure and chain A and resi 668
color fp_color, structure and chain A and resi 672
color fp_color, structure and chain A and resi 673
color fp_color, structure and chain A and resi 677
color fp_color, structure and chain A and resi 678
color fp_color, structure and chain A and resi 682
color fp_color, structure and chain A and resi 683
color fp_color, structure and chain A and resi 684
color fp_color, structure and chain A and resi 685
# TN (low CSP -- Allosteric): 35 residues
color tn_color, structure and chain A and resi 584
color tn_color, structure and chain A and resi 585
color tn_color, structure and chain A and resi 586
color tn_color, structure and chain A and resi 587
color tn_color, structure and chain A and resi 588
color tn_color, structure and chain A and resi 591
color tn_color, structure and chain A and resi 596
color tn_color, structure and chain A and resi 597
color tn_color, structure and chain A and resi 610
color tn_color, structure and chain A and resi 615
color tn_color, structure and chain A and resi 616
color tn_color, structure and chain A and resi 617
color tn_color, structure and chain A and resi 621
color tn_color, structure and chain A and resi 624
color tn_color, structure and chain A and resi 631
color tn_color, structure and chain A and resi 632
color tn_color, structure and chain A and resi 633
color tn_color, structure and chain A and resi 635
color tn_color, structure and chain A and resi 636
color tn_color, structure and chain A and resi 637
color tn_color, structure and chain A and resi 638
color tn_color, structure and chain A and resi 640
color tn_color, structure and chain A and resi 641
color tn_color, structure and chain A and resi 649
color tn_color, structure and chain A and resi 653
color tn_color, structure and chain A and resi 667
color tn_color, structure and chain A and resi 669
color tn_color, structure and chain A and resi 670
color tn_color, structure and chain A and resi 671
color tn_color, structure and chain A and resi 674
color tn_color, structure and chain A and resi 675
color tn_color, structure and chain A and resi 676
color tn_color, structure and chain A and resi 679
color tn_color, structure and chain A and resi 680
color tn_color, structure and chain A and resi 681
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain A and resi 608
color fn_color, structure and chain A and resi 625
color fn_color, structure and chain A and resi 659
color fn_color, structure and chain A and resi 662
color fn_color, structure and chain A and resi 663
color fn_color, structure and chain A and resi 666
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
# TP (Sig. CSP in Union Site): 10
# FP (Sig. CSP -- Allosteric): 41
# TN (low CSP -- Allosteric): 35
# FN (low CSP in Union Site): 6
# Residues without CSP data: 19
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2KOH/2KOH_csp.pdb
