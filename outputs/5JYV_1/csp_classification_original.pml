reinitialize
load ./outputs/5JYV_1/5JYV_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 14 residues
color tp_color, structure and chain A and resi 615
color tp_color, structure and chain A and resi 623
color tp_color, structure and chain A and resi 631
color tp_color, structure and chain A and resi 632
color tp_color, structure and chain A and resi 634
color tp_color, structure and chain A and resi 635
color tp_color, structure and chain A and resi 637
color tp_color, structure and chain A and resi 639
color tp_color, structure and chain A and resi 640
color tp_color, structure and chain A and resi 642
color tp_color, structure and chain A and resi 643
color tp_color, structure and chain A and resi 647
color tp_color, structure and chain A and resi 649
color tp_color, structure and chain A and resi 654
# FP (Sig. CSP -- Allosteric): 40 residues
color fp_color, structure and chain A and resi 616
color fp_color, structure and chain A and resi 617
color fp_color, structure and chain A and resi 618
color fp_color, structure and chain A and resi 620
color fp_color, structure and chain A and resi 621
color fp_color, structure and chain A and resi 622
color fp_color, structure and chain A and resi 624
color fp_color, structure and chain A and resi 625
color fp_color, structure and chain A and resi 626
color fp_color, structure and chain A and resi 627
color fp_color, structure and chain A and resi 633
color fp_color, structure and chain A and resi 636
color fp_color, structure and chain A and resi 645
color fp_color, structure and chain A and resi 646
color fp_color, structure and chain A and resi 648
color fp_color, structure and chain A and resi 651
color fp_color, structure and chain A and resi 652
color fp_color, structure and chain A and resi 657
color fp_color, structure and chain A and resi 658
color fp_color, structure and chain A and resi 659
color fp_color, structure and chain A and resi 662
color fp_color, structure and chain A and resi 666
color fp_color, structure and chain A and resi 667
color fp_color, structure and chain A and resi 670
color fp_color, structure and chain A and resi 673
color fp_color, structure and chain A and resi 677
color fp_color, structure and chain A and resi 680
color fp_color, structure and chain A and resi 684
color fp_color, structure and chain A and resi 689
color fp_color, structure and chain A and resi 695
color fp_color, structure and chain A and resi 699
color fp_color, structure and chain A and resi 700
color fp_color, structure and chain A and resi 702
color fp_color, structure and chain A and resi 706
color fp_color, structure and chain A and resi 708
color fp_color, structure and chain A and resi 709
color fp_color, structure and chain A and resi 713
color fp_color, structure and chain A and resi 714
color fp_color, structure and chain A and resi 718
color fp_color, structure and chain A and resi 720
# TN (low CSP -- Allosteric): 49 residues
color tn_color, structure and chain A and resi 619
color tn_color, structure and chain A and resi 629
color tn_color, structure and chain A and resi 638
color tn_color, structure and chain A and resi 655
color tn_color, structure and chain A and resi 660
color tn_color, structure and chain A and resi 661
color tn_color, structure and chain A and resi 663
color tn_color, structure and chain A and resi 665
color tn_color, structure and chain A and resi 668
color tn_color, structure and chain A and resi 669
color tn_color, structure and chain A and resi 671
color tn_color, structure and chain A and resi 672
color tn_color, structure and chain A and resi 674
color tn_color, structure and chain A and resi 675
color tn_color, structure and chain A and resi 676
color tn_color, structure and chain A and resi 678
color tn_color, structure and chain A and resi 679
color tn_color, structure and chain A and resi 681
color tn_color, structure and chain A and resi 682
color tn_color, structure and chain A and resi 685
color tn_color, structure and chain A and resi 686
color tn_color, structure and chain A and resi 687
color tn_color, structure and chain A and resi 688
color tn_color, structure and chain A and resi 690
color tn_color, structure and chain A and resi 691
color tn_color, structure and chain A and resi 692
color tn_color, structure and chain A and resi 693
color tn_color, structure and chain A and resi 694
color tn_color, structure and chain A and resi 696
color tn_color, structure and chain A and resi 697
color tn_color, structure and chain A and resi 698
color tn_color, structure and chain A and resi 701
color tn_color, structure and chain A and resi 703
color tn_color, structure and chain A and resi 704
color tn_color, structure and chain A and resi 711
color tn_color, structure and chain A and resi 712
color tn_color, structure and chain A and resi 715
color tn_color, structure and chain A and resi 716
color tn_color, structure and chain A and resi 717
color tn_color, structure and chain A and resi 719
color tn_color, structure and chain A and resi 721
color tn_color, structure and chain A and resi 722
color tn_color, structure and chain A and resi 723
color tn_color, structure and chain A and resi 724
color tn_color, structure and chain A and resi 725
color tn_color, structure and chain A and resi 726
color tn_color, structure and chain A and resi 727
color tn_color, structure and chain A and resi 728
color tn_color, structure and chain A and resi 729
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain A and resi 628
color fn_color, structure and chain A and resi 630
color fn_color, structure and chain A and resi 641
color fn_color, structure and chain A and resi 644
color fn_color, structure and chain A and resi 650
color fn_color, structure and chain A and resi 653
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
# TP (Sig. CSP in Union Site): 14
# FP (Sig. CSP -- Allosteric): 40
# TN (low CSP -- Allosteric): 49
# FN (low CSP in Union Site): 6
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/5JYV_1/5JYV_csp.pdb
