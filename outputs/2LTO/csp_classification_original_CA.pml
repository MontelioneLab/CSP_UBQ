reinitialize
load ./outputs/2LTO/2LTO_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 2 residues
color tp_color, structure and chain A and resi 591
color tp_color, structure and chain A and resi 596
# FP (Sig. CSP -- Allosteric): 15 residues
color fp_color, structure and chain A and resi 556
color fp_color, structure and chain A and resi 560
color fp_color, structure and chain A and resi 563
color fp_color, structure and chain A and resi 568
color fp_color, structure and chain A and resi 573
color fp_color, structure and chain A and resi 574
color fp_color, structure and chain A and resi 575
color fp_color, structure and chain A and resi 576
color fp_color, structure and chain A and resi 577
color fp_color, structure and chain A and resi 586
color fp_color, structure and chain A and resi 592
color fp_color, structure and chain A and resi 595
color fp_color, structure and chain A and resi 598
color fp_color, structure and chain A and resi 608
color fp_color, structure and chain A and resi 609
# TN (low CSP -- Allosteric): 29 residues
color tn_color, structure and chain A and resi 557
color tn_color, structure and chain A and resi 559
color tn_color, structure and chain A and resi 561
color tn_color, structure and chain A and resi 562
color tn_color, structure and chain A and resi 564
color tn_color, structure and chain A and resi 565
color tn_color, structure and chain A and resi 566
color tn_color, structure and chain A and resi 567
color tn_color, structure and chain A and resi 569
color tn_color, structure and chain A and resi 570
color tn_color, structure and chain A and resi 571
color tn_color, structure and chain A and resi 572
color tn_color, structure and chain A and resi 578
color tn_color, structure and chain A and resi 579
color tn_color, structure and chain A and resi 580
color tn_color, structure and chain A and resi 584
color tn_color, structure and chain A and resi 585
color tn_color, structure and chain A and resi 587
color tn_color, structure and chain A and resi 588
color tn_color, structure and chain A and resi 589
color tn_color, structure and chain A and resi 590
color tn_color, structure and chain A and resi 597
color tn_color, structure and chain A and resi 599
color tn_color, structure and chain A and resi 600
color tn_color, structure and chain A and resi 601
color tn_color, structure and chain A and resi 602
color tn_color, structure and chain A and resi 603
color tn_color, structure and chain A and resi 604
color tn_color, structure and chain A and resi 605
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
# TP (Sig. CSP in Union Site): 2
# FP (Sig. CSP -- Allosteric): 15
# TN (low CSP -- Allosteric): 29
# FN (low CSP in Union Site): 0
# Residues without CSP data: 12
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2LTO/2LTO_csp.pdb
