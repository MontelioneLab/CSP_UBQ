reinitialize
load ./outputs/2BN5/2BN5_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 5 residues
color tp_color, structure and chain A and resi 660
color tp_color, structure and chain A and resi 661
color tp_color, structure and chain A and resi 665
color tp_color, structure and chain A and resi 674
color tp_color, structure and chain A and resi 677
# FP (Sig. CSP -- Allosteric): 18 residues
color fp_color, structure and chain A and resi 653
color fp_color, structure and chain A and resi 654
color fp_color, structure and chain A and resi 655
color fp_color, structure and chain A and resi 656
color fp_color, structure and chain A and resi 657
color fp_color, structure and chain A and resi 658
color fp_color, structure and chain A and resi 662
color fp_color, structure and chain A and resi 664
color fp_color, structure and chain A and resi 666
color fp_color, structure and chain A and resi 667
color fp_color, structure and chain A and resi 668
color fp_color, structure and chain A and resi 670
color fp_color, structure and chain A and resi 671
color fp_color, structure and chain A and resi 672
color fp_color, structure and chain A and resi 673
color fp_color, structure and chain A and resi 675
color fp_color, structure and chain A and resi 676
color fp_color, structure and chain A and resi 678
# TN (low CSP -- Allosteric): 9 residues
color tn_color, structure and chain A and resi 652
color tn_color, structure and chain A and resi 659
color tn_color, structure and chain A and resi 663
color tn_color, structure and chain A and resi 669
color tn_color, structure and chain A and resi 679
color tn_color, structure and chain A and resi 680
color tn_color, structure and chain A and resi 681
color tn_color, structure and chain A and resi 682
color tn_color, structure and chain A and resi 683
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
# TP (Sig. CSP in Union Site): 5
# FP (Sig. CSP -- Allosteric): 18
# TN (low CSP -- Allosteric): 9
# FN (low CSP in Union Site): 0
# Residues without CSP data: 1
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2BN5/2BN5_csp.pdb
