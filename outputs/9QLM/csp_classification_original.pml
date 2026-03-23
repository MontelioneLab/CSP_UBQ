reinitialize
load ./outputs/9QLM/9QLM_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 16 residues
color tp_color, structure and chain A and resi 868
color tp_color, structure and chain A and resi 870
color tp_color, structure and chain A and resi 877
color tp_color, structure and chain A and resi 878
color tp_color, structure and chain A and resi 879
color tp_color, structure and chain A and resi 880
color tp_color, structure and chain A and resi 882
color tp_color, structure and chain A and resi 883
color tp_color, structure and chain A and resi 886
color tp_color, structure and chain A and resi 891
color tp_color, structure and chain A and resi 893
color tp_color, structure and chain A and resi 900
color tp_color, structure and chain A and resi 902
color tp_color, structure and chain A and resi 905
color tp_color, structure and chain A and resi 908
color tp_color, structure and chain A and resi 909
# FP (Sig. CSP -- Allosteric): 11 residues
color fp_color, structure and chain A and resi 857
color fp_color, structure and chain A and resi 858
color fp_color, structure and chain A and resi 859
color fp_color, structure and chain A and resi 860
color fp_color, structure and chain A and resi 862
color fp_color, structure and chain A and resi 865
color fp_color, structure and chain A and resi 867
color fp_color, structure and chain A and resi 869
color fp_color, structure and chain A and resi 889
color fp_color, structure and chain A and resi 890
color fp_color, structure and chain A and resi 906
# TN (low CSP -- Allosteric): 28 residues
color tn_color, structure and chain A and resi 854
color tn_color, structure and chain A and resi 855
color tn_color, structure and chain A and resi 856
color tn_color, structure and chain A and resi 861
color tn_color, structure and chain A and resi 863
color tn_color, structure and chain A and resi 864
color tn_color, structure and chain A and resi 866
color tn_color, structure and chain A and resi 872
color tn_color, structure and chain A and resi 874
color tn_color, structure and chain A and resi 875
color tn_color, structure and chain A and resi 887
color tn_color, structure and chain A and resi 892
color tn_color, structure and chain A and resi 897
color tn_color, structure and chain A and resi 898
color tn_color, structure and chain A and resi 899
color tn_color, structure and chain A and resi 901
color tn_color, structure and chain A and resi 910
color tn_color, structure and chain A and resi 913
color tn_color, structure and chain A and resi 915
color tn_color, structure and chain A and resi 916
color tn_color, structure and chain A and resi 917
color tn_color, structure and chain A and resi 918
color tn_color, structure and chain A and resi 919
color tn_color, structure and chain A and resi 920
color tn_color, structure and chain A and resi 921
color tn_color, structure and chain A and resi 922
color tn_color, structure and chain A and resi 923
color tn_color, structure and chain A and resi 924
# FN (low CSP in Binding Site): 8 residues
color fn_color, structure and chain A and resi 873
color fn_color, structure and chain A and resi 884
color fn_color, structure and chain A and resi 888
color fn_color, structure and chain A and resi 894
color fn_color, structure and chain A and resi 896
color fn_color, structure and chain A and resi 907
color fn_color, structure and chain A and resi 911
color fn_color, structure and chain A and resi 914
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
# TP (Sig. CSP in Union Site): 16
# FP (Sig. CSP -- Allosteric): 11
# TN (low CSP -- Allosteric): 28
# FN (low CSP in Union Site): 8
# Residues without CSP data: 12
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/9QLM/9QLM_csp.pdb
