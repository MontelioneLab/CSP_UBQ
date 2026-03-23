reinitialize
load ./outputs/2RSE_1/2RSE_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# FP (Sig. CSP -- Allosteric): 48 residues
color fp_color, structure and chain A and resi 903
color fp_color, structure and chain A and resi 904
color fp_color, structure and chain A and resi 911
color fp_color, structure and chain A and resi 913
color fp_color, structure and chain A and resi 916
color fp_color, structure and chain A and resi 924
color fp_color, structure and chain A and resi 926
color fp_color, structure and chain A and resi 927
color fp_color, structure and chain A and resi 928
color fp_color, structure and chain A and resi 932
color fp_color, structure and chain A and resi 933
color fp_color, structure and chain A and resi 936
color fp_color, structure and chain A and resi 937
color fp_color, structure and chain A and resi 940
color fp_color, structure and chain A and resi 942
color fp_color, structure and chain A and resi 944
color fp_color, structure and chain A and resi 945
color fp_color, structure and chain A and resi 949
color fp_color, structure and chain A and resi 952
color fp_color, structure and chain A and resi 953
color fp_color, structure and chain A and resi 954
color fp_color, structure and chain A and resi 955
color fp_color, structure and chain A and resi 956
color fp_color, structure and chain A and resi 959
color fp_color, structure and chain A and resi 960
color fp_color, structure and chain A and resi 962
color fp_color, structure and chain A and resi 966
color fp_color, structure and chain A and resi 967
color fp_color, structure and chain A and resi 968
color fp_color, structure and chain A and resi 969
color fp_color, structure and chain A and resi 974
color fp_color, structure and chain A and resi 975
color fp_color, structure and chain A and resi 976
color fp_color, structure and chain A and resi 977
color fp_color, structure and chain A and resi 978
color fp_color, structure and chain A and resi 980
color fp_color, structure and chain A and resi 981
color fp_color, structure and chain A and resi 991
color fp_color, structure and chain A and resi 995
color fp_color, structure and chain A and resi 996
color fp_color, structure and chain A and resi 997
color fp_color, structure and chain A and resi 998
color fp_color, structure and chain A and resi 999
color fp_color, structure and chain A and resi 1000
color fp_color, structure and chain A and resi 1001
color fp_color, structure and chain A and resi 1004
color fp_color, structure and chain A and resi 1005
color fp_color, structure and chain A and resi 1007
# TN (low CSP -- Allosteric): 46 residues
color tn_color, structure and chain A and resi 905
color tn_color, structure and chain A and resi 906
color tn_color, structure and chain A and resi 908
color tn_color, structure and chain A and resi 909
color tn_color, structure and chain A and resi 912
color tn_color, structure and chain A and resi 914
color tn_color, structure and chain A and resi 915
color tn_color, structure and chain A and resi 918
color tn_color, structure and chain A and resi 919
color tn_color, structure and chain A and resi 920
color tn_color, structure and chain A and resi 921
color tn_color, structure and chain A and resi 922
color tn_color, structure and chain A and resi 923
color tn_color, structure and chain A and resi 925
color tn_color, structure and chain A and resi 929
color tn_color, structure and chain A and resi 930
color tn_color, structure and chain A and resi 931
color tn_color, structure and chain A and resi 934
color tn_color, structure and chain A and resi 935
color tn_color, structure and chain A and resi 938
color tn_color, structure and chain A and resi 939
color tn_color, structure and chain A and resi 941
color tn_color, structure and chain A and resi 943
color tn_color, structure and chain A and resi 950
color tn_color, structure and chain A and resi 951
color tn_color, structure and chain A and resi 957
color tn_color, structure and chain A and resi 958
color tn_color, structure and chain A and resi 961
color tn_color, structure and chain A and resi 963
color tn_color, structure and chain A and resi 964
color tn_color, structure and chain A and resi 965
color tn_color, structure and chain A and resi 970
color tn_color, structure and chain A and resi 971
color tn_color, structure and chain A and resi 972
color tn_color, structure and chain A and resi 973
color tn_color, structure and chain A and resi 982
color tn_color, structure and chain A and resi 983
color tn_color, structure and chain A and resi 984
color tn_color, structure and chain A and resi 986
color tn_color, structure and chain A and resi 987
color tn_color, structure and chain A and resi 988
color tn_color, structure and chain A and resi 992
color tn_color, structure and chain A and resi 1002
color tn_color, structure and chain A and resi 1003
color tn_color, structure and chain A and resi 1006
color tn_color, structure and chain A and resi 1008
# FN (low CSP in Binding Site): 2 residues
color fn_color, structure and chain A and resi 947
color fn_color, structure and chain A and resi 948
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
# TP (Sig. CSP in Union Site): 0
# FP (Sig. CSP -- Allosteric): 48
# TN (low CSP -- Allosteric): 46
# FN (low CSP in Union Site): 2
# Residues without CSP data: 11
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2RSE_1/2RSE_csp.pdb
