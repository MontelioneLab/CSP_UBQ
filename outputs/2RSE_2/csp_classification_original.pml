reinitialize
load ./outputs/2RSE_2/2RSE_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 2 residues
color tp_color, structure and chain B and resi 1095
color tp_color, structure and chain B and resi 1099
# FP (Sig. CSP -- Allosteric): 35 residues
color fp_color, structure and chain B and resi 1011
color fp_color, structure and chain B and resi 1012
color fp_color, structure and chain B and resi 1015
color fp_color, structure and chain B and resi 1019
color fp_color, structure and chain B and resi 1020
color fp_color, structure and chain B and resi 1021
color fp_color, structure and chain B and resi 1022
color fp_color, structure and chain B and resi 1023
color fp_color, structure and chain B and resi 1024
color fp_color, structure and chain B and resi 1025
color fp_color, structure and chain B and resi 1026
color fp_color, structure and chain B and resi 1028
color fp_color, structure and chain B and resi 1030
color fp_color, structure and chain B and resi 1031
color fp_color, structure and chain B and resi 1032
color fp_color, structure and chain B and resi 1033
color fp_color, structure and chain B and resi 1034
color fp_color, structure and chain B and resi 1036
color fp_color, structure and chain B and resi 1037
color fp_color, structure and chain B and resi 1039
color fp_color, structure and chain B and resi 1041
color fp_color, structure and chain B and resi 1042
color fp_color, structure and chain B and resi 1046
color fp_color, structure and chain B and resi 1072
color fp_color, structure and chain B and resi 1073
color fp_color, structure and chain B and resi 1075
color fp_color, structure and chain B and resi 1076
color fp_color, structure and chain B and resi 1088
color fp_color, structure and chain B and resi 1090
color fp_color, structure and chain B and resi 1091
color fp_color, structure and chain B and resi 1093
color fp_color, structure and chain B and resi 1094
color fp_color, structure and chain B and resi 1096
color fp_color, structure and chain B and resi 1097
color fp_color, structure and chain B and resi 1102
# TN (low CSP -- Allosteric): 47 residues
color tn_color, structure and chain B and resi 1009
color tn_color, structure and chain B and resi 1010
color tn_color, structure and chain B and resi 1013
color tn_color, structure and chain B and resi 1014
color tn_color, structure and chain B and resi 1016
color tn_color, structure and chain B and resi 1018
color tn_color, structure and chain B and resi 1035
color tn_color, structure and chain B and resi 1038
color tn_color, structure and chain B and resi 1040
color tn_color, structure and chain B and resi 1044
color tn_color, structure and chain B and resi 1045
color tn_color, structure and chain B and resi 1047
color tn_color, structure and chain B and resi 1048
color tn_color, structure and chain B and resi 1049
color tn_color, structure and chain B and resi 1050
color tn_color, structure and chain B and resi 1051
color tn_color, structure and chain B and resi 1053
color tn_color, structure and chain B and resi 1054
color tn_color, structure and chain B and resi 1057
color tn_color, structure and chain B and resi 1058
color tn_color, structure and chain B and resi 1060
color tn_color, structure and chain B and resi 1061
color tn_color, structure and chain B and resi 1062
color tn_color, structure and chain B and resi 1063
color tn_color, structure and chain B and resi 1064
color tn_color, structure and chain B and resi 1065
color tn_color, structure and chain B and resi 1066
color tn_color, structure and chain B and resi 1067
color tn_color, structure and chain B and resi 1068
color tn_color, structure and chain B and resi 1069
color tn_color, structure and chain B and resi 1070
color tn_color, structure and chain B and resi 1071
color tn_color, structure and chain B and resi 1074
color tn_color, structure and chain B and resi 1077
color tn_color, structure and chain B and resi 1078
color tn_color, structure and chain B and resi 1079
color tn_color, structure and chain B and resi 1080
color tn_color, structure and chain B and resi 1081
color tn_color, structure and chain B and resi 1082
color tn_color, structure and chain B and resi 1083
color tn_color, structure and chain B and resi 1084
color tn_color, structure and chain B and resi 1086
color tn_color, structure and chain B and resi 1087
color tn_color, structure and chain B and resi 1089
color tn_color, structure and chain B and resi 1098
color tn_color, structure and chain B and resi 1100
color tn_color, structure and chain B and resi 1101
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
# TP (Sig. CSP in Union Site): 2
# FP (Sig. CSP -- Allosteric): 35
# TN (low CSP -- Allosteric): 47
# FN (low CSP in Union Site): 0
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2RSE_2/2RSE_csp.pdb
