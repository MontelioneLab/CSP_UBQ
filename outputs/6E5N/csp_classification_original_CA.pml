reinitialize
load ./outputs/6E5N/6E5N_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 9 residues
color tp_color, structure and chain B and resi 1054
color tp_color, structure and chain B and resi 1058
color tp_color, structure and chain B and resi 1088
color tp_color, structure and chain B and resi 1090
color tp_color, structure and chain B and resi 1091
color tp_color, structure and chain B and resi 1117
color tp_color, structure and chain B and resi 1120
color tp_color, structure and chain B and resi 1121
color tp_color, structure and chain B and resi 1124
# FP (Sig. CSP -- Allosteric): 19 residues
color fp_color, structure and chain B and resi 1056
color fp_color, structure and chain B and resi 1059
color fp_color, structure and chain B and resi 1060
color fp_color, structure and chain B and resi 1061
color fp_color, structure and chain B and resi 1063
color fp_color, structure and chain B and resi 1064
color fp_color, structure and chain B and resi 1079
color fp_color, structure and chain B and resi 1080
color fp_color, structure and chain B and resi 1089
color fp_color, structure and chain B and resi 1092
color fp_color, structure and chain B and resi 1094
color fp_color, structure and chain B and resi 1095
color fp_color, structure and chain B and resi 1112
color fp_color, structure and chain B and resi 1114
color fp_color, structure and chain B and resi 1118
color fp_color, structure and chain B and resi 1122
color fp_color, structure and chain B and resi 1125
color fp_color, structure and chain B and resi 1126
color fp_color, structure and chain B and resi 1128
# TN (low CSP -- Allosteric): 47 residues
color tn_color, structure and chain B and resi 1052
color tn_color, structure and chain B and resi 1057
color tn_color, structure and chain B and resi 1065
color tn_color, structure and chain B and resi 1066
color tn_color, structure and chain B and resi 1067
color tn_color, structure and chain B and resi 1068
color tn_color, structure and chain B and resi 1069
color tn_color, structure and chain B and resi 1071
color tn_color, structure and chain B and resi 1072
color tn_color, structure and chain B and resi 1073
color tn_color, structure and chain B and resi 1074
color tn_color, structure and chain B and resi 1075
color tn_color, structure and chain B and resi 1076
color tn_color, structure and chain B and resi 1077
color tn_color, structure and chain B and resi 1078
color tn_color, structure and chain B and resi 1081
color tn_color, structure and chain B and resi 1082
color tn_color, structure and chain B and resi 1083
color tn_color, structure and chain B and resi 1084
color tn_color, structure and chain B and resi 1085
color tn_color, structure and chain B and resi 1086
color tn_color, structure and chain B and resi 1093
color tn_color, structure and chain B and resi 1096
color tn_color, structure and chain B and resi 1097
color tn_color, structure and chain B and resi 1098
color tn_color, structure and chain B and resi 1099
color tn_color, structure and chain B and resi 1100
color tn_color, structure and chain B and resi 1101
color tn_color, structure and chain B and resi 1102
color tn_color, structure and chain B and resi 1103
color tn_color, structure and chain B and resi 1104
color tn_color, structure and chain B and resi 1105
color tn_color, structure and chain B and resi 1106
color tn_color, structure and chain B and resi 1107
color tn_color, structure and chain B and resi 1108
color tn_color, structure and chain B and resi 1109
color tn_color, structure and chain B and resi 1110
color tn_color, structure and chain B and resi 1111
color tn_color, structure and chain B and resi 1113
color tn_color, structure and chain B and resi 1115
color tn_color, structure and chain B and resi 1116
color tn_color, structure and chain B and resi 1119
color tn_color, structure and chain B and resi 1123
color tn_color, structure and chain B and resi 1127
color tn_color, structure and chain B and resi 1129
color tn_color, structure and chain B and resi 1130
color tn_color, structure and chain B and resi 1131
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain B and resi 1050
color fn_color, structure and chain B and resi 1053
color fn_color, structure and chain B and resi 1062
color fn_color, structure and chain B and resi 1087
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
# TP (Sig. CSP in Union Site): 9
# FP (Sig. CSP -- Allosteric): 19
# TN (low CSP -- Allosteric): 47
# FN (low CSP in Union Site): 4
# Residues without CSP data: 3
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/6E5N/6E5N_csp.pdb
