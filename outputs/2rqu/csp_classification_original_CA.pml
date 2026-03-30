reinitialize
load ./outputs/2rqu/2rqu_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 11 residues
color tp_color, structure and chain A and resi 1076
color tp_color, structure and chain A and resi 1077
color tp_color, structure and chain A and resi 1079
color tp_color, structure and chain A and resi 1082
color tp_color, structure and chain A and resi 1084
color tp_color, structure and chain A and resi 1100
color tp_color, structure and chain A and resi 1101
color tp_color, structure and chain A and resi 1103
color tp_color, structure and chain A and resi 1118
color tp_color, structure and chain A and resi 1122
color tp_color, structure and chain A and resi 1123
# FP (Sig. CSP -- Allosteric): 10 residues
color fp_color, structure and chain A and resi 1080
color fp_color, structure and chain A and resi 1083
color fp_color, structure and chain A and resi 1086
color fp_color, structure and chain A and resi 1087
color fp_color, structure and chain A and resi 1088
color fp_color, structure and chain A and resi 1092
color fp_color, structure and chain A and resi 1099
color fp_color, structure and chain A and resi 1102
color fp_color, structure and chain A and resi 1119
color fp_color, structure and chain A and resi 1129
# TN (low CSP -- Allosteric): 29 residues
color tn_color, structure and chain A and resi 1071
color tn_color, structure and chain A and resi 1072
color tn_color, structure and chain A and resi 1073
color tn_color, structure and chain A and resi 1074
color tn_color, structure and chain A and resi 1078
color tn_color, structure and chain A and resi 1089
color tn_color, structure and chain A and resi 1091
color tn_color, structure and chain A and resi 1093
color tn_color, structure and chain A and resi 1094
color tn_color, structure and chain A and resi 1095
color tn_color, structure and chain A and resi 1096
color tn_color, structure and chain A and resi 1097
color tn_color, structure and chain A and resi 1098
color tn_color, structure and chain A and resi 1105
color tn_color, structure and chain A and resi 1106
color tn_color, structure and chain A and resi 1107
color tn_color, structure and chain A and resi 1108
color tn_color, structure and chain A and resi 1109
color tn_color, structure and chain A and resi 1110
color tn_color, structure and chain A and resi 1114
color tn_color, structure and chain A and resi 1115
color tn_color, structure and chain A and resi 1116
color tn_color, structure and chain A and resi 1117
color tn_color, structure and chain A and resi 1121
color tn_color, structure and chain A and resi 1124
color tn_color, structure and chain A and resi 1125
color tn_color, structure and chain A and resi 1126
color tn_color, structure and chain A and resi 1127
color tn_color, structure and chain A and resi 1128
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain A and resi 1075
color fn_color, structure and chain A and resi 1081
color fn_color, structure and chain A and resi 1085
color fn_color, structure and chain A and resi 1090
color fn_color, structure and chain A and resi 1104
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
# TP (Sig. CSP in Union Site): 11
# FP (Sig. CSP -- Allosteric): 10
# TN (low CSP -- Allosteric): 29
# FN (low CSP in Union Site): 5
# Residues without CSP data: 6
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2rqu/2rqu_csp.pdb
