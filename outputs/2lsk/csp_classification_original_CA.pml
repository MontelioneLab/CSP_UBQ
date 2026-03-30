reinitialize
load ./outputs/2lsk/2lsk_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 9 residues
color tp_color, structure and chain A and resi 1159
color tp_color, structure and chain A and resi 1160
color tp_color, structure and chain A and resi 1161
color tp_color, structure and chain A and resi 1175
color tp_color, structure and chain A and resi 1179
color tp_color, structure and chain A and resi 1183
color tp_color, structure and chain A and resi 1185
color tp_color, structure and chain A and resi 1186
color tp_color, structure and chain A and resi 1189
# FP (Sig. CSP -- Allosteric): 31 residues
color fp_color, structure and chain A and resi 1162
color fp_color, structure and chain A and resi 1163
color fp_color, structure and chain A and resi 1164
color fp_color, structure and chain A and resi 1167
color fp_color, structure and chain A and resi 1168
color fp_color, structure and chain A and resi 1169
color fp_color, structure and chain A and resi 1170
color fp_color, structure and chain A and resi 1171
color fp_color, structure and chain A and resi 1173
color fp_color, structure and chain A and resi 1177
color fp_color, structure and chain A and resi 1178
color fp_color, structure and chain A and resi 1184
color fp_color, structure and chain A and resi 1187
color fp_color, structure and chain A and resi 1188
color fp_color, structure and chain A and resi 1191
color fp_color, structure and chain A and resi 1192
color fp_color, structure and chain A and resi 1196
color fp_color, structure and chain A and resi 1203
color fp_color, structure and chain A and resi 1208
color fp_color, structure and chain A and resi 1209
color fp_color, structure and chain A and resi 1210
color fp_color, structure and chain A and resi 1212
color fp_color, structure and chain A and resi 1213
color fp_color, structure and chain A and resi 1216
color fp_color, structure and chain A and resi 1219
color fp_color, structure and chain A and resi 1228
color fp_color, structure and chain A and resi 1229
color fp_color, structure and chain A and resi 1230
color fp_color, structure and chain A and resi 1231
color fp_color, structure and chain A and resi 1232
color fp_color, structure and chain A and resi 1245
# TN (low CSP -- Allosteric): 40 residues
color tn_color, structure and chain A and resi 1166
color tn_color, structure and chain A and resi 1193
color tn_color, structure and chain A and resi 1194
color tn_color, structure and chain A and resi 1195
color tn_color, structure and chain A and resi 1198
color tn_color, structure and chain A and resi 1199
color tn_color, structure and chain A and resi 1200
color tn_color, structure and chain A and resi 1201
color tn_color, structure and chain A and resi 1202
color tn_color, structure and chain A and resi 1204
color tn_color, structure and chain A and resi 1205
color tn_color, structure and chain A and resi 1206
color tn_color, structure and chain A and resi 1207
color tn_color, structure and chain A and resi 1211
color tn_color, structure and chain A and resi 1214
color tn_color, structure and chain A and resi 1215
color tn_color, structure and chain A and resi 1217
color tn_color, structure and chain A and resi 1218
color tn_color, structure and chain A and resi 1221
color tn_color, structure and chain A and resi 1222
color tn_color, structure and chain A and resi 1226
color tn_color, structure and chain A and resi 1227
color tn_color, structure and chain A and resi 1233
color tn_color, structure and chain A and resi 1234
color tn_color, structure and chain A and resi 1235
color tn_color, structure and chain A and resi 1236
color tn_color, structure and chain A and resi 1237
color tn_color, structure and chain A and resi 1238
color tn_color, structure and chain A and resi 1239
color tn_color, structure and chain A and resi 1240
color tn_color, structure and chain A and resi 1241
color tn_color, structure and chain A and resi 1242
color tn_color, structure and chain A and resi 1243
color tn_color, structure and chain A and resi 1244
color tn_color, structure and chain A and resi 1246
color tn_color, structure and chain A and resi 1247
color tn_color, structure and chain A and resi 1248
color tn_color, structure and chain A and resi 1249
color tn_color, structure and chain A and resi 1250
color tn_color, structure and chain A and resi 1251
# FN (low CSP in Binding Site): 1 residues
color fn_color, structure and chain A and resi 1190
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
# TP (Sig. CSP in Union Site): 9
# FP (Sig. CSP -- Allosteric): 31
# TN (low CSP -- Allosteric): 40
# FN (low CSP in Union Site): 1
# Residues without CSP data: 14
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2lsk/2lsk_csp.pdb
