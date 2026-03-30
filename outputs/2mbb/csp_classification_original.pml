reinitialize
load ./outputs/2mbb/2mbb_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 10 residues
color tp_color, structure and chain B and resi 206
color tp_color, structure and chain B and resi 207
color tp_color, structure and chain B and resi 208
color tp_color, structure and chain B and resi 209
color tp_color, structure and chain B and resi 210
color tp_color, structure and chain B and resi 244
color tp_color, structure and chain B and resi 247
color tp_color, structure and chain B and resi 270
color tp_color, structure and chain B and resi 272
color tp_color, structure and chain B and resi 273
# FP (Sig. CSP -- Allosteric): 17 residues
color fp_color, structure and chain B and resi 202
color fp_color, structure and chain B and resi 203
color fp_color, structure and chain B and resi 205
color fp_color, structure and chain B and resi 211
color fp_color, structure and chain B and resi 213
color fp_color, structure and chain B and resi 214
color fp_color, structure and chain B and resi 217
color fp_color, structure and chain B and resi 218
color fp_color, structure and chain B and resi 229
color fp_color, structure and chain B and resi 243
color fp_color, structure and chain B and resi 248
color fp_color, structure and chain B and resi 249
color fp_color, structure and chain B and resi 250
color fp_color, structure and chain B and resi 257
color fp_color, structure and chain B and resi 261
color fp_color, structure and chain B and resi 267
color fp_color, structure and chain B and resi 269
# TN (low CSP -- Allosteric): 36 residues
color tn_color, structure and chain B and resi 204
color tn_color, structure and chain B and resi 212
color tn_color, structure and chain B and resi 215
color tn_color, structure and chain B and resi 216
color tn_color, structure and chain B and resi 220
color tn_color, structure and chain B and resi 221
color tn_color, structure and chain B and resi 222
color tn_color, structure and chain B and resi 223
color tn_color, structure and chain B and resi 225
color tn_color, structure and chain B and resi 226
color tn_color, structure and chain B and resi 227
color tn_color, structure and chain B and resi 228
color tn_color, structure and chain B and resi 230
color tn_color, structure and chain B and resi 231
color tn_color, structure and chain B and resi 232
color tn_color, structure and chain B and resi 233
color tn_color, structure and chain B and resi 235
color tn_color, structure and chain B and resi 239
color tn_color, structure and chain B and resi 241
color tn_color, structure and chain B and resi 245
color tn_color, structure and chain B and resi 246
color tn_color, structure and chain B and resi 251
color tn_color, structure and chain B and resi 252
color tn_color, structure and chain B and resi 254
color tn_color, structure and chain B and resi 255
color tn_color, structure and chain B and resi 256
color tn_color, structure and chain B and resi 258
color tn_color, structure and chain B and resi 259
color tn_color, structure and chain B and resi 260
color tn_color, structure and chain B and resi 262
color tn_color, structure and chain B and resi 263
color tn_color, structure and chain B and resi 264
color tn_color, structure and chain B and resi 265
color tn_color, structure and chain B and resi 266
color tn_color, structure and chain B and resi 275
color tn_color, structure and chain B and resi 276
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain B and resi 234
color fn_color, structure and chain B and resi 236
color fn_color, structure and chain B and resi 240
color fn_color, structure and chain B and resi 268
color fn_color, structure and chain B and resi 274
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
# TP (Sig. CSP in Union Site): 10
# FP (Sig. CSP -- Allosteric): 17
# TN (low CSP -- Allosteric): 36
# FN (low CSP in Union Site): 5
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2mbb/2mbb_csp.pdb
