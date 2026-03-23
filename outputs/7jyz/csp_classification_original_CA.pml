reinitialize
load ./outputs/7jyz/7jyz_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 14 residues
color tp_color, structure and chain B and resi 230
color tp_color, structure and chain B and resi 234
color tp_color, structure and chain B and resi 237
color tp_color, structure and chain B and resi 247
color tp_color, structure and chain B and resi 248
color tp_color, structure and chain B and resi 251
color tp_color, structure and chain B and resi 268
color tp_color, structure and chain B and resi 269
color tp_color, structure and chain B and resi 270
color tp_color, structure and chain B and resi 271
color tp_color, structure and chain B and resi 272
color tp_color, structure and chain B and resi 273
color tp_color, structure and chain B and resi 274
color tp_color, structure and chain B and resi 275
# FP (Sig. CSP -- Allosteric): 23 residues
color fp_color, structure and chain B and resi 210
color fp_color, structure and chain B and resi 224
color fp_color, structure and chain B and resi 229
color fp_color, structure and chain B and resi 231
color fp_color, structure and chain B and resi 236
color fp_color, structure and chain B and resi 240
color fp_color, structure and chain B and resi 242
color fp_color, structure and chain B and resi 245
color fp_color, structure and chain B and resi 249
color fp_color, structure and chain B and resi 252
color fp_color, structure and chain B and resi 253
color fp_color, structure and chain B and resi 254
color fp_color, structure and chain B and resi 255
color fp_color, structure and chain B and resi 258
color fp_color, structure and chain B and resi 259
color fp_color, structure and chain B and resi 261
color fp_color, structure and chain B and resi 262
color fp_color, structure and chain B and resi 263
color fp_color, structure and chain B and resi 264
color fp_color, structure and chain B and resi 265
color fp_color, structure and chain B and resi 266
color fp_color, structure and chain B and resi 283
color fp_color, structure and chain B and resi 286
# TN (low CSP -- Allosteric): 42 residues
color tn_color, structure and chain B and resi 211
color tn_color, structure and chain B and resi 212
color tn_color, structure and chain B and resi 213
color tn_color, structure and chain B and resi 214
color tn_color, structure and chain B and resi 215
color tn_color, structure and chain B and resi 216
color tn_color, structure and chain B and resi 217
color tn_color, structure and chain B and resi 218
color tn_color, structure and chain B and resi 219
color tn_color, structure and chain B and resi 222
color tn_color, structure and chain B and resi 223
color tn_color, structure and chain B and resi 225
color tn_color, structure and chain B and resi 226
color tn_color, structure and chain B and resi 228
color tn_color, structure and chain B and resi 232
color tn_color, structure and chain B and resi 235
color tn_color, structure and chain B and resi 238
color tn_color, structure and chain B and resi 239
color tn_color, structure and chain B and resi 241
color tn_color, structure and chain B and resi 243
color tn_color, structure and chain B and resi 246
color tn_color, structure and chain B and resi 250
color tn_color, structure and chain B and resi 256
color tn_color, structure and chain B and resi 257
color tn_color, structure and chain B and resi 276
color tn_color, structure and chain B and resi 277
color tn_color, structure and chain B and resi 278
color tn_color, structure and chain B and resi 280
color tn_color, structure and chain B and resi 281
color tn_color, structure and chain B and resi 282
color tn_color, structure and chain B and resi 284
color tn_color, structure and chain B and resi 285
color tn_color, structure and chain B and resi 287
color tn_color, structure and chain B and resi 288
color tn_color, structure and chain B and resi 289
color tn_color, structure and chain B and resi 290
color tn_color, structure and chain B and resi 291
color tn_color, structure and chain B and resi 292
color tn_color, structure and chain B and resi 293
color tn_color, structure and chain B and resi 294
color tn_color, structure and chain B and resi 295
color tn_color, structure and chain B and resi 296
# FN (low CSP in Binding Site): 1 residues
color fn_color, structure and chain B and resi 233
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
# TP (Sig. CSP in Union Site): 14
# FP (Sig. CSP -- Allosteric): 23
# TN (low CSP -- Allosteric): 42
# FN (low CSP in Union Site): 1
# Residues without CSP data: 16
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/7jyz/7jyz_csp.pdb
