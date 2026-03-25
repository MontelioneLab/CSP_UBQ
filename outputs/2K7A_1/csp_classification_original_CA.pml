reinitialize
load ./outputs/2K7A_1/2K7A_csp.pdb, structure
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
color tp_color, structure and chain B and resi 280
color tp_color, structure and chain B and resi 281
color tp_color, structure and chain B and resi 282
color tp_color, structure and chain B and resi 283
color tp_color, structure and chain B and resi 290
color tp_color, structure and chain B and resi 328
color tp_color, structure and chain B and resi 330
color tp_color, structure and chain B and resi 331
color tp_color, structure and chain B and resi 332
# FP (Sig. CSP -- Allosteric): 31 residues
color fp_color, structure and chain B and resi 242
color fp_color, structure and chain B and resi 243
color fp_color, structure and chain B and resi 245
color fp_color, structure and chain B and resi 257
color fp_color, structure and chain B and resi 258
color fp_color, structure and chain B and resi 261
color fp_color, structure and chain B and resi 265
color fp_color, structure and chain B and resi 268
color fp_color, structure and chain B and resi 269
color fp_color, structure and chain B and resi 277
color fp_color, structure and chain B and resi 278
color fp_color, structure and chain B and resi 284
color fp_color, structure and chain B and resi 285
color fp_color, structure and chain B and resi 286
color fp_color, structure and chain B and resi 288
color fp_color, structure and chain B and resi 289
color fp_color, structure and chain B and resi 293
color fp_color, structure and chain B and resi 299
color fp_color, structure and chain B and resi 305
color fp_color, structure and chain B and resi 310
color fp_color, structure and chain B and resi 311
color fp_color, structure and chain B and resi 322
color fp_color, structure and chain B and resi 323
color fp_color, structure and chain B and resi 324
color fp_color, structure and chain B and resi 325
color fp_color, structure and chain B and resi 327
color fp_color, structure and chain B and resi 333
color fp_color, structure and chain B and resi 334
color fp_color, structure and chain B and resi 337
color fp_color, structure and chain B and resi 338
color fp_color, structure and chain B and resi 339
# TN (low CSP -- Allosteric): 55 residues
color tn_color, structure and chain B and resi 234
color tn_color, structure and chain B and resi 235
color tn_color, structure and chain B and resi 236
color tn_color, structure and chain B and resi 237
color tn_color, structure and chain B and resi 238
color tn_color, structure and chain B and resi 239
color tn_color, structure and chain B and resi 240
color tn_color, structure and chain B and resi 241
color tn_color, structure and chain B and resi 244
color tn_color, structure and chain B and resi 246
color tn_color, structure and chain B and resi 247
color tn_color, structure and chain B and resi 248
color tn_color, structure and chain B and resi 249
color tn_color, structure and chain B and resi 250
color tn_color, structure and chain B and resi 251
color tn_color, structure and chain B and resi 252
color tn_color, structure and chain B and resi 253
color tn_color, structure and chain B and resi 254
color tn_color, structure and chain B and resi 255
color tn_color, structure and chain B and resi 256
color tn_color, structure and chain B and resi 259
color tn_color, structure and chain B and resi 260
color tn_color, structure and chain B and resi 262
color tn_color, structure and chain B and resi 263
color tn_color, structure and chain B and resi 264
color tn_color, structure and chain B and resi 266
color tn_color, structure and chain B and resi 267
color tn_color, structure and chain B and resi 271
color tn_color, structure and chain B and resi 272
color tn_color, structure and chain B and resi 273
color tn_color, structure and chain B and resi 274
color tn_color, structure and chain B and resi 275
color tn_color, structure and chain B and resi 276
color tn_color, structure and chain B and resi 279
color tn_color, structure and chain B and resi 291
color tn_color, structure and chain B and resi 294
color tn_color, structure and chain B and resi 296
color tn_color, structure and chain B and resi 297
color tn_color, structure and chain B and resi 298
color tn_color, structure and chain B and resi 300
color tn_color, structure and chain B and resi 302
color tn_color, structure and chain B and resi 303
color tn_color, structure and chain B and resi 304
color tn_color, structure and chain B and resi 306
color tn_color, structure and chain B and resi 308
color tn_color, structure and chain B and resi 312
color tn_color, structure and chain B and resi 313
color tn_color, structure and chain B and resi 314
color tn_color, structure and chain B and resi 315
color tn_color, structure and chain B and resi 317
color tn_color, structure and chain B and resi 318
color tn_color, structure and chain B and resi 319
color tn_color, structure and chain B and resi 320
color tn_color, structure and chain B and resi 321
color tn_color, structure and chain B and resi 326
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain B and resi 292
color fn_color, structure and chain B and resi 295
color fn_color, structure and chain B and resi 307
color fn_color, structure and chain B and resi 309
color fn_color, structure and chain B and resi 329
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
# FP (Sig. CSP -- Allosteric): 31
# TN (low CSP -- Allosteric): 55
# FN (low CSP in Union Site): 5
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2K7A_1/2K7A_csp.pdb
