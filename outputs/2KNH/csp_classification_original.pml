reinitialize
load ./outputs/2KNH/2KNH_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 7 residues
color tp_color, structure and chain A and resi 273
color tp_color, structure and chain A and resi 277
color tp_color, structure and chain A and resi 280
color tp_color, structure and chain A and resi 284
color tp_color, structure and chain A and resi 326
color tp_color, structure and chain A and resi 329
color tp_color, structure and chain A and resi 332
# FP (Sig. CSP -- Allosteric): 31 residues
color fp_color, structure and chain A and resi 269
color fp_color, structure and chain A and resi 271
color fp_color, structure and chain A and resi 274
color fp_color, structure and chain A and resi 275
color fp_color, structure and chain A and resi 279
color fp_color, structure and chain A and resi 281
color fp_color, structure and chain A and resi 282
color fp_color, structure and chain A and resi 283
color fp_color, structure and chain A and resi 285
color fp_color, structure and chain A and resi 301
color fp_color, structure and chain A and resi 304
color fp_color, structure and chain A and resi 305
color fp_color, structure and chain A and resi 311
color fp_color, structure and chain A and resi 321
color fp_color, structure and chain A and resi 328
color fp_color, structure and chain A and resi 330
color fp_color, structure and chain A and resi 333
color fp_color, structure and chain A and resi 342
color fp_color, structure and chain A and resi 343
color fp_color, structure and chain A and resi 345
color fp_color, structure and chain A and resi 348
color fp_color, structure and chain A and resi 350
color fp_color, structure and chain A and resi 351
color fp_color, structure and chain A and resi 353
color fp_color, structure and chain A and resi 354
color fp_color, structure and chain A and resi 357
color fp_color, structure and chain A and resi 358
color fp_color, structure and chain A and resi 360
color fp_color, structure and chain A and resi 361
color fp_color, structure and chain A and resi 362
color fp_color, structure and chain A and resi 364
# TN (low CSP -- Allosteric): 46 residues
color tn_color, structure and chain A and resi 270
color tn_color, structure and chain A and resi 272
color tn_color, structure and chain A and resi 278
color tn_color, structure and chain A and resi 286
color tn_color, structure and chain A and resi 287
color tn_color, structure and chain A and resi 288
color tn_color, structure and chain A and resi 289
color tn_color, structure and chain A and resi 291
color tn_color, structure and chain A and resi 292
color tn_color, structure and chain A and resi 293
color tn_color, structure and chain A and resi 294
color tn_color, structure and chain A and resi 295
color tn_color, structure and chain A and resi 296
color tn_color, structure and chain A and resi 297
color tn_color, structure and chain A and resi 298
color tn_color, structure and chain A and resi 299
color tn_color, structure and chain A and resi 300
color tn_color, structure and chain A and resi 302
color tn_color, structure and chain A and resi 303
color tn_color, structure and chain A and resi 306
color tn_color, structure and chain A and resi 307
color tn_color, structure and chain A and resi 308
color tn_color, structure and chain A and resi 309
color tn_color, structure and chain A and resi 310
color tn_color, structure and chain A and resi 312
color tn_color, structure and chain A and resi 313
color tn_color, structure and chain A and resi 314
color tn_color, structure and chain A and resi 315
color tn_color, structure and chain A and resi 316
color tn_color, structure and chain A and resi 317
color tn_color, structure and chain A and resi 318
color tn_color, structure and chain A and resi 320
color tn_color, structure and chain A and resi 334
color tn_color, structure and chain A and resi 335
color tn_color, structure and chain A and resi 336
color tn_color, structure and chain A and resi 337
color tn_color, structure and chain A and resi 339
color tn_color, structure and chain A and resi 340
color tn_color, structure and chain A and resi 341
color tn_color, structure and chain A and resi 344
color tn_color, structure and chain A and resi 346
color tn_color, structure and chain A and resi 347
color tn_color, structure and chain A and resi 349
color tn_color, structure and chain A and resi 352
color tn_color, structure and chain A and resi 356
color tn_color, structure and chain A and resi 363
# FN (low CSP in Binding Site): 1 residues
color fn_color, structure and chain A and resi 276
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
# TP (Sig. CSP in Union Site): 7
# FP (Sig. CSP -- Allosteric): 31
# TN (low CSP -- Allosteric): 46
# FN (low CSP in Union Site): 1
# Residues without CSP data: 18
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2KNH/2KNH_csp.pdb
