reinitialize
load ./outputs/2lgg/2lgg_csp.pdb, structure
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
color tp_color, structure and chain A and resi 341
color tp_color, structure and chain A and resi 342
color tp_color, structure and chain A and resi 343
color tp_color, structure and chain A and resi 344
color tp_color, structure and chain A and resi 345
color tp_color, structure and chain A and resi 346
color tp_color, structure and chain A and resi 347
color tp_color, structure and chain A and resi 349
color tp_color, structure and chain A and resi 354
color tp_color, structure and chain A and resi 369
color tp_color, structure and chain A and resi 370
# FP (Sig. CSP -- Allosteric): 11 residues
color fp_color, structure and chain A and resi 332
color fp_color, structure and chain A and resi 338
color fp_color, structure and chain A and resi 339
color fp_color, structure and chain A and resi 348
color fp_color, structure and chain A and resi 351
color fp_color, structure and chain A and resi 352
color fp_color, structure and chain A and resi 353
color fp_color, structure and chain A and resi 355
color fp_color, structure and chain A and resi 363
color fp_color, structure and chain A and resi 368
color fp_color, structure and chain A and resi 371
# TN (low CSP -- Allosteric): 22 residues
color tn_color, structure and chain A and resi 312
color tn_color, structure and chain A and resi 314
color tn_color, structure and chain A and resi 316
color tn_color, structure and chain A and resi 317
color tn_color, structure and chain A and resi 319
color tn_color, structure and chain A and resi 324
color tn_color, structure and chain A and resi 327
color tn_color, structure and chain A and resi 328
color tn_color, structure and chain A and resi 333
color tn_color, structure and chain A and resi 335
color tn_color, structure and chain A and resi 336
color tn_color, structure and chain A and resi 356
color tn_color, structure and chain A and resi 358
color tn_color, structure and chain A and resi 359
color tn_color, structure and chain A and resi 362
color tn_color, structure and chain A and resi 364
color tn_color, structure and chain A and resi 365
color tn_color, structure and chain A and resi 372
color tn_color, structure and chain A and resi 375
color tn_color, structure and chain A and resi 377
color tn_color, structure and chain A and resi 378
color tn_color, structure and chain A and resi 379
# FN (low CSP in Binding Site): 17 residues
color fn_color, structure and chain A and resi 315
color fn_color, structure and chain A and resi 318
color fn_color, structure and chain A and resi 320
color fn_color, structure and chain A and resi 321
color fn_color, structure and chain A and resi 322
color fn_color, structure and chain A and resi 323
color fn_color, structure and chain A and resi 325
color fn_color, structure and chain A and resi 326
color fn_color, structure and chain A and resi 329
color fn_color, structure and chain A and resi 330
color fn_color, structure and chain A and resi 331
color fn_color, structure and chain A and resi 334
color fn_color, structure and chain A and resi 337
color fn_color, structure and chain A and resi 357
color fn_color, structure and chain A and resi 367
color fn_color, structure and chain A and resi 373
color fn_color, structure and chain A and resi 376
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
# FP (Sig. CSP -- Allosteric): 11
# TN (low CSP -- Allosteric): 22
# FN (low CSP in Union Site): 17
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2lgg/2lgg_csp.pdb
