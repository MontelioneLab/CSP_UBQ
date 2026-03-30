reinitialize
load ./outputs/2M3M/2M3M_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 12 residues
color tp_color, structure and chain A and resi 327
color tp_color, structure and chain A and resi 330
color tp_color, structure and chain A and resi 331
color tp_color, structure and chain A and resi 332
color tp_color, structure and chain A and resi 338
color tp_color, structure and chain A and resi 339
color tp_color, structure and chain A and resi 354
color tp_color, structure and chain A and resi 384
color tp_color, structure and chain A and resi 385
color tp_color, structure and chain A and resi 389
color tp_color, structure and chain A and resi 391
color tp_color, structure and chain A and resi 392
# FP (Sig. CSP -- Allosteric): 27 residues
color fp_color, structure and chain A and resi 319
color fp_color, structure and chain A and resi 335
color fp_color, structure and chain A and resi 336
color fp_color, structure and chain A and resi 337
color fp_color, structure and chain A and resi 341
color fp_color, structure and chain A and resi 346
color fp_color, structure and chain A and resi 348
color fp_color, structure and chain A and resi 349
color fp_color, structure and chain A and resi 353
color fp_color, structure and chain A and resi 355
color fp_color, structure and chain A and resi 356
color fp_color, structure and chain A and resi 358
color fp_color, structure and chain A and resi 360
color fp_color, structure and chain A and resi 371
color fp_color, structure and chain A and resi 372
color fp_color, structure and chain A and resi 374
color fp_color, structure and chain A and resi 377
color fp_color, structure and chain A and resi 378
color fp_color, structure and chain A and resi 379
color fp_color, structure and chain A and resi 380
color fp_color, structure and chain A and resi 381
color fp_color, structure and chain A and resi 390
color fp_color, structure and chain A and resi 395
color fp_color, structure and chain A and resi 402
color fp_color, structure and chain A and resi 403
color fp_color, structure and chain A and resi 404
color fp_color, structure and chain A and resi 408
# TN (low CSP -- Allosteric): 39 residues
color tn_color, structure and chain A and resi 320
color tn_color, structure and chain A and resi 321
color tn_color, structure and chain A and resi 322
color tn_color, structure and chain A and resi 323
color tn_color, structure and chain A and resi 324
color tn_color, structure and chain A and resi 325
color tn_color, structure and chain A and resi 340
color tn_color, structure and chain A and resi 342
color tn_color, structure and chain A and resi 344
color tn_color, structure and chain A and resi 345
color tn_color, structure and chain A and resi 347
color tn_color, structure and chain A and resi 350
color tn_color, structure and chain A and resi 357
color tn_color, structure and chain A and resi 359
color tn_color, structure and chain A and resi 361
color tn_color, structure and chain A and resi 362
color tn_color, structure and chain A and resi 363
color tn_color, structure and chain A and resi 364
color tn_color, structure and chain A and resi 365
color tn_color, structure and chain A and resi 366
color tn_color, structure and chain A and resi 367
color tn_color, structure and chain A and resi 368
color tn_color, structure and chain A and resi 369
color tn_color, structure and chain A and resi 370
color tn_color, structure and chain A and resi 373
color tn_color, structure and chain A and resi 375
color tn_color, structure and chain A and resi 382
color tn_color, structure and chain A and resi 383
color tn_color, structure and chain A and resi 386
color tn_color, structure and chain A and resi 387
color tn_color, structure and chain A and resi 393
color tn_color, structure and chain A and resi 394
color tn_color, structure and chain A and resi 396
color tn_color, structure and chain A and resi 397
color tn_color, structure and chain A and resi 398
color tn_color, structure and chain A and resi 399
color tn_color, structure and chain A and resi 400
color tn_color, structure and chain A and resi 401
color tn_color, structure and chain A and resi 406
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain A and resi 328
color fn_color, structure and chain A and resi 333
color fn_color, structure and chain A and resi 334
color fn_color, structure and chain A and resi 351
color fn_color, structure and chain A and resi 352
color fn_color, structure and chain A and resi 388
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
# TP (Sig. CSP in Union Site): 12
# FP (Sig. CSP -- Allosteric): 27
# TN (low CSP -- Allosteric): 39
# FN (low CSP in Union Site): 6
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2M3M/2M3M_csp.pdb
