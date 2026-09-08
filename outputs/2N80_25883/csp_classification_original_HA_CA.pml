reinitialize
load ./outputs/2N80_25883/2N80_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain A and resi 343
color tp_color, structure and chain A and resi 346
color tp_color, structure and chain A and resi 349
color tp_color, structure and chain A and resi 350
color tp_color, structure and chain A and resi 352
color tp_color, structure and chain A and resi 353
color tp_color, structure and chain A and resi 354
color tp_color, structure and chain A and resi 355
# FP (Sig. CSP -- Allosteric): 36 residues
color fp_color, structure and chain A and resi 334
color fp_color, structure and chain A and resi 336
color fp_color, structure and chain A and resi 337
color fp_color, structure and chain A and resi 344
color fp_color, structure and chain A and resi 345
color fp_color, structure and chain A and resi 347
color fp_color, structure and chain A and resi 348
color fp_color, structure and chain A and resi 351
color fp_color, structure and chain A and resi 358
color fp_color, structure and chain A and resi 359
color fp_color, structure and chain A and resi 360
color fp_color, structure and chain A and resi 361
color fp_color, structure and chain A and resi 363
color fp_color, structure and chain A and resi 372
color fp_color, structure and chain A and resi 373
color fp_color, structure and chain A and resi 374
color fp_color, structure and chain A and resi 376
color fp_color, structure and chain A and resi 377
color fp_color, structure and chain A and resi 379
color fp_color, structure and chain A and resi 383
color fp_color, structure and chain A and resi 386
color fp_color, structure and chain A and resi 388
color fp_color, structure and chain A and resi 391
color fp_color, structure and chain A and resi 394
color fp_color, structure and chain A and resi 395
color fp_color, structure and chain A and resi 396
color fp_color, structure and chain A and resi 397
color fp_color, structure and chain A and resi 398
color fp_color, structure and chain A and resi 399
color fp_color, structure and chain A and resi 402
color fp_color, structure and chain A and resi 408
color fp_color, structure and chain A and resi 414
color fp_color, structure and chain A and resi 417
color fp_color, structure and chain A and resi 418
color fp_color, structure and chain A and resi 421
color fp_color, structure and chain A and resi 422
# TN (low CSP -- Allosteric): 31 residues
color tn_color, structure and chain A and resi 335
color tn_color, structure and chain A and resi 338
color tn_color, structure and chain A and resi 357
color tn_color, structure and chain A and resi 362
color tn_color, structure and chain A and resi 365
color tn_color, structure and chain A and resi 366
color tn_color, structure and chain A and resi 367
color tn_color, structure and chain A and resi 368
color tn_color, structure and chain A and resi 371
color tn_color, structure and chain A and resi 375
color tn_color, structure and chain A and resi 378
color tn_color, structure and chain A and resi 380
color tn_color, structure and chain A and resi 384
color tn_color, structure and chain A and resi 385
color tn_color, structure and chain A and resi 387
color tn_color, structure and chain A and resi 389
color tn_color, structure and chain A and resi 390
color tn_color, structure and chain A and resi 392
color tn_color, structure and chain A and resi 393
color tn_color, structure and chain A and resi 400
color tn_color, structure and chain A and resi 401
color tn_color, structure and chain A and resi 403
color tn_color, structure and chain A and resi 404
color tn_color, structure and chain A and resi 405
color tn_color, structure and chain A and resi 406
color tn_color, structure and chain A and resi 407
color tn_color, structure and chain A and resi 411
color tn_color, structure and chain A and resi 415
color tn_color, structure and chain A and resi 423
color tn_color, structure and chain A and resi 424
color tn_color, structure and chain A and resi 427
# FN (low CSP in Binding Site): 8 residues
color fn_color, structure and chain A and resi 342
color fn_color, structure and chain A and resi 409
color fn_color, structure and chain A and resi 410
color fn_color, structure and chain A and resi 412
color fn_color, structure and chain A and resi 413
color fn_color, structure and chain A and resi 416
color fn_color, structure and chain A and resi 419
color fn_color, structure and chain A and resi 420
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 36
# TN (low CSP -- Allosteric): 31
# FN (low CSP in Union Site): 8
# Residues without CSP data: 11
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/GTM_Research/01_Active_Projects/new_CSP_UBQ/CSP_UBQ/outputs/2N80_25883/2N80_csp.pdb
