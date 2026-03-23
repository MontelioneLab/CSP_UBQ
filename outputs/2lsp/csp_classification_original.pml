reinitialize
load ./outputs/2lsp/2lsp_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain B and resi 380
color tp_color, structure and chain B and resi 387
color tp_color, structure and chain B and resi 428
color tp_color, structure and chain B and resi 429
color tp_color, structure and chain B and resi 431
color tp_color, structure and chain B and resi 433
color tp_color, structure and chain B and resi 437
color tp_color, structure and chain B and resi 438
# FP (Sig. CSP -- Allosteric): 37 residues
color fp_color, structure and chain B and resi 353
color fp_color, structure and chain B and resi 354
color fp_color, structure and chain B and resi 355
color fp_color, structure and chain B and resi 356
color fp_color, structure and chain B and resi 357
color fp_color, structure and chain B and resi 358
color fp_color, structure and chain B and resi 360
color fp_color, structure and chain B and resi 361
color fp_color, structure and chain B and resi 362
color fp_color, structure and chain B and resi 363
color fp_color, structure and chain B and resi 364
color fp_color, structure and chain B and resi 372
color fp_color, structure and chain B and resi 381
color fp_color, structure and chain B and resi 388
color fp_color, structure and chain B and resi 390
color fp_color, structure and chain B and resi 391
color fp_color, structure and chain B and resi 395
color fp_color, structure and chain B and resi 396
color fp_color, structure and chain B and resi 402
color fp_color, structure and chain B and resi 408
color fp_color, structure and chain B and resi 413
color fp_color, structure and chain B and resi 414
color fp_color, structure and chain B and resi 415
color fp_color, structure and chain B and resi 422
color fp_color, structure and chain B and resi 427
color fp_color, structure and chain B and resi 430
color fp_color, structure and chain B and resi 446
color fp_color, structure and chain B and resi 447
color fp_color, structure and chain B and resi 449
color fp_color, structure and chain B and resi 450
color fp_color, structure and chain B and resi 453
color fp_color, structure and chain B and resi 454
color fp_color, structure and chain B and resi 455
color fp_color, structure and chain B and resi 456
color fp_color, structure and chain B and resi 457
color fp_color, structure and chain B and resi 459
color fp_color, structure and chain B and resi 460
# TN (low CSP -- Allosteric): 45 residues
color tn_color, structure and chain B and resi 359
color tn_color, structure and chain B and resi 365
color tn_color, structure and chain B and resi 366
color tn_color, structure and chain B and resi 368
color tn_color, structure and chain B and resi 369
color tn_color, structure and chain B and resi 370
color tn_color, structure and chain B and resi 371
color tn_color, structure and chain B and resi 373
color tn_color, structure and chain B and resi 377
color tn_color, structure and chain B and resi 378
color tn_color, structure and chain B and resi 382
color tn_color, structure and chain B and resi 383
color tn_color, structure and chain B and resi 384
color tn_color, structure and chain B and resi 394
color tn_color, structure and chain B and resi 398
color tn_color, structure and chain B and resi 399
color tn_color, structure and chain B and resi 400
color tn_color, structure and chain B and resi 401
color tn_color, structure and chain B and resi 403
color tn_color, structure and chain B and resi 404
color tn_color, structure and chain B and resi 405
color tn_color, structure and chain B and resi 406
color tn_color, structure and chain B and resi 407
color tn_color, structure and chain B and resi 409
color tn_color, structure and chain B and resi 410
color tn_color, structure and chain B and resi 411
color tn_color, structure and chain B and resi 412
color tn_color, structure and chain B and resi 416
color tn_color, structure and chain B and resi 417
color tn_color, structure and chain B and resi 418
color tn_color, structure and chain B and resi 419
color tn_color, structure and chain B and resi 420
color tn_color, structure and chain B and resi 421
color tn_color, structure and chain B and resi 423
color tn_color, structure and chain B and resi 424
color tn_color, structure and chain B and resi 425
color tn_color, structure and chain B and resi 426
color tn_color, structure and chain B and resi 436
color tn_color, structure and chain B and resi 440
color tn_color, structure and chain B and resi 441
color tn_color, structure and chain B and resi 442
color tn_color, structure and chain B and resi 443
color tn_color, structure and chain B and resi 444
color tn_color, structure and chain B and resi 445
color tn_color, structure and chain B and resi 451
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain B and resi 374
color fn_color, structure and chain B and resi 376
color fn_color, structure and chain B and resi 392
color fn_color, structure and chain B and resi 393
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 37
# TN (low CSP -- Allosteric): 45
# FN (low CSP in Union Site): 4
# Residues without CSP data: 34
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2lsp/2lsp_csp.pdb
