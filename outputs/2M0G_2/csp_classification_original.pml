reinitialize
load ./outputs/2M0G_2/2M0G_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 16 residues
color tp_color, structure and chain B and resi 373
color tp_color, structure and chain B and resi 394
color tp_color, structure and chain B and resi 397
color tp_color, structure and chain B and resi 400
color tp_color, structure and chain B and resi 404
color tp_color, structure and chain B and resi 405
color tp_color, structure and chain B and resi 449
color tp_color, structure and chain B and resi 451
color tp_color, structure and chain B and resi 452
color tp_color, structure and chain B and resi 454
color tp_color, structure and chain B and resi 455
color tp_color, structure and chain B and resi 456
color tp_color, structure and chain B and resi 457
color tp_color, structure and chain B and resi 458
color tp_color, structure and chain B and resi 459
color tp_color, structure and chain B and resi 460
# FP (Sig. CSP -- Allosteric): 27 residues
color fp_color, structure and chain B and resi 375
color fp_color, structure and chain B and resi 378
color fp_color, structure and chain B and resi 379
color fp_color, structure and chain B and resi 380
color fp_color, structure and chain B and resi 381
color fp_color, structure and chain B and resi 383
color fp_color, structure and chain B and resi 384
color fp_color, structure and chain B and resi 385
color fp_color, structure and chain B and resi 387
color fp_color, structure and chain B and resi 388
color fp_color, structure and chain B and resi 390
color fp_color, structure and chain B and resi 391
color fp_color, structure and chain B and resi 393
color fp_color, structure and chain B and resi 395
color fp_color, structure and chain B and resi 398
color fp_color, structure and chain B and resi 399
color fp_color, structure and chain B and resi 402
color fp_color, structure and chain B and resi 403
color fp_color, structure and chain B and resi 406
color fp_color, structure and chain B and resi 407
color fp_color, structure and chain B and resi 408
color fp_color, structure and chain B and resi 431
color fp_color, structure and chain B and resi 432
color fp_color, structure and chain B and resi 447
color fp_color, structure and chain B and resi 473
color fp_color, structure and chain B and resi 474
color fp_color, structure and chain B and resi 475
# TN (low CSP -- Allosteric): 47 residues
color tn_color, structure and chain B and resi 376
color tn_color, structure and chain B and resi 377
color tn_color, structure and chain B and resi 389
color tn_color, structure and chain B and resi 392
color tn_color, structure and chain B and resi 409
color tn_color, structure and chain B and resi 410
color tn_color, structure and chain B and resi 411
color tn_color, structure and chain B and resi 412
color tn_color, structure and chain B and resi 413
color tn_color, structure and chain B and resi 414
color tn_color, structure and chain B and resi 415
color tn_color, structure and chain B and resi 416
color tn_color, structure and chain B and resi 417
color tn_color, structure and chain B and resi 419
color tn_color, structure and chain B and resi 421
color tn_color, structure and chain B and resi 422
color tn_color, structure and chain B and resi 423
color tn_color, structure and chain B and resi 424
color tn_color, structure and chain B and resi 425
color tn_color, structure and chain B and resi 426
color tn_color, structure and chain B and resi 428
color tn_color, structure and chain B and resi 429
color tn_color, structure and chain B and resi 430
color tn_color, structure and chain B and resi 433
color tn_color, structure and chain B and resi 434
color tn_color, structure and chain B and resi 435
color tn_color, structure and chain B and resi 436
color tn_color, structure and chain B and resi 437
color tn_color, structure and chain B and resi 438
color tn_color, structure and chain B and resi 439
color tn_color, structure and chain B and resi 440
color tn_color, structure and chain B and resi 441
color tn_color, structure and chain B and resi 442
color tn_color, structure and chain B and resi 443
color tn_color, structure and chain B and resi 444
color tn_color, structure and chain B and resi 445
color tn_color, structure and chain B and resi 446
color tn_color, structure and chain B and resi 448
color tn_color, structure and chain B and resi 463
color tn_color, structure and chain B and resi 464
color tn_color, structure and chain B and resi 465
color tn_color, structure and chain B and resi 467
color tn_color, structure and chain B and resi 468
color tn_color, structure and chain B and resi 469
color tn_color, structure and chain B and resi 470
color tn_color, structure and chain B and resi 471
color tn_color, structure and chain B and resi 472
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain B and resi 382
color fn_color, structure and chain B and resi 396
color fn_color, structure and chain B and resi 401
color fn_color, structure and chain B and resi 450
color fn_color, structure and chain B and resi 462
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
# TP (Sig. CSP in Union Site): 16
# FP (Sig. CSP -- Allosteric): 27
# TN (low CSP -- Allosteric): 47
# FN (low CSP in Union Site): 5
# Residues without CSP data: 9
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2M0G_2/2M0G_csp.pdb
