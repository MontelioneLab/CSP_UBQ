reinitialize
load ./outputs/2N83_2/2N83_csp.pdb, structure
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
color tp_color, structure and chain B and resi 466
color tp_color, structure and chain B and resi 467
color tp_color, structure and chain B and resi 468
color tp_color, structure and chain B and resi 469
color tp_color, structure and chain B and resi 472
color tp_color, structure and chain B and resi 496
color tp_color, structure and chain B and resi 498
color tp_color, structure and chain B and resi 499
color tp_color, structure and chain B and resi 504
# FP (Sig. CSP -- Allosteric): 25 residues
color fp_color, structure and chain B and resi 436
color fp_color, structure and chain B and resi 442
color fp_color, structure and chain B and resi 443
color fp_color, structure and chain B and resi 452
color fp_color, structure and chain B and resi 456
color fp_color, structure and chain B and resi 460
color fp_color, structure and chain B and resi 463
color fp_color, structure and chain B and resi 464
color fp_color, structure and chain B and resi 465
color fp_color, structure and chain B and resi 471
color fp_color, structure and chain B and resi 475
color fp_color, structure and chain B and resi 476
color fp_color, structure and chain B and resi 477
color fp_color, structure and chain B and resi 491
color fp_color, structure and chain B and resi 494
color fp_color, structure and chain B and resi 495
color fp_color, structure and chain B and resi 501
color fp_color, structure and chain B and resi 502
color fp_color, structure and chain B and resi 505
color fp_color, structure and chain B and resi 508
color fp_color, structure and chain B and resi 511
color fp_color, structure and chain B and resi 513
color fp_color, structure and chain B and resi 515
color fp_color, structure and chain B and resi 520
color fp_color, structure and chain B and resi 523
# TN (low CSP -- Allosteric): 54 residues
color tn_color, structure and chain B and resi 435
color tn_color, structure and chain B and resi 438
color tn_color, structure and chain B and resi 439
color tn_color, structure and chain B and resi 440
color tn_color, structure and chain B and resi 441
color tn_color, structure and chain B and resi 444
color tn_color, structure and chain B and resi 445
color tn_color, structure and chain B and resi 446
color tn_color, structure and chain B and resi 447
color tn_color, structure and chain B and resi 448
color tn_color, structure and chain B and resi 449
color tn_color, structure and chain B and resi 450
color tn_color, structure and chain B and resi 451
color tn_color, structure and chain B and resi 453
color tn_color, structure and chain B and resi 454
color tn_color, structure and chain B and resi 455
color tn_color, structure and chain B and resi 457
color tn_color, structure and chain B and resi 458
color tn_color, structure and chain B and resi 459
color tn_color, structure and chain B and resi 461
color tn_color, structure and chain B and resi 462
color tn_color, structure and chain B and resi 473
color tn_color, structure and chain B and resi 474
color tn_color, structure and chain B and resi 478
color tn_color, structure and chain B and resi 479
color tn_color, structure and chain B and resi 480
color tn_color, structure and chain B and resi 482
color tn_color, structure and chain B and resi 483
color tn_color, structure and chain B and resi 484
color tn_color, structure and chain B and resi 485
color tn_color, structure and chain B and resi 486
color tn_color, structure and chain B and resi 487
color tn_color, structure and chain B and resi 488
color tn_color, structure and chain B and resi 489
color tn_color, structure and chain B and resi 490
color tn_color, structure and chain B and resi 492
color tn_color, structure and chain B and resi 493
color tn_color, structure and chain B and resi 506
color tn_color, structure and chain B and resi 507
color tn_color, structure and chain B and resi 509
color tn_color, structure and chain B and resi 510
color tn_color, structure and chain B and resi 512
color tn_color, structure and chain B and resi 514
color tn_color, structure and chain B and resi 516
color tn_color, structure and chain B and resi 517
color tn_color, structure and chain B and resi 518
color tn_color, structure and chain B and resi 522
color tn_color, structure and chain B and resi 525
color tn_color, structure and chain B and resi 527
color tn_color, structure and chain B and resi 535
color tn_color, structure and chain B and resi 536
color tn_color, structure and chain B and resi 537
color tn_color, structure and chain B and resi 538
color tn_color, structure and chain B and resi 539
# FN (low CSP in Binding Site): 13 residues
color fn_color, structure and chain B and resi 437
color fn_color, structure and chain B and resi 470
color fn_color, structure and chain B and resi 497
color fn_color, structure and chain B and resi 500
color fn_color, structure and chain B and resi 503
color fn_color, structure and chain B and resi 524
color fn_color, structure and chain B and resi 526
color fn_color, structure and chain B and resi 528
color fn_color, structure and chain B and resi 529
color fn_color, structure and chain B and resi 531
color fn_color, structure and chain B and resi 532
color fn_color, structure and chain B and resi 533
color fn_color, structure and chain B and resi 534
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
# FP (Sig. CSP -- Allosteric): 25
# TN (low CSP -- Allosteric): 54
# FN (low CSP in Union Site): 13
# Residues without CSP data: 4
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2N83_2/2N83_csp.pdb
