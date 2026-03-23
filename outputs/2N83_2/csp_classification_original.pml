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
# TP (Sig. CSP in Binding Site): 16 residues
color tp_color, structure and chain B and resi 437
color tp_color, structure and chain B and resi 469
color tp_color, structure and chain B and resi 470
color tp_color, structure and chain B and resi 472
color tp_color, structure and chain B and resi 497
color tp_color, structure and chain B and resi 498
color tp_color, structure and chain B and resi 500
color tp_color, structure and chain B and resi 503
color tp_color, structure and chain B and resi 524
color tp_color, structure and chain B and resi 526
color tp_color, structure and chain B and resi 528
color tp_color, structure and chain B and resi 529
color tp_color, structure and chain B and resi 531
color tp_color, structure and chain B and resi 532
color tp_color, structure and chain B and resi 533
color tp_color, structure and chain B and resi 534
# FP (Sig. CSP -- Allosteric): 61 residues
color fp_color, structure and chain B and resi 435
color fp_color, structure and chain B and resi 436
color fp_color, structure and chain B and resi 438
color fp_color, structure and chain B and resi 439
color fp_color, structure and chain B and resi 440
color fp_color, structure and chain B and resi 441
color fp_color, structure and chain B and resi 444
color fp_color, structure and chain B and resi 445
color fp_color, structure and chain B and resi 446
color fp_color, structure and chain B and resi 447
color fp_color, structure and chain B and resi 448
color fp_color, structure and chain B and resi 449
color fp_color, structure and chain B and resi 450
color fp_color, structure and chain B and resi 451
color fp_color, structure and chain B and resi 453
color fp_color, structure and chain B and resi 454
color fp_color, structure and chain B and resi 455
color fp_color, structure and chain B and resi 457
color fp_color, structure and chain B and resi 458
color fp_color, structure and chain B and resi 459
color fp_color, structure and chain B and resi 461
color fp_color, structure and chain B and resi 462
color fp_color, structure and chain B and resi 471
color fp_color, structure and chain B and resi 473
color fp_color, structure and chain B and resi 474
color fp_color, structure and chain B and resi 475
color fp_color, structure and chain B and resi 476
color fp_color, structure and chain B and resi 478
color fp_color, structure and chain B and resi 479
color fp_color, structure and chain B and resi 480
color fp_color, structure and chain B and resi 482
color fp_color, structure and chain B and resi 483
color fp_color, structure and chain B and resi 484
color fp_color, structure and chain B and resi 485
color fp_color, structure and chain B and resi 486
color fp_color, structure and chain B and resi 487
color fp_color, structure and chain B and resi 488
color fp_color, structure and chain B and resi 489
color fp_color, structure and chain B and resi 490
color fp_color, structure and chain B and resi 491
color fp_color, structure and chain B and resi 492
color fp_color, structure and chain B and resi 493
color fp_color, structure and chain B and resi 494
color fp_color, structure and chain B and resi 505
color fp_color, structure and chain B and resi 506
color fp_color, structure and chain B and resi 507
color fp_color, structure and chain B and resi 509
color fp_color, structure and chain B and resi 510
color fp_color, structure and chain B and resi 512
color fp_color, structure and chain B and resi 514
color fp_color, structure and chain B and resi 516
color fp_color, structure and chain B and resi 517
color fp_color, structure and chain B and resi 518
color fp_color, structure and chain B and resi 522
color fp_color, structure and chain B and resi 525
color fp_color, structure and chain B and resi 527
color fp_color, structure and chain B and resi 535
color fp_color, structure and chain B and resi 536
color fp_color, structure and chain B and resi 537
color fp_color, structure and chain B and resi 538
color fp_color, structure and chain B and resi 539
# TN (low CSP -- Allosteric): 18 residues
color tn_color, structure and chain B and resi 442
color tn_color, structure and chain B and resi 443
color tn_color, structure and chain B and resi 452
color tn_color, structure and chain B and resi 456
color tn_color, structure and chain B and resi 460
color tn_color, structure and chain B and resi 463
color tn_color, structure and chain B and resi 464
color tn_color, structure and chain B and resi 465
color tn_color, structure and chain B and resi 477
color tn_color, structure and chain B and resi 495
color tn_color, structure and chain B and resi 501
color tn_color, structure and chain B and resi 502
color tn_color, structure and chain B and resi 508
color tn_color, structure and chain B and resi 511
color tn_color, structure and chain B and resi 513
color tn_color, structure and chain B and resi 515
color tn_color, structure and chain B and resi 520
color tn_color, structure and chain B and resi 523
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain B and resi 466
color fn_color, structure and chain B and resi 467
color fn_color, structure and chain B and resi 468
color fn_color, structure and chain B and resi 496
color fn_color, structure and chain B and resi 499
color fn_color, structure and chain B and resi 504
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
# FP (Sig. CSP -- Allosteric): 61
# TN (low CSP -- Allosteric): 18
# FN (low CSP in Union Site): 6
# Residues without CSP data: 4
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2N83_2/2N83_csp.pdb
