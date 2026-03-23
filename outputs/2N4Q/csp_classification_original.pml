reinitialize
load ./outputs/2N4Q/2N4Q_csp.pdb, structure
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
color tp_color, structure and chain B and resi 508
color tp_color, structure and chain B and resi 524
color tp_color, structure and chain B and resi 527
color tp_color, structure and chain B and resi 542
color tp_color, structure and chain B and resi 543
color tp_color, structure and chain B and resi 545
color tp_color, structure and chain B and resi 546
color tp_color, structure and chain B and resi 547
# FP (Sig. CSP -- Allosteric): 22 residues
color fp_color, structure and chain B and resi 501
color fp_color, structure and chain B and resi 502
color fp_color, structure and chain B and resi 512
color fp_color, structure and chain B and resi 514
color fp_color, structure and chain B and resi 515
color fp_color, structure and chain B and resi 519
color fp_color, structure and chain B and resi 520
color fp_color, structure and chain B and resi 522
color fp_color, structure and chain B and resi 525
color fp_color, structure and chain B and resi 526
color fp_color, structure and chain B and resi 528
color fp_color, structure and chain B and resi 530
color fp_color, structure and chain B and resi 531
color fp_color, structure and chain B and resi 532
color fp_color, structure and chain B and resi 536
color fp_color, structure and chain B and resi 537
color fp_color, structure and chain B and resi 538
color fp_color, structure and chain B and resi 539
color fp_color, structure and chain B and resi 540
color fp_color, structure and chain B and resi 541
color fp_color, structure and chain B and resi 548
color fp_color, structure and chain B and resi 551
# TN (low CSP -- Allosteric): 34 residues
color tn_color, structure and chain B and resi 503
color tn_color, structure and chain B and resi 505
color tn_color, structure and chain B and resi 506
color tn_color, structure and chain B and resi 507
color tn_color, structure and chain B and resi 509
color tn_color, structure and chain B and resi 510
color tn_color, structure and chain B and resi 513
color tn_color, structure and chain B and resi 516
color tn_color, structure and chain B and resi 517
color tn_color, structure and chain B and resi 518
color tn_color, structure and chain B and resi 523
color tn_color, structure and chain B and resi 529
color tn_color, structure and chain B and resi 533
color tn_color, structure and chain B and resi 534
color tn_color, structure and chain B and resi 535
color tn_color, structure and chain B and resi 549
color tn_color, structure and chain B and resi 550
color tn_color, structure and chain B and resi 552
color tn_color, structure and chain B and resi 553
color tn_color, structure and chain B and resi 554
color tn_color, structure and chain B and resi 555
color tn_color, structure and chain B and resi 556
color tn_color, structure and chain B and resi 557
color tn_color, structure and chain B and resi 558
color tn_color, structure and chain B and resi 559
color tn_color, structure and chain B and resi 560
color tn_color, structure and chain B and resi 561
color tn_color, structure and chain B and resi 562
color tn_color, structure and chain B and resi 563
color tn_color, structure and chain B and resi 564
color tn_color, structure and chain B and resi 565
color tn_color, structure and chain B and resi 566
color tn_color, structure and chain B and resi 567
color tn_color, structure and chain B and resi 568
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain B and resi 504
color fn_color, structure and chain B and resi 511
color fn_color, structure and chain B and resi 544
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
# FP (Sig. CSP -- Allosteric): 22
# TN (low CSP -- Allosteric): 34
# FN (low CSP in Union Site): 3
# Residues without CSP data: 3
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2N4Q/2N4Q_csp.pdb
