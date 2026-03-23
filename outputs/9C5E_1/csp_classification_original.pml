reinitialize
load ./outputs/9C5E_1/9C5E_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 9 residues
color tp_color, structure and chain A and resi 422
color tp_color, structure and chain A and resi 424
color tp_color, structure and chain A and resi 439
color tp_color, structure and chain A and resi 441
color tp_color, structure and chain A and resi 444
color tp_color, structure and chain A and resi 445
color tp_color, structure and chain A and resi 446
color tp_color, structure and chain A and resi 464
color tp_color, structure and chain A and resi 479
# FP (Sig. CSP -- Allosteric): 17 residues
color fp_color, structure and chain A and resi 418
color fp_color, structure and chain A and resi 420
color fp_color, structure and chain A and resi 427
color fp_color, structure and chain A and resi 431
color fp_color, structure and chain A and resi 433
color fp_color, structure and chain A and resi 438
color fp_color, structure and chain A and resi 443
color fp_color, structure and chain A and resi 447
color fp_color, structure and chain A and resi 448
color fp_color, structure and chain A and resi 457
color fp_color, structure and chain A and resi 460
color fp_color, structure and chain A and resi 462
color fp_color, structure and chain A and resi 463
color fp_color, structure and chain A and resi 477
color fp_color, structure and chain A and resi 478
color fp_color, structure and chain A and resi 482
color fp_color, structure and chain A and resi 502
# TN (low CSP -- Allosteric): 13 residues
color tn_color, structure and chain A and resi 419
color tn_color, structure and chain A and resi 423
color tn_color, structure and chain A and resi 430
color tn_color, structure and chain A and resi 434
color tn_color, structure and chain A and resi 440
color tn_color, structure and chain A and resi 458
color tn_color, structure and chain A and resi 461
color tn_color, structure and chain A and resi 465
color tn_color, structure and chain A and resi 466
color tn_color, structure and chain A and resi 470
color tn_color, structure and chain A and resi 471
color tn_color, structure and chain A and resi 472
color tn_color, structure and chain A and resi 474
# FN (low CSP in Binding Site): 15 residues
color fn_color, structure and chain A and resi 421
color fn_color, structure and chain A and resi 425
color fn_color, structure and chain A and resi 426
color fn_color, structure and chain A and resi 428
color fn_color, structure and chain A and resi 429
color fn_color, structure and chain A and resi 436
color fn_color, structure and chain A and resi 449
color fn_color, structure and chain A and resi 451
color fn_color, structure and chain A and resi 452
color fn_color, structure and chain A and resi 453
color fn_color, structure and chain A and resi 454
color fn_color, structure and chain A and resi 455
color fn_color, structure and chain A and resi 456
color fn_color, structure and chain A and resi 459
color fn_color, structure and chain A and resi 475
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
# TP (Sig. CSP in Union Site): 9
# FP (Sig. CSP -- Allosteric): 17
# TN (low CSP -- Allosteric): 13
# FN (low CSP in Union Site): 15
# Residues without CSP data: 20
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/9C5E_1/9C5E_csp.pdb
