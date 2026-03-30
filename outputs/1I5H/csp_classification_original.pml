reinitialize
load ./outputs/1I5H/1I5H_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain W
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 3 residues
color tp_color, structure and chain W and resi 483
color tp_color, structure and chain W and resi 485
color tp_color, structure and chain W and resi 486
# FP (Sig. CSP -- Allosteric): 5 residues
color fp_color, structure and chain W and resi 472
color fp_color, structure and chain W and resi 481
color fp_color, structure and chain W and resi 484
color fp_color, structure and chain W and resi 488
color fp_color, structure and chain W and resi 494
# TN (low CSP -- Allosteric): 18 residues
color tn_color, structure and chain W and resi 458
color tn_color, structure and chain W and resi 459
color tn_color, structure and chain W and resi 461
color tn_color, structure and chain W and resi 464
color tn_color, structure and chain W and resi 465
color tn_color, structure and chain W and resi 467
color tn_color, structure and chain W and resi 468
color tn_color, structure and chain W and resi 469
color tn_color, structure and chain W and resi 473
color tn_color, structure and chain W and resi 474
color tn_color, structure and chain W and resi 475
color tn_color, structure and chain W and resi 477
color tn_color, structure and chain W and resi 479
color tn_color, structure and chain W and resi 480
color tn_color, structure and chain W and resi 489
color tn_color, structure and chain W and resi 491
color tn_color, structure and chain W and resi 492
color tn_color, structure and chain W and resi 493
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain W and resi 470
color fn_color, structure and chain W and resi 476
color fn_color, structure and chain W and resi 478
color fn_color, structure and chain W and resi 487
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
# Receptor chain: W
# Ligand chain: B
# TP (Sig. CSP in Union Site): 3
# FP (Sig. CSP -- Allosteric): 5
# TN (low CSP -- Allosteric): 18
# FN (low CSP in Union Site): 4
# Residues without CSP data: 20
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/1I5H/1I5H_csp.pdb
