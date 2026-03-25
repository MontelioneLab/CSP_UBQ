reinitialize
load ./outputs/2KWI_2/2KWI_csp.pdb, structure
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
color tp_color, structure and chain B and resi 421
color tp_color, structure and chain B and resi 422
color tp_color, structure and chain B and resi 426
color tp_color, structure and chain B and resi 429
color tp_color, structure and chain B and resi 430
color tp_color, structure and chain B and resi 433
color tp_color, structure and chain B and resi 434
color tp_color, structure and chain B and resi 437
# FP (Sig. CSP -- Allosteric): 14 residues
color fp_color, structure and chain B and resi 411
color fp_color, structure and chain B and resi 412
color fp_color, structure and chain B and resi 428
color fp_color, structure and chain B and resi 431
color fp_color, structure and chain B and resi 432
color fp_color, structure and chain B and resi 435
color fp_color, structure and chain B and resi 436
color fp_color, structure and chain B and resi 438
color fp_color, structure and chain B and resi 439
color fp_color, structure and chain B and resi 441
color fp_color, structure and chain B and resi 442
color fp_color, structure and chain B and resi 443
color fp_color, structure and chain B and resi 444
color fp_color, structure and chain B and resi 445
# TN (low CSP -- Allosteric): 23 residues
color tn_color, structure and chain B and resi 396
color tn_color, structure and chain B and resi 397
color tn_color, structure and chain B and resi 398
color tn_color, structure and chain B and resi 399
color tn_color, structure and chain B and resi 400
color tn_color, structure and chain B and resi 401
color tn_color, structure and chain B and resi 402
color tn_color, structure and chain B and resi 403
color tn_color, structure and chain B and resi 404
color tn_color, structure and chain B and resi 405
color tn_color, structure and chain B and resi 406
color tn_color, structure and chain B and resi 407
color tn_color, structure and chain B and resi 408
color tn_color, structure and chain B and resi 410
color tn_color, structure and chain B and resi 414
color tn_color, structure and chain B and resi 415
color tn_color, structure and chain B and resi 418
color tn_color, structure and chain B and resi 419
color tn_color, structure and chain B and resi 420
color tn_color, structure and chain B and resi 423
color tn_color, structure and chain B and resi 424
color tn_color, structure and chain B and resi 425
color tn_color, structure and chain B and resi 446
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain B and resi 409
color fn_color, structure and chain B and resi 413
color fn_color, structure and chain B and resi 416
color fn_color, structure and chain B and resi 417
color fn_color, structure and chain B and resi 427
color fn_color, structure and chain B and resi 440
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
# FP (Sig. CSP -- Allosteric): 14
# TN (low CSP -- Allosteric): 23
# FN (low CSP in Union Site): 6
# Residues without CSP data: 5
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2KWI_2/2KWI_csp.pdb
