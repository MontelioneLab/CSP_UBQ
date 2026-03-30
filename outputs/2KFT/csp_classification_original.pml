reinitialize
load ./outputs/2KFT/2KFT_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 14 residues
color tp_color, structure and chain A and resi 297
color tp_color, structure and chain A and resi 298
color tp_color, structure and chain A and resi 299
color tp_color, structure and chain A and resi 302
color tp_color, structure and chain A and resi 304
color tp_color, structure and chain A and resi 305
color tp_color, structure and chain A and resi 306
color tp_color, structure and chain A and resi 311
color tp_color, structure and chain A and resi 312
color tp_color, structure and chain A and resi 314
color tp_color, structure and chain A and resi 333
color tp_color, structure and chain A and resi 335
color tp_color, structure and chain A and resi 337
color tp_color, structure and chain A and resi 340
# FP (Sig. CSP -- Allosteric): 14 residues
color fp_color, structure and chain A and resi 296
color fp_color, structure and chain A and resi 300
color fp_color, structure and chain A and resi 303
color fp_color, structure and chain A and resi 316
color fp_color, structure and chain A and resi 318
color fp_color, structure and chain A and resi 320
color fp_color, structure and chain A and resi 329
color fp_color, structure and chain A and resi 332
color fp_color, structure and chain A and resi 336
color fp_color, structure and chain A and resi 341
color fp_color, structure and chain A and resi 344
color fp_color, structure and chain A and resi 345
color fp_color, structure and chain A and resi 346
color fp_color, structure and chain A and resi 347
# TN (low CSP -- Allosteric): 11 residues
color tn_color, structure and chain A and resi 313
color tn_color, structure and chain A and resi 317
color tn_color, structure and chain A and resi 321
color tn_color, structure and chain A and resi 323
color tn_color, structure and chain A and resi 324
color tn_color, structure and chain A and resi 327
color tn_color, structure and chain A and resi 328
color tn_color, structure and chain A and resi 338
color tn_color, structure and chain A and resi 339
color tn_color, structure and chain A and resi 342
color tn_color, structure and chain A and resi 343
# FN (low CSP in Binding Site): 7 residues
color fn_color, structure and chain A and resi 301
color fn_color, structure and chain A and resi 307
color fn_color, structure and chain A and resi 308
color fn_color, structure and chain A and resi 309
color fn_color, structure and chain A and resi 319
color fn_color, structure and chain A and resi 322
color fn_color, structure and chain A and resi 330
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
# TP (Sig. CSP in Union Site): 14
# FP (Sig. CSP -- Allosteric): 14
# TN (low CSP -- Allosteric): 11
# FN (low CSP in Union Site): 7
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2KFT/2KFT_csp.pdb
