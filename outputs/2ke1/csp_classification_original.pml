reinitialize
load ./outputs/2ke1/2ke1_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 19 residues
color tp_color, structure and chain A and resi 295
color tp_color, structure and chain A and resi 296
color tp_color, structure and chain A and resi 297
color tp_color, structure and chain A and resi 298
color tp_color, structure and chain A and resi 299
color tp_color, structure and chain A and resi 304
color tp_color, structure and chain A and resi 305
color tp_color, structure and chain A and resi 306
color tp_color, structure and chain A and resi 308
color tp_color, structure and chain A and resi 309
color tp_color, structure and chain A and resi 310
color tp_color, structure and chain A and resi 312
color tp_color, structure and chain A and resi 314
color tp_color, structure and chain A and resi 319
color tp_color, structure and chain A and resi 330
color tp_color, structure and chain A and resi 332
color tp_color, structure and chain A and resi 333
color tp_color, structure and chain A and resi 334
color tp_color, structure and chain A and resi 337
# FP (Sig. CSP -- Allosteric): 8 residues
color fp_color, structure and chain A and resi 293
color fp_color, structure and chain A and resi 294
color fp_color, structure and chain A and resi 300
color fp_color, structure and chain A and resi 313
color fp_color, structure and chain A and resi 316
color fp_color, structure and chain A and resi 317
color fp_color, structure and chain A and resi 320
color fp_color, structure and chain A and resi 336
# TN (low CSP -- Allosteric): 25 residues
color tn_color, structure and chain A and resi 291
color tn_color, structure and chain A and resi 292
color tn_color, structure and chain A and resi 303
color tn_color, structure and chain A and resi 318
color tn_color, structure and chain A and resi 321
color tn_color, structure and chain A and resi 323
color tn_color, structure and chain A and resi 324
color tn_color, structure and chain A and resi 327
color tn_color, structure and chain A and resi 328
color tn_color, structure and chain A and resi 329
color tn_color, structure and chain A and resi 338
color tn_color, structure and chain A and resi 339
color tn_color, structure and chain A and resi 341
color tn_color, structure and chain A and resi 342
color tn_color, structure and chain A and resi 343
color tn_color, structure and chain A and resi 344
color tn_color, structure and chain A and resi 345
color tn_color, structure and chain A and resi 346
color tn_color, structure and chain A and resi 347
color tn_color, structure and chain A and resi 348
color tn_color, structure and chain A and resi 349
color tn_color, structure and chain A and resi 351
color tn_color, structure and chain A and resi 352
color tn_color, structure and chain A and resi 353
color tn_color, structure and chain A and resi 354
# FN (low CSP in Binding Site): 7 residues
color fn_color, structure and chain A and resi 301
color fn_color, structure and chain A and resi 302
color fn_color, structure and chain A and resi 307
color fn_color, structure and chain A and resi 311
color fn_color, structure and chain A and resi 322
color fn_color, structure and chain A and resi 335
color fn_color, structure and chain A and resi 340
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
# TP (Sig. CSP in Union Site): 19
# FP (Sig. CSP -- Allosteric): 8
# TN (low CSP -- Allosteric): 25
# FN (low CSP in Union Site): 7
# Residues without CSP data: 7
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2ke1/2ke1_csp.pdb
