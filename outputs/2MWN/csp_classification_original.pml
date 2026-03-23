reinitialize
load ./outputs/2MWN/2MWN_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 3 residues
color tp_color, structure and chain B and resi 316
color tp_color, structure and chain B and resi 324
color tp_color, structure and chain B and resi 325
# FP (Sig. CSP -- Allosteric): 20 residues
color fp_color, structure and chain B and resi 310
color fp_color, structure and chain B and resi 311
color fp_color, structure and chain B and resi 312
color fp_color, structure and chain B and resi 313
color fp_color, structure and chain B and resi 314
color fp_color, structure and chain B and resi 315
color fp_color, structure and chain B and resi 317
color fp_color, structure and chain B and resi 318
color fp_color, structure and chain B and resi 319
color fp_color, structure and chain B and resi 320
color fp_color, structure and chain B and resi 321
color fp_color, structure and chain B and resi 322
color fp_color, structure and chain B and resi 326
color fp_color, structure and chain B and resi 328
color fp_color, structure and chain B and resi 329
color fp_color, structure and chain B and resi 364
color fp_color, structure and chain B and resi 382
color fp_color, structure and chain B and resi 393
color fp_color, structure and chain B and resi 398
color fp_color, structure and chain B and resi 400
# TN (low CSP -- Allosteric): 56 residues
color tn_color, structure and chain B and resi 330
color tn_color, structure and chain B and resi 331
color tn_color, structure and chain B and resi 332
color tn_color, structure and chain B and resi 333
color tn_color, structure and chain B and resi 334
color tn_color, structure and chain B and resi 335
color tn_color, structure and chain B and resi 336
color tn_color, structure and chain B and resi 337
color tn_color, structure and chain B and resi 338
color tn_color, structure and chain B and resi 339
color tn_color, structure and chain B and resi 340
color tn_color, structure and chain B and resi 341
color tn_color, structure and chain B and resi 342
color tn_color, structure and chain B and resi 343
color tn_color, structure and chain B and resi 344
color tn_color, structure and chain B and resi 345
color tn_color, structure and chain B and resi 346
color tn_color, structure and chain B and resi 347
color tn_color, structure and chain B and resi 348
color tn_color, structure and chain B and resi 349
color tn_color, structure and chain B and resi 350
color tn_color, structure and chain B and resi 351
color tn_color, structure and chain B and resi 352
color tn_color, structure and chain B and resi 353
color tn_color, structure and chain B and resi 354
color tn_color, structure and chain B and resi 355
color tn_color, structure and chain B and resi 356
color tn_color, structure and chain B and resi 357
color tn_color, structure and chain B and resi 359
color tn_color, structure and chain B and resi 361
color tn_color, structure and chain B and resi 366
color tn_color, structure and chain B and resi 368
color tn_color, structure and chain B and resi 370
color tn_color, structure and chain B and resi 371
color tn_color, structure and chain B and resi 372
color tn_color, structure and chain B and resi 373
color tn_color, structure and chain B and resi 374
color tn_color, structure and chain B and resi 375
color tn_color, structure and chain B and resi 376
color tn_color, structure and chain B and resi 377
color tn_color, structure and chain B and resi 378
color tn_color, structure and chain B and resi 383
color tn_color, structure and chain B and resi 384
color tn_color, structure and chain B and resi 385
color tn_color, structure and chain B and resi 386
color tn_color, structure and chain B and resi 387
color tn_color, structure and chain B and resi 388
color tn_color, structure and chain B and resi 389
color tn_color, structure and chain B and resi 390
color tn_color, structure and chain B and resi 391
color tn_color, structure and chain B and resi 392
color tn_color, structure and chain B and resi 394
color tn_color, structure and chain B and resi 395
color tn_color, structure and chain B and resi 396
color tn_color, structure and chain B and resi 397
color tn_color, structure and chain B and resi 399
# FN (low CSP in Binding Site): 9 residues
color fn_color, structure and chain B and resi 358
color fn_color, structure and chain B and resi 360
color fn_color, structure and chain B and resi 362
color fn_color, structure and chain B and resi 365
color fn_color, structure and chain B and resi 367
color fn_color, structure and chain B and resi 369
color fn_color, structure and chain B and resi 379
color fn_color, structure and chain B and resi 380
color fn_color, structure and chain B and resi 381
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
# TP (Sig. CSP in Union Site): 3
# FP (Sig. CSP -- Allosteric): 20
# TN (low CSP -- Allosteric): 56
# FN (low CSP in Union Site): 9
# Residues without CSP data: 5
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2MWN/2MWN_csp.pdb
