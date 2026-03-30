reinitialize
load ./outputs/2H7D/2H7D_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain A and resi 318
color tp_color, structure and chain A and resi 324
color tp_color, structure and chain A and resi 325
color tp_color, structure and chain A and resi 357
color tp_color, structure and chain A and resi 358
color tp_color, structure and chain A and resi 359
color tp_color, structure and chain A and resi 364
color tp_color, structure and chain A and resi 380
# FP (Sig. CSP -- Allosteric): 31 residues
color fp_color, structure and chain A and resi 306
color fp_color, structure and chain A and resi 307
color fp_color, structure and chain A and resi 308
color fp_color, structure and chain A and resi 309
color fp_color, structure and chain A and resi 310
color fp_color, structure and chain A and resi 311
color fp_color, structure and chain A and resi 312
color fp_color, structure and chain A and resi 313
color fp_color, structure and chain A and resi 314
color fp_color, structure and chain A and resi 315
color fp_color, structure and chain A and resi 316
color fp_color, structure and chain A and resi 317
color fp_color, structure and chain A and resi 319
color fp_color, structure and chain A and resi 320
color fp_color, structure and chain A and resi 321
color fp_color, structure and chain A and resi 322
color fp_color, structure and chain A and resi 326
color fp_color, structure and chain A and resi 328
color fp_color, structure and chain A and resi 329
color fp_color, structure and chain A and resi 334
color fp_color, structure and chain A and resi 344
color fp_color, structure and chain A and resi 345
color fp_color, structure and chain A and resi 361
color fp_color, structure and chain A and resi 370
color fp_color, structure and chain A and resi 374
color fp_color, structure and chain A and resi 393
color fp_color, structure and chain A and resi 394
color fp_color, structure and chain A and resi 395
color fp_color, structure and chain A and resi 398
color fp_color, structure and chain A and resi 400
color fp_color, structure and chain A and resi 404
# TN (low CSP -- Allosteric): 51 residues
color tn_color, structure and chain A and resi 330
color tn_color, structure and chain A and resi 331
color tn_color, structure and chain A and resi 332
color tn_color, structure and chain A and resi 333
color tn_color, structure and chain A and resi 335
color tn_color, structure and chain A and resi 336
color tn_color, structure and chain A and resi 337
color tn_color, structure and chain A and resi 338
color tn_color, structure and chain A and resi 339
color tn_color, structure and chain A and resi 340
color tn_color, structure and chain A and resi 341
color tn_color, structure and chain A and resi 342
color tn_color, structure and chain A and resi 343
color tn_color, structure and chain A and resi 346
color tn_color, structure and chain A and resi 347
color tn_color, structure and chain A and resi 348
color tn_color, structure and chain A and resi 349
color tn_color, structure and chain A and resi 350
color tn_color, structure and chain A and resi 351
color tn_color, structure and chain A and resi 352
color tn_color, structure and chain A and resi 353
color tn_color, structure and chain A and resi 355
color tn_color, structure and chain A and resi 356
color tn_color, structure and chain A and resi 365
color tn_color, structure and chain A and resi 366
color tn_color, structure and chain A and resi 367
color tn_color, structure and chain A and resi 371
color tn_color, structure and chain A and resi 372
color tn_color, structure and chain A and resi 373
color tn_color, structure and chain A and resi 375
color tn_color, structure and chain A and resi 376
color tn_color, structure and chain A and resi 377
color tn_color, structure and chain A and resi 378
color tn_color, structure and chain A and resi 382
color tn_color, structure and chain A and resi 383
color tn_color, structure and chain A and resi 384
color tn_color, structure and chain A and resi 385
color tn_color, structure and chain A and resi 386
color tn_color, structure and chain A and resi 387
color tn_color, structure and chain A and resi 388
color tn_color, structure and chain A and resi 389
color tn_color, structure and chain A and resi 390
color tn_color, structure and chain A and resi 391
color tn_color, structure and chain A and resi 392
color tn_color, structure and chain A and resi 396
color tn_color, structure and chain A and resi 397
color tn_color, structure and chain A and resi 399
color tn_color, structure and chain A and resi 401
color tn_color, structure and chain A and resi 402
color tn_color, structure and chain A and resi 403
color tn_color, structure and chain A and resi 405
# FN (low CSP in Binding Site): 7 residues
color fn_color, structure and chain A and resi 354
color fn_color, structure and chain A and resi 360
color fn_color, structure and chain A and resi 362
color fn_color, structure and chain A and resi 368
color fn_color, structure and chain A and resi 369
color fn_color, structure and chain A and resi 379
color fn_color, structure and chain A and resi 381
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 31
# TN (low CSP -- Allosteric): 51
# FN (low CSP in Union Site): 7
# Residues without CSP data: 4
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2H7D/2H7D_csp.pdb
