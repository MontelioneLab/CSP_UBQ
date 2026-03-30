reinitialize
load ./outputs/2RR4/2RR4_csp.pdb, structure
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
color tp_color, structure and chain A and resi 252
color tp_color, structure and chain A and resi 255
color tp_color, structure and chain A and resi 256
color tp_color, structure and chain A and resi 257
color tp_color, structure and chain A and resi 258
color tp_color, structure and chain A and resi 259
color tp_color, structure and chain A and resi 260
color tp_color, structure and chain A and resi 267
color tp_color, structure and chain A and resi 281
color tp_color, structure and chain A and resi 283
color tp_color, structure and chain A and resi 296
color tp_color, structure and chain A and resi 303
color tp_color, structure and chain A and resi 306
color tp_color, structure and chain A and resi 307
# FP (Sig. CSP -- Allosteric): 15 residues
color fp_color, structure and chain A and resi 261
color fp_color, structure and chain A and resi 268
color fp_color, structure and chain A and resi 269
color fp_color, structure and chain A and resi 270
color fp_color, structure and chain A and resi 271
color fp_color, structure and chain A and resi 272
color fp_color, structure and chain A and resi 273
color fp_color, structure and chain A and resi 274
color fp_color, structure and chain A and resi 275
color fp_color, structure and chain A and resi 278
color fp_color, structure and chain A and resi 279
color fp_color, structure and chain A and resi 289
color fp_color, structure and chain A and resi 295
color fp_color, structure and chain A and resi 300
color fp_color, structure and chain A and resi 301
# TN (low CSP -- Allosteric): 24 residues
color tn_color, structure and chain A and resi 245
color tn_color, structure and chain A and resi 246
color tn_color, structure and chain A and resi 247
color tn_color, structure and chain A and resi 248
color tn_color, structure and chain A and resi 249
color tn_color, structure and chain A and resi 250
color tn_color, structure and chain A and resi 251
color tn_color, structure and chain A and resi 263
color tn_color, structure and chain A and resi 265
color tn_color, structure and chain A and resi 266
color tn_color, structure and chain A and resi 277
color tn_color, structure and chain A and resi 284
color tn_color, structure and chain A and resi 286
color tn_color, structure and chain A and resi 287
color tn_color, structure and chain A and resi 288
color tn_color, structure and chain A and resi 290
color tn_color, structure and chain A and resi 291
color tn_color, structure and chain A and resi 292
color tn_color, structure and chain A and resi 293
color tn_color, structure and chain A and resi 294
color tn_color, structure and chain A and resi 297
color tn_color, structure and chain A and resi 298
color tn_color, structure and chain A and resi 302
color tn_color, structure and chain A and resi 304
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain A and resi 264
color fn_color, structure and chain A and resi 282
color fn_color, structure and chain A and resi 285
color fn_color, structure and chain A and resi 305
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
# FP (Sig. CSP -- Allosteric): 15
# TN (low CSP -- Allosteric): 24
# FN (low CSP in Union Site): 4
# Residues without CSP data: 12
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2RR4/2RR4_csp.pdb
