reinitialize
load ./outputs/7SFT/7SFT_csp.pdb, structure
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
color tp_color, structure and chain A and resi 2270
color tp_color, structure and chain A and resi 2271
color tp_color, structure and chain A and resi 2272
color tp_color, structure and chain A and resi 2274
color tp_color, structure and chain A and resi 2275
color tp_color, structure and chain A and resi 2276
color tp_color, structure and chain A and resi 2277
color tp_color, structure and chain A and resi 2280
# FP (Sig. CSP -- Allosteric): 22 residues
color fp_color, structure and chain A and resi 2236
color fp_color, structure and chain A and resi 2237
color fp_color, structure and chain A and resi 2238
color fp_color, structure and chain A and resi 2239
color fp_color, structure and chain A and resi 2240
color fp_color, structure and chain A and resi 2241
color fp_color, structure and chain A and resi 2242
color fp_color, structure and chain A and resi 2243
color fp_color, structure and chain A and resi 2244
color fp_color, structure and chain A and resi 2246
color fp_color, structure and chain A and resi 2247
color fp_color, structure and chain A and resi 2248
color fp_color, structure and chain A and resi 2260
color fp_color, structure and chain A and resi 2281
color fp_color, structure and chain A and resi 2284
color fp_color, structure and chain A and resi 2290
color fp_color, structure and chain A and resi 2292
color fp_color, structure and chain A and resi 2300
color fp_color, structure and chain A and resi 2303
color fp_color, structure and chain A and resi 2310
color fp_color, structure and chain A and resi 2312
color fp_color, structure and chain A and resi 2329
# TN (low CSP -- Allosteric): 32 residues
color tn_color, structure and chain A and resi 2251
color tn_color, structure and chain A and resi 2252
color tn_color, structure and chain A and resi 2253
color tn_color, structure and chain A and resi 2254
color tn_color, structure and chain A and resi 2255
color tn_color, structure and chain A and resi 2257
color tn_color, structure and chain A and resi 2258
color tn_color, structure and chain A and resi 2262
color tn_color, structure and chain A and resi 2266
color tn_color, structure and chain A and resi 2282
color tn_color, structure and chain A and resi 2295
color tn_color, structure and chain A and resi 2296
color tn_color, structure and chain A and resi 2297
color tn_color, structure and chain A and resi 2298
color tn_color, structure and chain A and resi 2299
color tn_color, structure and chain A and resi 2305
color tn_color, structure and chain A and resi 2306
color tn_color, structure and chain A and resi 2307
color tn_color, structure and chain A and resi 2308
color tn_color, structure and chain A and resi 2309
color tn_color, structure and chain A and resi 2313
color tn_color, structure and chain A and resi 2314
color tn_color, structure and chain A and resi 2315
color tn_color, structure and chain A and resi 2316
color tn_color, structure and chain A and resi 2318
color tn_color, structure and chain A and resi 2319
color tn_color, structure and chain A and resi 2321
color tn_color, structure and chain A and resi 2322
color tn_color, structure and chain A and resi 2323
color tn_color, structure and chain A and resi 2325
color tn_color, structure and chain A and resi 2326
color tn_color, structure and chain A and resi 2327
# FN (low CSP in Binding Site): 10 residues
color fn_color, structure and chain A and resi 2263
color fn_color, structure and chain A and resi 2267
color fn_color, structure and chain A and resi 2268
color fn_color, structure and chain A and resi 2269
color fn_color, structure and chain A and resi 2273
color fn_color, structure and chain A and resi 2279
color fn_color, structure and chain A and resi 2285
color fn_color, structure and chain A and resi 2286
color fn_color, structure and chain A and resi 2291
color fn_color, structure and chain A and resi 2311
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
# FP (Sig. CSP -- Allosteric): 22
# TN (low CSP -- Allosteric): 32
# FN (low CSP in Union Site): 10
# Residues without CSP data: 23
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/7SFT/7SFT_csp.pdb
