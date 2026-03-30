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
# TP (Sig. CSP in Binding Site): 2 residues
color tp_color, structure and chain A and resi 2269
color tp_color, structure and chain A and resi 2270
# FP (Sig. CSP -- Allosteric): 15 residues
color fp_color, structure and chain A and resi 2237
color fp_color, structure and chain A and resi 2242
color fp_color, structure and chain A and resi 2243
color fp_color, structure and chain A and resi 2244
color fp_color, structure and chain A and resi 2246
color fp_color, structure and chain A and resi 2255
color fp_color, structure and chain A and resi 2257
color fp_color, structure and chain A and resi 2260
color fp_color, structure and chain A and resi 2266
color fp_color, structure and chain A and resi 2281
color fp_color, structure and chain A and resi 2290
color fp_color, structure and chain A and resi 2297
color fp_color, structure and chain A and resi 2298
color fp_color, structure and chain A and resi 2303
color fp_color, structure and chain A and resi 2329
# TN (low CSP -- Allosteric): 21 residues
color tn_color, structure and chain A and resi 2236
color tn_color, structure and chain A and resi 2238
color tn_color, structure and chain A and resi 2239
color tn_color, structure and chain A and resi 2247
color tn_color, structure and chain A and resi 2248
color tn_color, structure and chain A and resi 2251
color tn_color, structure and chain A and resi 2252
color tn_color, structure and chain A and resi 2253
color tn_color, structure and chain A and resi 2254
color tn_color, structure and chain A and resi 2292
color tn_color, structure and chain A and resi 2295
color tn_color, structure and chain A and resi 2299
color tn_color, structure and chain A and resi 2306
color tn_color, structure and chain A and resi 2307
color tn_color, structure and chain A and resi 2312
color tn_color, structure and chain A and resi 2313
color tn_color, structure and chain A and resi 2318
color tn_color, structure and chain A and resi 2323
color tn_color, structure and chain A and resi 2325
color tn_color, structure and chain A and resi 2326
color tn_color, structure and chain A and resi 2327
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain A and resi 2267
color fn_color, structure and chain A and resi 2271
color fn_color, structure and chain A and resi 2277
color fn_color, structure and chain A and resi 2291
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
# TP (Sig. CSP in Union Site): 2
# FP (Sig. CSP -- Allosteric): 15
# TN (low CSP -- Allosteric): 21
# FN (low CSP in Union Site): 4
# Residues without CSP data: 53
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/7SFT/7SFT_csp.pdb
