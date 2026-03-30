reinitialize
load ./outputs/2law/2law_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 6 residues
color tp_color, structure and chain A and resi 239
color tp_color, structure and chain A and resi 247
color tp_color, structure and chain A and resi 249
color tp_color, structure and chain A and resi 251
color tp_color, structure and chain A and resi 254
color tp_color, structure and chain A and resi 255
# FP (Sig. CSP -- Allosteric): 5 residues
color fp_color, structure and chain A and resi 243
color fp_color, structure and chain A and resi 248
color fp_color, structure and chain A and resi 250
color fp_color, structure and chain A and resi 259
color fp_color, structure and chain A and resi 262
# TN (low CSP -- Allosteric): 15 residues
color tn_color, structure and chain A and resi 230
color tn_color, structure and chain A and resi 232
color tn_color, structure and chain A and resi 234
color tn_color, structure and chain A and resi 235
color tn_color, structure and chain A and resi 236
color tn_color, structure and chain A and resi 238
color tn_color, structure and chain A and resi 240
color tn_color, structure and chain A and resi 244
color tn_color, structure and chain A and resi 245
color tn_color, structure and chain A and resi 246
color tn_color, structure and chain A and resi 252
color tn_color, structure and chain A and resi 253
color tn_color, structure and chain A and resi 257
color tn_color, structure and chain A and resi 260
color tn_color, structure and chain A and resi 263
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain A and resi 237
color fn_color, structure and chain A and resi 241
color fn_color, structure and chain A and resi 256
color fn_color, structure and chain A and resi 258
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
# TP (Sig. CSP in Union Site): 6
# FP (Sig. CSP -- Allosteric): 5
# TN (low CSP -- Allosteric): 15
# FN (low CSP in Union Site): 4
# Residues without CSP data: 4
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2law/2law_csp.pdb
