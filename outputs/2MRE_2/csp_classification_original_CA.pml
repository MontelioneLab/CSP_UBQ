reinitialize
load ./outputs/2MRE_2/2MRE_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 7 residues
color tp_color, structure and chain B and resi 202
color tp_color, structure and chain B and resi 203
color tp_color, structure and chain B and resi 207
color tp_color, structure and chain B and resi 216
color tp_color, structure and chain B and resi 217
color tp_color, structure and chain B and resi 223
color tp_color, structure and chain B and resi 224
# FP (Sig. CSP -- Allosteric): 3 residues
color fp_color, structure and chain B and resi 211
color fp_color, structure and chain B and resi 222
color fp_color, structure and chain B and resi 225
# TN (low CSP -- Allosteric): 9 residues
color tn_color, structure and chain B and resi 199
color tn_color, structure and chain B and resi 200
color tn_color, structure and chain B and resi 201
color tn_color, structure and chain B and resi 208
color tn_color, structure and chain B and resi 209
color tn_color, structure and chain B and resi 210
color tn_color, structure and chain B and resi 214
color tn_color, structure and chain B and resi 215
color tn_color, structure and chain B and resi 218
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain B and resi 204
color fn_color, structure and chain B and resi 206
color fn_color, structure and chain B and resi 213
color fn_color, structure and chain B and resi 219
color fn_color, structure and chain B and resi 220
color fn_color, structure and chain B and resi 221
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
# TP (Sig. CSP in Union Site): 7
# FP (Sig. CSP -- Allosteric): 3
# TN (low CSP -- Allosteric): 9
# FN (low CSP in Union Site): 6
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2MRE_2/2MRE_csp.pdb
