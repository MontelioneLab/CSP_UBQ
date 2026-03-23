reinitialize
load ./outputs/1cf4/1cf4_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 13 residues
color tp_color, structure and chain A and resi 7
color tp_color, structure and chain A and resi 9
color tp_color, structure and chain A and resi 17
color tp_color, structure and chain A and resi 36
color tp_color, structure and chain A and resi 42
color tp_color, structure and chain A and resi 43
color tp_color, structure and chain A and resi 44
color tp_color, structure and chain A and resi 45
color tp_color, structure and chain A and resi 46
color tp_color, structure and chain A and resi 70
color tp_color, structure and chain A and resi 173
color tp_color, structure and chain A and resi 174
color tp_color, structure and chain A and resi 177
# FP (Sig. CSP -- Allosteric): 11 residues
color fp_color, structure and chain A and resi 3
color fp_color, structure and chain A and resi 4
color fp_color, structure and chain A and resi 20
color fp_color, structure and chain A and resi 21
color fp_color, structure and chain A and resi 24
color fp_color, structure and chain A and resi 25
color fp_color, structure and chain A and resi 52
color fp_color, structure and chain A and resi 55
color fp_color, structure and chain A and resi 111
color fp_color, structure and chain A and resi 112
color fp_color, structure and chain A and resi 138
# TN (low CSP -- Allosteric): 33 residues
color tn_color, structure and chain A and resi 19
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 75
color tn_color, structure and chain A and resi 77
color tn_color, structure and chain A and resi 79
color tn_color, structure and chain A and resi 80
color tn_color, structure and chain A and resi 84
color tn_color, structure and chain A and resi 93
color tn_color, structure and chain A and resi 98
color tn_color, structure and chain A and resi 101
color tn_color, structure and chain A and resi 102
color tn_color, structure and chain A and resi 108
color tn_color, structure and chain A and resi 117
color tn_color, structure and chain A and resi 119
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 126
color tn_color, structure and chain A and resi 129
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 137
color tn_color, structure and chain A and resi 141
color tn_color, structure and chain A and resi 142
color tn_color, structure and chain A and resi 145
color tn_color, structure and chain A and resi 146
color tn_color, structure and chain A and resi 149
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 155
color tn_color, structure and chain A and resi 161
color tn_color, structure and chain A and resi 165
color tn_color, structure and chain A and resi 168
color tn_color, structure and chain A and resi 172
color tn_color, structure and chain A and resi 175
color tn_color, structure and chain A and resi 176
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain A and resi 8
color fn_color, structure and chain A and resi 85
color fn_color, structure and chain A and resi 113
color fn_color, structure and chain A and resi 115
color fn_color, structure and chain A and resi 159
color fn_color, structure and chain A and resi 160
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
# TP (Sig. CSP in Union Site): 13
# FP (Sig. CSP -- Allosteric): 11
# TN (low CSP -- Allosteric): 33
# FN (low CSP in Union Site): 6
# Residues without CSP data: 121
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/1cf4/1cf4_csp.pdb
