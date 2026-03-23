reinitialize
load ./outputs/6OQJ_2/6OQJ_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 6 residues
color tp_color, structure and chain B and resi 121
color tp_color, structure and chain B and resi 128
color tp_color, structure and chain B and resi 129
color tp_color, structure and chain B and resi 133
color tp_color, structure and chain B and resi 136
color tp_color, structure and chain B and resi 139
# FP (Sig. CSP -- Allosteric): 16 residues
color fp_color, structure and chain B and resi 115
color fp_color, structure and chain B and resi 116
color fp_color, structure and chain B and resi 119
color fp_color, structure and chain B and resi 123
color fp_color, structure and chain B and resi 124
color fp_color, structure and chain B and resi 127
color fp_color, structure and chain B and resi 131
color fp_color, structure and chain B and resi 134
color fp_color, structure and chain B and resi 135
color fp_color, structure and chain B and resi 137
color fp_color, structure and chain B and resi 138
color fp_color, structure and chain B and resi 142
color fp_color, structure and chain B and resi 143
color fp_color, structure and chain B and resi 144
color fp_color, structure and chain B and resi 145
color fp_color, structure and chain B and resi 146
# TN (low CSP -- Allosteric): 23 residues
color tn_color, structure and chain B and resi 103
color tn_color, structure and chain B and resi 104
color tn_color, structure and chain B and resi 105
color tn_color, structure and chain B and resi 106
color tn_color, structure and chain B and resi 107
color tn_color, structure and chain B and resi 108
color tn_color, structure and chain B and resi 109
color tn_color, structure and chain B and resi 110
color tn_color, structure and chain B and resi 111
color tn_color, structure and chain B and resi 112
color tn_color, structure and chain B and resi 113
color tn_color, structure and chain B and resi 114
color tn_color, structure and chain B and resi 118
color tn_color, structure and chain B and resi 120
color tn_color, structure and chain B and resi 122
color tn_color, structure and chain B and resi 125
color tn_color, structure and chain B and resi 140
color tn_color, structure and chain B and resi 141
color tn_color, structure and chain B and resi 147
color tn_color, structure and chain B and resi 148
color tn_color, structure and chain B and resi 150
color tn_color, structure and chain B and resi 151
color tn_color, structure and chain B and resi 152
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain B and resi 117
color fn_color, structure and chain B and resi 126
color fn_color, structure and chain B and resi 130
color fn_color, structure and chain B and resi 132
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
# TP (Sig. CSP in Union Site): 6
# FP (Sig. CSP -- Allosteric): 16
# TN (low CSP -- Allosteric): 23
# FN (low CSP in Union Site): 4
# Residues without CSP data: 3
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/6OQJ_2/6OQJ_csp.pdb
