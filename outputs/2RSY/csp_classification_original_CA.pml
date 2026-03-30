reinitialize
load ./outputs/2RSY/2RSY_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 12 residues
color tp_color, structure and chain A and resi 88
color tp_color, structure and chain A and resi 89
color tp_color, structure and chain A and resi 109
color tp_color, structure and chain A and resi 110
color tp_color, structure and chain A and resi 111
color tp_color, structure and chain A and resi 123
color tp_color, structure and chain A and resi 126
color tp_color, structure and chain A and resi 128
color tp_color, structure and chain A and resi 129
color tp_color, structure and chain A and resi 130
color tp_color, structure and chain A and resi 159
color tp_color, structure and chain A and resi 163
# FP (Sig. CSP -- Allosteric): 34 residues
color fp_color, structure and chain A and resi 77
color fp_color, structure and chain A and resi 78
color fp_color, structure and chain A and resi 79
color fp_color, structure and chain A and resi 80
color fp_color, structure and chain A and resi 85
color fp_color, structure and chain A and resi 86
color fp_color, structure and chain A and resi 87
color fp_color, structure and chain A and resi 90
color fp_color, structure and chain A and resi 92
color fp_color, structure and chain A and resi 93
color fp_color, structure and chain A and resi 94
color fp_color, structure and chain A and resi 95
color fp_color, structure and chain A and resi 97
color fp_color, structure and chain A and resi 102
color fp_color, structure and chain A and resi 108
color fp_color, structure and chain A and resi 112
color fp_color, structure and chain A and resi 114
color fp_color, structure and chain A and resi 120
color fp_color, structure and chain A and resi 121
color fp_color, structure and chain A and resi 122
color fp_color, structure and chain A and resi 124
color fp_color, structure and chain A and resi 133
color fp_color, structure and chain A and resi 134
color fp_color, structure and chain A and resi 136
color fp_color, structure and chain A and resi 144
color fp_color, structure and chain A and resi 155
color fp_color, structure and chain A and resi 157
color fp_color, structure and chain A and resi 158
color fp_color, structure and chain A and resi 164
color fp_color, structure and chain A and resi 165
color fp_color, structure and chain A and resi 167
color fp_color, structure and chain A and resi 168
color fp_color, structure and chain A and resi 172
color fp_color, structure and chain A and resi 173
# TN (low CSP -- Allosteric): 37 residues
color tn_color, structure and chain A and resi 82
color tn_color, structure and chain A and resi 83
color tn_color, structure and chain A and resi 84
color tn_color, structure and chain A and resi 91
color tn_color, structure and chain A and resi 96
color tn_color, structure and chain A and resi 100
color tn_color, structure and chain A and resi 101
color tn_color, structure and chain A and resi 103
color tn_color, structure and chain A and resi 104
color tn_color, structure and chain A and resi 105
color tn_color, structure and chain A and resi 106
color tn_color, structure and chain A and resi 115
color tn_color, structure and chain A and resi 116
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 119
color tn_color, structure and chain A and resi 131
color tn_color, structure and chain A and resi 132
color tn_color, structure and chain A and resi 137
color tn_color, structure and chain A and resi 138
color tn_color, structure and chain A and resi 139
color tn_color, structure and chain A and resi 142
color tn_color, structure and chain A and resi 143
color tn_color, structure and chain A and resi 145
color tn_color, structure and chain A and resi 146
color tn_color, structure and chain A and resi 147
color tn_color, structure and chain A and resi 148
color tn_color, structure and chain A and resi 149
color tn_color, structure and chain A and resi 150
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 153
color tn_color, structure and chain A and resi 154
color tn_color, structure and chain A and resi 156
color tn_color, structure and chain A and resi 160
color tn_color, structure and chain A and resi 166
color tn_color, structure and chain A and resi 169
color tn_color, structure and chain A and resi 171
# FN (low CSP in Binding Site): 8 residues
color fn_color, structure and chain A and resi 107
color fn_color, structure and chain A and resi 117
color fn_color, structure and chain A and resi 125
color fn_color, structure and chain A and resi 127
color fn_color, structure and chain A and resi 140
color fn_color, structure and chain A and resi 141
color fn_color, structure and chain A and resi 161
color fn_color, structure and chain A and resi 162
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
# TP (Sig. CSP in Union Site): 12
# FP (Sig. CSP -- Allosteric): 34
# TN (low CSP -- Allosteric): 37
# FN (low CSP in Union Site): 8
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2RSY/2RSY_csp.pdb
