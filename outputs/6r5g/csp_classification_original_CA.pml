reinitialize
load ./outputs/6r5g/6r5g_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 15 residues
color tp_color, structure and chain A and resi 119
color tp_color, structure and chain A and resi 140
color tp_color, structure and chain A and resi 148
color tp_color, structure and chain A and resi 169
color tp_color, structure and chain A and resi 170
color tp_color, structure and chain A and resi 171
color tp_color, structure and chain A and resi 181
color tp_color, structure and chain A and resi 182
color tp_color, structure and chain A and resi 183
color tp_color, structure and chain A and resi 184
color tp_color, structure and chain A and resi 197
color tp_color, structure and chain A and resi 203
color tp_color, structure and chain A and resi 204
color tp_color, structure and chain A and resi 205
color tp_color, structure and chain A and resi 210
# FP (Sig. CSP -- Allosteric): 23 residues
color fp_color, structure and chain A and resi 108
color fp_color, structure and chain A and resi 116
color fp_color, structure and chain A and resi 117
color fp_color, structure and chain A and resi 122
color fp_color, structure and chain A and resi 133
color fp_color, structure and chain A and resi 139
color fp_color, structure and chain A and resi 143
color fp_color, structure and chain A and resi 145
color fp_color, structure and chain A and resi 146
color fp_color, structure and chain A and resi 167
color fp_color, structure and chain A and resi 172
color fp_color, structure and chain A and resi 173
color fp_color, structure and chain A and resi 175
color fp_color, structure and chain A and resi 179
color fp_color, structure and chain A and resi 180
color fp_color, structure and chain A and resi 185
color fp_color, structure and chain A and resi 198
color fp_color, structure and chain A and resi 199
color fp_color, structure and chain A and resi 208
color fp_color, structure and chain A and resi 209
color fp_color, structure and chain A and resi 211
color fp_color, structure and chain A and resi 212
color fp_color, structure and chain A and resi 213
# TN (low CSP -- Allosteric): 66 residues
color tn_color, structure and chain A and resi 104
color tn_color, structure and chain A and resi 105
color tn_color, structure and chain A and resi 106
color tn_color, structure and chain A and resi 109
color tn_color, structure and chain A and resi 110
color tn_color, structure and chain A and resi 111
color tn_color, structure and chain A and resi 112
color tn_color, structure and chain A and resi 113
color tn_color, structure and chain A and resi 114
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 121
color tn_color, structure and chain A and resi 123
color tn_color, structure and chain A and resi 124
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 126
color tn_color, structure and chain A and resi 127
color tn_color, structure and chain A and resi 128
color tn_color, structure and chain A and resi 129
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 131
color tn_color, structure and chain A and resi 132
color tn_color, structure and chain A and resi 134
color tn_color, structure and chain A and resi 135
color tn_color, structure and chain A and resi 136
color tn_color, structure and chain A and resi 137
color tn_color, structure and chain A and resi 147
color tn_color, structure and chain A and resi 149
color tn_color, structure and chain A and resi 150
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 153
color tn_color, structure and chain A and resi 154
color tn_color, structure and chain A and resi 155
color tn_color, structure and chain A and resi 156
color tn_color, structure and chain A and resi 157
color tn_color, structure and chain A and resi 158
color tn_color, structure and chain A and resi 159
color tn_color, structure and chain A and resi 160
color tn_color, structure and chain A and resi 162
color tn_color, structure and chain A and resi 163
color tn_color, structure and chain A and resi 164
color tn_color, structure and chain A and resi 165
color tn_color, structure and chain A and resi 174
color tn_color, structure and chain A and resi 176
color tn_color, structure and chain A and resi 177
color tn_color, structure and chain A and resi 178
color tn_color, structure and chain A and resi 186
color tn_color, structure and chain A and resi 187
color tn_color, structure and chain A and resi 188
color tn_color, structure and chain A and resi 189
color tn_color, structure and chain A and resi 190
color tn_color, structure and chain A and resi 191
color tn_color, structure and chain A and resi 192
color tn_color, structure and chain A and resi 193
color tn_color, structure and chain A and resi 194
color tn_color, structure and chain A and resi 195
color tn_color, structure and chain A and resi 196
color tn_color, structure and chain A and resi 200
color tn_color, structure and chain A and resi 206
color tn_color, structure and chain A and resi 207
color tn_color, structure and chain A and resi 214
color tn_color, structure and chain A and resi 216
color tn_color, structure and chain A and resi 217
color tn_color, structure and chain A and resi 218
color tn_color, structure and chain A and resi 219
color tn_color, structure and chain A and resi 220
# FN (low CSP in Binding Site): 5 residues
color fn_color, structure and chain A and resi 120
color fn_color, structure and chain A and resi 138
color fn_color, structure and chain A and resi 166
color fn_color, structure and chain A and resi 168
color fn_color, structure and chain A and resi 202
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
# TP (Sig. CSP in Union Site): 15
# FP (Sig. CSP -- Allosteric): 23
# TN (low CSP -- Allosteric): 66
# FN (low CSP in Union Site): 5
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/6r5g/6r5g_csp.pdb
