reinitialize
load ./outputs/2JXC/2JXC_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 16 residues
color tp_color, structure and chain A and resi 152
color tp_color, structure and chain A and resi 155
color tp_color, structure and chain A and resi 158
color tp_color, structure and chain A and resi 159
color tp_color, structure and chain A and resi 162
color tp_color, structure and chain A and resi 163
color tp_color, structure and chain A and resi 165
color tp_color, structure and chain A and resi 166
color tp_color, structure and chain A and resi 167
color tp_color, structure and chain A and resi 169
color tp_color, structure and chain A and resi 170
color tp_color, structure and chain A and resi 173
color tp_color, structure and chain A and resi 176
color tp_color, structure and chain A and resi 192
color tp_color, structure and chain A and resi 193
color tp_color, structure and chain A and resi 197
# FP (Sig. CSP -- Allosteric): 14 residues
color fp_color, structure and chain A and resi 132
color fp_color, structure and chain A and resi 144
color fp_color, structure and chain A and resi 149
color fp_color, structure and chain A and resi 151
color fp_color, structure and chain A and resi 154
color fp_color, structure and chain A and resi 160
color fp_color, structure and chain A and resi 164
color fp_color, structure and chain A and resi 171
color fp_color, structure and chain A and resi 190
color fp_color, structure and chain A and resi 194
color fp_color, structure and chain A and resi 195
color fp_color, structure and chain A and resi 213
color fp_color, structure and chain A and resi 214
color fp_color, structure and chain A and resi 215
# TN (low CSP -- Allosteric): 42 residues
color tn_color, structure and chain A and resi 122
color tn_color, structure and chain A and resi 123
color tn_color, structure and chain A and resi 124
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 127
color tn_color, structure and chain A and resi 128
color tn_color, structure and chain A and resi 129
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 131
color tn_color, structure and chain A and resi 133
color tn_color, structure and chain A and resi 134
color tn_color, structure and chain A and resi 135
color tn_color, structure and chain A and resi 136
color tn_color, structure and chain A and resi 137
color tn_color, structure and chain A and resi 138
color tn_color, structure and chain A and resi 140
color tn_color, structure and chain A and resi 142
color tn_color, structure and chain A and resi 143
color tn_color, structure and chain A and resi 145
color tn_color, structure and chain A and resi 146
color tn_color, structure and chain A and resi 147
color tn_color, structure and chain A and resi 150
color tn_color, structure and chain A and resi 168
color tn_color, structure and chain A and resi 172
color tn_color, structure and chain A and resi 174
color tn_color, structure and chain A and resi 180
color tn_color, structure and chain A and resi 181
color tn_color, structure and chain A and resi 182
color tn_color, structure and chain A and resi 183
color tn_color, structure and chain A and resi 185
color tn_color, structure and chain A and resi 186
color tn_color, structure and chain A and resi 187
color tn_color, structure and chain A and resi 188
color tn_color, structure and chain A and resi 191
color tn_color, structure and chain A and resi 199
color tn_color, structure and chain A and resi 201
color tn_color, structure and chain A and resi 203
color tn_color, structure and chain A and resi 204
color tn_color, structure and chain A and resi 205
color tn_color, structure and chain A and resi 208
color tn_color, structure and chain A and resi 209
color tn_color, structure and chain A and resi 210
# FN (low CSP in Binding Site): 10 residues
color fn_color, structure and chain A and resi 139
color fn_color, structure and chain A and resi 148
color fn_color, structure and chain A and resi 175
color fn_color, structure and chain A and resi 177
color fn_color, structure and chain A and resi 178
color fn_color, structure and chain A and resi 179
color fn_color, structure and chain A and resi 184
color fn_color, structure and chain A and resi 189
color fn_color, structure and chain A and resi 196
color fn_color, structure and chain A and resi 198
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
# TP (Sig. CSP in Union Site): 16
# FP (Sig. CSP -- Allosteric): 14
# TN (low CSP -- Allosteric): 42
# FN (low CSP in Union Site): 10
# Residues without CSP data: 13
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2JXC/2JXC_csp.pdb
