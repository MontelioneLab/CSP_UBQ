reinitialize
load ./outputs/2m0v/2m0v_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 7 residues
color tp_color, structure and chain A and resi 163
color tp_color, structure and chain A and resi 165
color tp_color, structure and chain A and resi 166
color tp_color, structure and chain A and resi 168
color tp_color, structure and chain A and resi 169
color tp_color, structure and chain A and resi 180
color tp_color, structure and chain A and resi 212
# FP (Sig. CSP -- Allosteric): 42 residues
color fp_color, structure and chain A and resi 156
color fp_color, structure and chain A and resi 157
color fp_color, structure and chain A and resi 160
color fp_color, structure and chain A and resi 171
color fp_color, structure and chain A and resi 173
color fp_color, structure and chain A and resi 174
color fp_color, structure and chain A and resi 176
color fp_color, structure and chain A and resi 178
color fp_color, structure and chain A and resi 179
color fp_color, structure and chain A and resi 181
color fp_color, structure and chain A and resi 182
color fp_color, structure and chain A and resi 188
color fp_color, structure and chain A and resi 189
color fp_color, structure and chain A and resi 190
color fp_color, structure and chain A and resi 191
color fp_color, structure and chain A and resi 193
color fp_color, structure and chain A and resi 194
color fp_color, structure and chain A and resi 196
color fp_color, structure and chain A and resi 198
color fp_color, structure and chain A and resi 199
color fp_color, structure and chain A and resi 200
color fp_color, structure and chain A and resi 202
color fp_color, structure and chain A and resi 203
color fp_color, structure and chain A and resi 205
color fp_color, structure and chain A and resi 207
color fp_color, structure and chain A and resi 210
color fp_color, structure and chain A and resi 211
color fp_color, structure and chain A and resi 213
color fp_color, structure and chain A and resi 214
color fp_color, structure and chain A and resi 215
color fp_color, structure and chain A and resi 216
color fp_color, structure and chain A and resi 217
color fp_color, structure and chain A and resi 218
color fp_color, structure and chain A and resi 220
color fp_color, structure and chain A and resi 221
color fp_color, structure and chain A and resi 222
color fp_color, structure and chain A and resi 223
color fp_color, structure and chain A and resi 225
color fp_color, structure and chain A and resi 226
color fp_color, structure and chain A and resi 227
color fp_color, structure and chain A and resi 237
color fp_color, structure and chain A and resi 248
# TN (low CSP -- Allosteric): 60 residues
color tn_color, structure and chain A and resi 145
color tn_color, structure and chain A and resi 147
color tn_color, structure and chain A and resi 148
color tn_color, structure and chain A and resi 149
color tn_color, structure and chain A and resi 150
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 153
color tn_color, structure and chain A and resi 154
color tn_color, structure and chain A and resi 155
color tn_color, structure and chain A and resi 158
color tn_color, structure and chain A and resi 159
color tn_color, structure and chain A and resi 162
color tn_color, structure and chain A and resi 183
color tn_color, structure and chain A and resi 185
color tn_color, structure and chain A and resi 186
color tn_color, structure and chain A and resi 192
color tn_color, structure and chain A and resi 195
color tn_color, structure and chain A and resi 197
color tn_color, structure and chain A and resi 201
color tn_color, structure and chain A and resi 204
color tn_color, structure and chain A and resi 206
color tn_color, structure and chain A and resi 208
color tn_color, structure and chain A and resi 209
color tn_color, structure and chain A and resi 224
color tn_color, structure and chain A and resi 228
color tn_color, structure and chain A and resi 229
color tn_color, structure and chain A and resi 230
color tn_color, structure and chain A and resi 231
color tn_color, structure and chain A and resi 232
color tn_color, structure and chain A and resi 233
color tn_color, structure and chain A and resi 234
color tn_color, structure and chain A and resi 235
color tn_color, structure and chain A and resi 236
color tn_color, structure and chain A and resi 238
color tn_color, structure and chain A and resi 239
color tn_color, structure and chain A and resi 240
color tn_color, structure and chain A and resi 241
color tn_color, structure and chain A and resi 242
color tn_color, structure and chain A and resi 243
color tn_color, structure and chain A and resi 244
color tn_color, structure and chain A and resi 245
color tn_color, structure and chain A and resi 247
color tn_color, structure and chain A and resi 249
color tn_color, structure and chain A and resi 250
color tn_color, structure and chain A and resi 251
color tn_color, structure and chain A and resi 253
color tn_color, structure and chain A and resi 255
color tn_color, structure and chain A and resi 257
color tn_color, structure and chain A and resi 259
color tn_color, structure and chain A and resi 260
color tn_color, structure and chain A and resi 261
color tn_color, structure and chain A and resi 262
color tn_color, structure and chain A and resi 263
color tn_color, structure and chain A and resi 264
color tn_color, structure and chain A and resi 265
color tn_color, structure and chain A and resi 266
color tn_color, structure and chain A and resi 267
color tn_color, structure and chain A and resi 268
color tn_color, structure and chain A and resi 269
color tn_color, structure and chain A and resi 270
# FN (low CSP in Binding Site): 2 residues
color fn_color, structure and chain A and resi 167
color fn_color, structure and chain A and resi 219
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
# TP (Sig. CSP in Union Site): 7
# FP (Sig. CSP -- Allosteric): 42
# TN (low CSP -- Allosteric): 60
# FN (low CSP in Union Site): 2
# Residues without CSP data: 17
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2m0v/2m0v_csp.pdb
