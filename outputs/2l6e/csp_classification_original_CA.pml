reinitialize
load ./outputs/2l6e/2l6e_csp.pdb, structure
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
color tp_color, structure and chain A and resi 169
color tp_color, structure and chain A and resi 172
color tp_color, structure and chain A and resi 173
color tp_color, structure and chain A and resi 179
color tp_color, structure and chain A and resi 183
color tp_color, structure and chain A and resi 186
color tp_color, structure and chain A and resi 187
color tp_color, structure and chain A and resi 190
color tp_color, structure and chain A and resi 209
color tp_color, structure and chain A and resi 211
color tp_color, structure and chain A and resi 212
color tp_color, structure and chain A and resi 214
color tp_color, structure and chain A and resi 215
# FP (Sig. CSP -- Allosteric): 21 residues
color fp_color, structure and chain A and resi 145
color fp_color, structure and chain A and resi 148
color fp_color, structure and chain A and resi 150
color fp_color, structure and chain A and resi 151
color fp_color, structure and chain A and resi 167
color fp_color, structure and chain A and resi 168
color fp_color, structure and chain A and resi 170
color fp_color, structure and chain A and resi 171
color fp_color, structure and chain A and resi 174
color fp_color, structure and chain A and resi 175
color fp_color, structure and chain A and resi 180
color fp_color, structure and chain A and resi 181
color fp_color, structure and chain A and resi 182
color fp_color, structure and chain A and resi 184
color fp_color, structure and chain A and resi 185
color fp_color, structure and chain A and resi 188
color fp_color, structure and chain A and resi 189
color fp_color, structure and chain A and resi 191
color fp_color, structure and chain A and resi 193
color fp_color, structure and chain A and resi 206
color fp_color, structure and chain A and resi 227
# TN (low CSP -- Allosteric): 41 residues
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 153
color tn_color, structure and chain A and resi 154
color tn_color, structure and chain A and resi 155
color tn_color, structure and chain A and resi 156
color tn_color, structure and chain A and resi 158
color tn_color, structure and chain A and resi 159
color tn_color, structure and chain A and resi 161
color tn_color, structure and chain A and resi 162
color tn_color, structure and chain A and resi 163
color tn_color, structure and chain A and resi 164
color tn_color, structure and chain A and resi 176
color tn_color, structure and chain A and resi 177
color tn_color, structure and chain A and resi 178
color tn_color, structure and chain A and resi 192
color tn_color, structure and chain A and resi 194
color tn_color, structure and chain A and resi 195
color tn_color, structure and chain A and resi 197
color tn_color, structure and chain A and resi 198
color tn_color, structure and chain A and resi 199
color tn_color, structure and chain A and resi 200
color tn_color, structure and chain A and resi 201
color tn_color, structure and chain A and resi 202
color tn_color, structure and chain A and resi 203
color tn_color, structure and chain A and resi 204
color tn_color, structure and chain A and resi 205
color tn_color, structure and chain A and resi 208
color tn_color, structure and chain A and resi 213
color tn_color, structure and chain A and resi 216
color tn_color, structure and chain A and resi 217
color tn_color, structure and chain A and resi 218
color tn_color, structure and chain A and resi 219
color tn_color, structure and chain A and resi 220
color tn_color, structure and chain A and resi 221
color tn_color, structure and chain A and resi 222
color tn_color, structure and chain A and resi 223
color tn_color, structure and chain A and resi 225
color tn_color, structure and chain A and resi 226
color tn_color, structure and chain A and resi 228
color tn_color, structure and chain A and resi 229
color tn_color, structure and chain A and resi 230
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain A and resi 165
color fn_color, structure and chain A and resi 166
color fn_color, structure and chain A and resi 210
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
# FP (Sig. CSP -- Allosteric): 21
# TN (low CSP -- Allosteric): 41
# FN (low CSP in Union Site): 3
# Residues without CSP data: 16
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2l6e/2l6e_csp.pdb
