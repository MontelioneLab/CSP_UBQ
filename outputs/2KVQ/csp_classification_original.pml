reinitialize
load ./outputs/2KVQ/2KVQ_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain E
set cartoon_tube_radius, 0.45, structure and chain E
color gray, structure and chain G
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 3 residues
color tp_color, structure and chain G and resi 164
color tp_color, structure and chain G and resi 167
color tp_color, structure and chain G and resi 171
# FP (Sig. CSP -- Allosteric): 16 residues
color fp_color, structure and chain G and resi 129
color fp_color, structure and chain G and resi 138
color fp_color, structure and chain G and resi 139
color fp_color, structure and chain G and resi 143
color fp_color, structure and chain G and resi 144
color fp_color, structure and chain G and resi 145
color fp_color, structure and chain G and resi 146
color fp_color, structure and chain G and resi 156
color fp_color, structure and chain G and resi 159
color fp_color, structure and chain G and resi 160
color fp_color, structure and chain G and resi 161
color fp_color, structure and chain G and resi 162
color fp_color, structure and chain G and resi 163
color fp_color, structure and chain G and resi 173
color fp_color, structure and chain G and resi 174
color fp_color, structure and chain G and resi 178
# TN (low CSP -- Allosteric): 27 residues
color tn_color, structure and chain G and resi 127
color tn_color, structure and chain G and resi 128
color tn_color, structure and chain G and resi 131
color tn_color, structure and chain G and resi 132
color tn_color, structure and chain G and resi 133
color tn_color, structure and chain G and resi 134
color tn_color, structure and chain G and resi 135
color tn_color, structure and chain G and resi 136
color tn_color, structure and chain G and resi 137
color tn_color, structure and chain G and resi 142
color tn_color, structure and chain G and resi 147
color tn_color, structure and chain G and resi 148
color tn_color, structure and chain G and resi 149
color tn_color, structure and chain G and resi 150
color tn_color, structure and chain G and resi 151
color tn_color, structure and chain G and resi 152
color tn_color, structure and chain G and resi 153
color tn_color, structure and chain G and resi 154
color tn_color, structure and chain G and resi 155
color tn_color, structure and chain G and resi 157
color tn_color, structure and chain G and resi 158
color tn_color, structure and chain G and resi 175
color tn_color, structure and chain G and resi 176
color tn_color, structure and chain G and resi 177
color tn_color, structure and chain G and resi 179
color tn_color, structure and chain G and resi 180
color tn_color, structure and chain G and resi 181
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
# Receptor chain: G
# Ligand chain: E
# TP (Sig. CSP in Union Site): 3
# FP (Sig. CSP -- Allosteric): 16
# TN (low CSP -- Allosteric): 27
# FN (low CSP in Union Site): 0
# Residues without CSP data: 13
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2KVQ/2KVQ_csp.pdb
