reinitialize
load ./outputs/2RMK/2RMK_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 5 residues
color tp_color, structure and chain A and resi 63
color tp_color, structure and chain A and resi 116
color tp_color, structure and chain A and resi 118
color tp_color, structure and chain A and resi 185
color tp_color, structure and chain A and resi 188
# FP (Sig. CSP -- Allosteric): 41 residues
color fp_color, structure and chain A and resi 12
color fp_color, structure and chain A and resi 19
color fp_color, structure and chain A and resi 20
color fp_color, structure and chain A and resi 21
color fp_color, structure and chain A and resi 27
color fp_color, structure and chain A and resi 30
color fp_color, structure and chain A and resi 42
color fp_color, structure and chain A and resi 43
color fp_color, structure and chain A and resi 51
color fp_color, structure and chain A and resi 56
color fp_color, structure and chain A and resi 79
color fp_color, structure and chain A and resi 82
color fp_color, structure and chain A and resi 89
color fp_color, structure and chain A and resi 96
color fp_color, structure and chain A and resi 101
color fp_color, structure and chain A and resi 102
color fp_color, structure and chain A and resi 112
color fp_color, structure and chain A and resi 117
color fp_color, structure and chain A and resi 119
color fp_color, structure and chain A and resi 124
color fp_color, structure and chain A and resi 125
color fp_color, structure and chain A and resi 126
color fp_color, structure and chain A and resi 129
color fp_color, structure and chain A and resi 138
color fp_color, structure and chain A and resi 139
color fp_color, structure and chain A and resi 142
color fp_color, structure and chain A and resi 143
color fp_color, structure and chain A and resi 144
color fp_color, structure and chain A and resi 145
color fp_color, structure and chain A and resi 146
color fp_color, structure and chain A and resi 149
color fp_color, structure and chain A and resi 150
color fp_color, structure and chain A and resi 153
color fp_color, structure and chain A and resi 155
color fp_color, structure and chain A and resi 156
color fp_color, structure and chain A and resi 165
color fp_color, structure and chain A and resi 170
color fp_color, structure and chain A and resi 172
color fp_color, structure and chain A and resi 177
color fp_color, structure and chain A and resi 182
color fp_color, structure and chain A and resi 183
# TN (low CSP -- Allosteric): 60 residues
color tn_color, structure and chain A and resi 10
color tn_color, structure and chain A and resi 11
color tn_color, structure and chain A and resi 24
color tn_color, structure and chain A and resi 25
color tn_color, structure and chain A and resi 26
color tn_color, structure and chain A and resi 45
color tn_color, structure and chain A and resi 46
color tn_color, structure and chain A and resi 47
color tn_color, structure and chain A and resi 48
color tn_color, structure and chain A and resi 49
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 72
color tn_color, structure and chain A and resi 78
color tn_color, structure and chain A and resi 84
color tn_color, structure and chain A and resi 85
color tn_color, structure and chain A and resi 86
color tn_color, structure and chain A and resi 90
color tn_color, structure and chain A and resi 91
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 93
color tn_color, structure and chain A and resi 94
color tn_color, structure and chain A and resi 95
color tn_color, structure and chain A and resi 98
color tn_color, structure and chain A and resi 110
color tn_color, structure and chain A and resi 111
color tn_color, structure and chain A and resi 113
color tn_color, structure and chain A and resi 114
color tn_color, structure and chain A and resi 115
color tn_color, structure and chain A and resi 120
color tn_color, structure and chain A and resi 121
color tn_color, structure and chain A and resi 122
color tn_color, structure and chain A and resi 123
color tn_color, structure and chain A and resi 127
color tn_color, structure and chain A and resi 128
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 131
color tn_color, structure and chain A and resi 132
color tn_color, structure and chain A and resi 133
color tn_color, structure and chain A and resi 134
color tn_color, structure and chain A and resi 135
color tn_color, structure and chain A and resi 137
color tn_color, structure and chain A and resi 141
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 161
color tn_color, structure and chain A and resi 162
color tn_color, structure and chain A and resi 163
color tn_color, structure and chain A and resi 164
color tn_color, structure and chain A and resi 166
color tn_color, structure and chain A and resi 167
color tn_color, structure and chain A and resi 168
color tn_color, structure and chain A and resi 169
color tn_color, structure and chain A and resi 173
color tn_color, structure and chain A and resi 174
color tn_color, structure and chain A and resi 175
color tn_color, structure and chain A and resi 176
color tn_color, structure and chain A and resi 184
# FN (low CSP in Binding Site): 7 residues
color fn_color, structure and chain A and resi 15
color fn_color, structure and chain A and resi 18
color fn_color, structure and chain A and resi 28
color fn_color, structure and chain A and resi 70
color fn_color, structure and chain A and resi 158
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
# TP (Sig. CSP in Union Site): 5
# FP (Sig. CSP -- Allosteric): 41
# TN (low CSP -- Allosteric): 60
# FN (low CSP in Union Site): 7
# Residues without CSP data: 79
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2RMK/2RMK_csp.pdb
