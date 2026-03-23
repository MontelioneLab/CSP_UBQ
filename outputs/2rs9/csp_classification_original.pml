reinitialize
load ./outputs/2rs9/2rs9_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# FP (Sig. CSP -- Allosteric): 33 residues
color fp_color, structure and chain B and resi 44
color fp_color, structure and chain B and resi 45
color fp_color, structure and chain B and resi 46
color fp_color, structure and chain B and resi 47
color fp_color, structure and chain B and resi 49
color fp_color, structure and chain B and resi 50
color fp_color, structure and chain B and resi 51
color fp_color, structure and chain B and resi 52
color fp_color, structure and chain B and resi 79
color fp_color, structure and chain B and resi 81
color fp_color, structure and chain B and resi 82
color fp_color, structure and chain B and resi 83
color fp_color, structure and chain B and resi 84
color fp_color, structure and chain B and resi 85
color fp_color, structure and chain B and resi 86
color fp_color, structure and chain B and resi 88
color fp_color, structure and chain B and resi 90
color fp_color, structure and chain B and resi 91
color fp_color, structure and chain B and resi 92
color fp_color, structure and chain B and resi 95
color fp_color, structure and chain B and resi 97
color fp_color, structure and chain B and resi 98
color fp_color, structure and chain B and resi 108
color fp_color, structure and chain B and resi 111
color fp_color, structure and chain B and resi 112
color fp_color, structure and chain B and resi 114
color fp_color, structure and chain B and resi 115
color fp_color, structure and chain B and resi 116
color fp_color, structure and chain B and resi 117
color fp_color, structure and chain B and resi 120
color fp_color, structure and chain B and resi 144
color fp_color, structure and chain B and resi 160
color fp_color, structure and chain B and resi 161
# TN (low CSP -- Allosteric): 60 residues
color tn_color, structure and chain B and resi 53
color tn_color, structure and chain B and resi 54
color tn_color, structure and chain B and resi 55
color tn_color, structure and chain B and resi 56
color tn_color, structure and chain B and resi 57
color tn_color, structure and chain B and resi 58
color tn_color, structure and chain B and resi 59
color tn_color, structure and chain B and resi 60
color tn_color, structure and chain B and resi 61
color tn_color, structure and chain B and resi 62
color tn_color, structure and chain B and resi 63
color tn_color, structure and chain B and resi 64
color tn_color, structure and chain B and resi 65
color tn_color, structure and chain B and resi 66
color tn_color, structure and chain B and resi 67
color tn_color, structure and chain B and resi 68
color tn_color, structure and chain B and resi 69
color tn_color, structure and chain B and resi 70
color tn_color, structure and chain B and resi 74
color tn_color, structure and chain B and resi 75
color tn_color, structure and chain B and resi 76
color tn_color, structure and chain B and resi 94
color tn_color, structure and chain B and resi 96
color tn_color, structure and chain B and resi 100
color tn_color, structure and chain B and resi 101
color tn_color, structure and chain B and resi 102
color tn_color, structure and chain B and resi 103
color tn_color, structure and chain B and resi 104
color tn_color, structure and chain B and resi 105
color tn_color, structure and chain B and resi 106
color tn_color, structure and chain B and resi 107
color tn_color, structure and chain B and resi 109
color tn_color, structure and chain B and resi 110
color tn_color, structure and chain B and resi 118
color tn_color, structure and chain B and resi 121
color tn_color, structure and chain B and resi 125
color tn_color, structure and chain B and resi 128
color tn_color, structure and chain B and resi 130
color tn_color, structure and chain B and resi 132
color tn_color, structure and chain B and resi 133
color tn_color, structure and chain B and resi 134
color tn_color, structure and chain B and resi 135
color tn_color, structure and chain B and resi 136
color tn_color, structure and chain B and resi 137
color tn_color, structure and chain B and resi 138
color tn_color, structure and chain B and resi 139
color tn_color, structure and chain B and resi 140
color tn_color, structure and chain B and resi 141
color tn_color, structure and chain B and resi 142
color tn_color, structure and chain B and resi 143
color tn_color, structure and chain B and resi 145
color tn_color, structure and chain B and resi 146
color tn_color, structure and chain B and resi 147
color tn_color, structure and chain B and resi 148
color tn_color, structure and chain B and resi 149
color tn_color, structure and chain B and resi 151
color tn_color, structure and chain B and resi 152
color tn_color, structure and chain B and resi 153
color tn_color, structure and chain B and resi 154
color tn_color, structure and chain B and resi 155
# FN (low CSP in Binding Site): 6 residues
color fn_color, structure and chain B and resi 77
color fn_color, structure and chain B and resi 119
color fn_color, structure and chain B and resi 122
color fn_color, structure and chain B and resi 123
color fn_color, structure and chain B and resi 124
color fn_color, structure and chain B and resi 127
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
# TP (Sig. CSP in Union Site): 0
# FP (Sig. CSP -- Allosteric): 33
# TN (low CSP -- Allosteric): 60
# FN (low CSP in Union Site): 6
# Residues without CSP data: 22
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2rs9/2rs9_csp.pdb
