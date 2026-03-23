reinitialize
load ./outputs/2YKA/2YKA_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 14 residues
color tp_color, structure and chain A and resi 83
color tp_color, structure and chain A and resi 85
color tp_color, structure and chain A and resi 86
color tp_color, structure and chain A and resi 90
color tp_color, structure and chain A and resi 91
color tp_color, structure and chain A and resi 94
color tp_color, structure and chain A and resi 97
color tp_color, structure and chain A and resi 134
color tp_color, structure and chain A and resi 135
color tp_color, structure and chain A and resi 137
color tp_color, structure and chain A and resi 138
color tp_color, structure and chain A and resi 140
color tp_color, structure and chain A and resi 141
color tp_color, structure and chain A and resi 145
# FP (Sig. CSP -- Allosteric): 36 residues
color fp_color, structure and chain A and resi 78
color fp_color, structure and chain A and resi 79
color fp_color, structure and chain A and resi 80
color fp_color, structure and chain A and resi 81
color fp_color, structure and chain A and resi 84
color fp_color, structure and chain A and resi 88
color fp_color, structure and chain A and resi 89
color fp_color, structure and chain A and resi 92
color fp_color, structure and chain A and resi 93
color fp_color, structure and chain A and resi 95
color fp_color, structure and chain A and resi 96
color fp_color, structure and chain A and resi 98
color fp_color, structure and chain A and resi 101
color fp_color, structure and chain A and resi 102
color fp_color, structure and chain A and resi 104
color fp_color, structure and chain A and resi 105
color fp_color, structure and chain A and resi 106
color fp_color, structure and chain A and resi 107
color fp_color, structure and chain A and resi 114
color fp_color, structure and chain A and resi 115
color fp_color, structure and chain A and resi 117
color fp_color, structure and chain A and resi 119
color fp_color, structure and chain A and resi 120
color fp_color, structure and chain A and resi 124
color fp_color, structure and chain A and resi 128
color fp_color, structure and chain A and resi 130
color fp_color, structure and chain A and resi 131
color fp_color, structure and chain A and resi 132
color fp_color, structure and chain A and resi 133
color fp_color, structure and chain A and resi 136
color fp_color, structure and chain A and resi 142
color fp_color, structure and chain A and resi 143
color fp_color, structure and chain A and resi 146
color fp_color, structure and chain A and resi 147
color fp_color, structure and chain A and resi 148
color fp_color, structure and chain A and resi 163
# TN (low CSP -- Allosteric): 60 residues
color tn_color, structure and chain A and resi 43
color tn_color, structure and chain A and resi 44
color tn_color, structure and chain A and resi 47
color tn_color, structure and chain A and resi 48
color tn_color, structure and chain A and resi 49
color tn_color, structure and chain A and resi 50
color tn_color, structure and chain A and resi 51
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 56
color tn_color, structure and chain A and resi 57
color tn_color, structure and chain A and resi 58
color tn_color, structure and chain A and resi 59
color tn_color, structure and chain A and resi 60
color tn_color, structure and chain A and resi 61
color tn_color, structure and chain A and resi 62
color tn_color, structure and chain A and resi 63
color tn_color, structure and chain A and resi 64
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 66
color tn_color, structure and chain A and resi 67
color tn_color, structure and chain A and resi 68
color tn_color, structure and chain A and resi 69
color tn_color, structure and chain A and resi 70
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 72
color tn_color, structure and chain A and resi 73
color tn_color, structure and chain A and resi 74
color tn_color, structure and chain A and resi 75
color tn_color, structure and chain A and resi 76
color tn_color, structure and chain A and resi 77
color tn_color, structure and chain A and resi 82
color tn_color, structure and chain A and resi 99
color tn_color, structure and chain A and resi 100
color tn_color, structure and chain A and resi 103
color tn_color, structure and chain A and resi 108
color tn_color, structure and chain A and resi 109
color tn_color, structure and chain A and resi 110
color tn_color, structure and chain A and resi 111
color tn_color, structure and chain A and resi 112
color tn_color, structure and chain A and resi 113
color tn_color, structure and chain A and resi 116
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 121
color tn_color, structure and chain A and resi 122
color tn_color, structure and chain A and resi 123
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 126
color tn_color, structure and chain A and resi 127
color tn_color, structure and chain A and resi 129
color tn_color, structure and chain A and resi 149
color tn_color, structure and chain A and resi 150
color tn_color, structure and chain A and resi 151
color tn_color, structure and chain A and resi 152
color tn_color, structure and chain A and resi 153
color tn_color, structure and chain A and resi 154
color tn_color, structure and chain A and resi 155
color tn_color, structure and chain A and resi 156
color tn_color, structure and chain A and resi 157
# FN (low CSP in Binding Site): 1 residues
color fn_color, structure and chain A and resi 87
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
# TP (Sig. CSP in Union Site): 14
# FP (Sig. CSP -- Allosteric): 36
# TN (low CSP -- Allosteric): 60
# FN (low CSP in Union Site): 1
# Residues without CSP data: 13
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2YKA/2YKA_csp.pdb
