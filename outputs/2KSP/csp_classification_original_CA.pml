reinitialize
load ./outputs/2KSP/2KSP_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 10 residues
color tp_color, structure and chain A and resi 73
color tp_color, structure and chain A and resi 76
color tp_color, structure and chain A and resi 77
color tp_color, structure and chain A and resi 86
color tp_color, structure and chain A and resi 87
color tp_color, structure and chain A and resi 90
color tp_color, structure and chain A and resi 91
color tp_color, structure and chain A and resi 94
color tp_color, structure and chain A and resi 98
color tp_color, structure and chain A and resi 99
# FP (Sig. CSP -- Allosteric): 38 residues
color fp_color, structure and chain A and resi 41
color fp_color, structure and chain A and resi 42
color fp_color, structure and chain A and resi 44
color fp_color, structure and chain A and resi 47
color fp_color, structure and chain A and resi 48
color fp_color, structure and chain A and resi 52
color fp_color, structure and chain A and resi 58
color fp_color, structure and chain A and resi 60
color fp_color, structure and chain A and resi 61
color fp_color, structure and chain A and resi 63
color fp_color, structure and chain A and resi 66
color fp_color, structure and chain A and resi 68
color fp_color, structure and chain A and resi 69
color fp_color, structure and chain A and resi 70
color fp_color, structure and chain A and resi 71
color fp_color, structure and chain A and resi 72
color fp_color, structure and chain A and resi 74
color fp_color, structure and chain A and resi 75
color fp_color, structure and chain A and resi 78
color fp_color, structure and chain A and resi 79
color fp_color, structure and chain A and resi 80
color fp_color, structure and chain A and resi 81
color fp_color, structure and chain A and resi 85
color fp_color, structure and chain A and resi 88
color fp_color, structure and chain A and resi 89
color fp_color, structure and chain A and resi 93
color fp_color, structure and chain A and resi 95
color fp_color, structure and chain A and resi 100
color fp_color, structure and chain A and resi 110
color fp_color, structure and chain A and resi 113
color fp_color, structure and chain A and resi 118
color fp_color, structure and chain A and resi 119
color fp_color, structure and chain A and resi 120
color fp_color, structure and chain A and resi 122
color fp_color, structure and chain A and resi 126
color fp_color, structure and chain A and resi 129
color fp_color, structure and chain A and resi 131
color fp_color, structure and chain A and resi 137
# TN (low CSP -- Allosteric): 36 residues
color tn_color, structure and chain A and resi 43
color tn_color, structure and chain A and resi 45
color tn_color, structure and chain A and resi 49
color tn_color, structure and chain A and resi 50
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 56
color tn_color, structure and chain A and resi 57
color tn_color, structure and chain A and resi 59
color tn_color, structure and chain A and resi 64
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 67
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 101
color tn_color, structure and chain A and resi 102
color tn_color, structure and chain A and resi 103
color tn_color, structure and chain A and resi 104
color tn_color, structure and chain A and resi 106
color tn_color, structure and chain A and resi 107
color tn_color, structure and chain A and resi 108
color tn_color, structure and chain A and resi 109
color tn_color, structure and chain A and resi 111
color tn_color, structure and chain A and resi 112
color tn_color, structure and chain A and resi 114
color tn_color, structure and chain A and resi 115
color tn_color, structure and chain A and resi 116
color tn_color, structure and chain A and resi 117
color tn_color, structure and chain A and resi 121
color tn_color, structure and chain A and resi 124
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 134
color tn_color, structure and chain A and resi 135
color tn_color, structure and chain A and resi 136
color tn_color, structure and chain A and resi 139
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain A and resi 96
color fn_color, structure and chain A and resi 97
color fn_color, structure and chain A and resi 105
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
# TP (Sig. CSP in Union Site): 10
# FP (Sig. CSP -- Allosteric): 38
# TN (low CSP -- Allosteric): 36
# FN (low CSP in Union Site): 3
# Residues without CSP data: 18
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2KSP/2KSP_csp.pdb
