reinitialize
load ./outputs/2M0G_1/2M0G_csp.pdb, structure
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
color tp_color, structure and chain A and resi 27
color tp_color, structure and chain A and resi 33
color tp_color, structure and chain A and resi 38
color tp_color, structure and chain A and resi 40
color tp_color, structure and chain A and resi 49
color tp_color, structure and chain A and resi 53
color tp_color, structure and chain A and resi 57
# FP (Sig. CSP -- Allosteric): 26 residues
color fp_color, structure and chain A and resi 5
color fp_color, structure and chain A and resi 32
color fp_color, structure and chain A and resi 35
color fp_color, structure and chain A and resi 36
color fp_color, structure and chain A and resi 44
color fp_color, structure and chain A and resi 45
color fp_color, structure and chain A and resi 46
color fp_color, structure and chain A and resi 50
color fp_color, structure and chain A and resi 51
color fp_color, structure and chain A and resi 52
color fp_color, structure and chain A and resi 54
color fp_color, structure and chain A and resi 55
color fp_color, structure and chain A and resi 56
color fp_color, structure and chain A and resi 59
color fp_color, structure and chain A and resi 67
color fp_color, structure and chain A and resi 93
color fp_color, structure and chain A and resi 94
color fp_color, structure and chain A and resi 97
color fp_color, structure and chain A and resi 101
color fp_color, structure and chain A and resi 103
color fp_color, structure and chain A and resi 105
color fp_color, structure and chain A and resi 108
color fp_color, structure and chain A and resi 109
color fp_color, structure and chain A and resi 112
color fp_color, structure and chain A and resi 114
color fp_color, structure and chain A and resi 115
# TN (low CSP -- Allosteric): 49 residues
color tn_color, structure and chain A and resi 31
color tn_color, structure and chain A and resi 43
color tn_color, structure and chain A and resi 47
color tn_color, structure and chain A and resi 48
color tn_color, structure and chain A and resi 58
color tn_color, structure and chain A and resi 62
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 66
color tn_color, structure and chain A and resi 68
color tn_color, structure and chain A and resi 69
color tn_color, structure and chain A and resi 70
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 72
color tn_color, structure and chain A and resi 75
color tn_color, structure and chain A and resi 77
color tn_color, structure and chain A and resi 78
color tn_color, structure and chain A and resi 79
color tn_color, structure and chain A and resi 80
color tn_color, structure and chain A and resi 82
color tn_color, structure and chain A and resi 84
color tn_color, structure and chain A and resi 86
color tn_color, structure and chain A and resi 87
color tn_color, structure and chain A and resi 88
color tn_color, structure and chain A and resi 89
color tn_color, structure and chain A and resi 91
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 95
color tn_color, structure and chain A and resi 96
color tn_color, structure and chain A and resi 98
color tn_color, structure and chain A and resi 102
color tn_color, structure and chain A and resi 104
color tn_color, structure and chain A and resi 106
color tn_color, structure and chain A and resi 107
color tn_color, structure and chain A and resi 110
color tn_color, structure and chain A and resi 111
color tn_color, structure and chain A and resi 113
color tn_color, structure and chain A and resi 116
color tn_color, structure and chain A and resi 117
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 119
color tn_color, structure and chain A and resi 120
color tn_color, structure and chain A and resi 122
color tn_color, structure and chain A and resi 123
color tn_color, structure and chain A and resi 124
color tn_color, structure and chain A and resi 127
color tn_color, structure and chain A and resi 128
color tn_color, structure and chain A and resi 129
color tn_color, structure and chain A and resi 130
color tn_color, structure and chain A and resi 133
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
# FP (Sig. CSP -- Allosteric): 26
# TN (low CSP -- Allosteric): 49
# FN (low CSP in Union Site): 0
# Residues without CSP data: 63
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2M0G_1/2M0G_csp.pdb
