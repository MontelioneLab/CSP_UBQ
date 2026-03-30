reinitialize
load ./outputs/2l75/2l75_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain A and resi 91
color tp_color, structure and chain A and resi 98
color tp_color, structure and chain A and resi 102
color tp_color, structure and chain A and resi 103
color tp_color, structure and chain A and resi 104
color tp_color, structure and chain A and resi 124
color tp_color, structure and chain A and resi 129
color tp_color, structure and chain A and resi 134
# FP (Sig. CSP -- Allosteric): 10 residues
color fp_color, structure and chain A and resi 110
color fp_color, structure and chain A and resi 111
color fp_color, structure and chain A and resi 114
color fp_color, structure and chain A and resi 115
color fp_color, structure and chain A and resi 117
color fp_color, structure and chain A and resi 118
color fp_color, structure and chain A and resi 121
color fp_color, structure and chain A and resi 123
color fp_color, structure and chain A and resi 130
color fp_color, structure and chain A and resi 139
# TN (low CSP -- Allosteric): 12 residues
color tn_color, structure and chain A and resi 87
color tn_color, structure and chain A and resi 90
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 94
color tn_color, structure and chain A and resi 95
color tn_color, structure and chain A and resi 107
color tn_color, structure and chain A and resi 112
color tn_color, structure and chain A and resi 126
color tn_color, structure and chain A and resi 133
color tn_color, structure and chain A and resi 140
color tn_color, structure and chain A and resi 141
color tn_color, structure and chain A and resi 142
# FN (low CSP in Binding Site): 15 residues
color fn_color, structure and chain A and resi 89
color fn_color, structure and chain A and resi 93
color fn_color, structure and chain A and resi 96
color fn_color, structure and chain A and resi 97
color fn_color, structure and chain A and resi 99
color fn_color, structure and chain A and resi 100
color fn_color, structure and chain A and resi 101
color fn_color, structure and chain A and resi 105
color fn_color, structure and chain A and resi 106
color fn_color, structure and chain A and resi 108
color fn_color, structure and chain A and resi 113
color fn_color, structure and chain A and resi 116
color fn_color, structure and chain A and resi 127
color fn_color, structure and chain A and resi 128
color fn_color, structure and chain A and resi 131
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 10
# TN (low CSP -- Allosteric): 12
# FN (low CSP in Union Site): 15
# Residues without CSP data: 11
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2l75/2l75_csp.pdb
