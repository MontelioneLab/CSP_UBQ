reinitialize
load ./outputs/1l8c/1l8c_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 23 residues
color tp_color, structure and chain A and resi 9
color tp_color, structure and chain A and resi 14
color tp_color, structure and chain A and resi 17
color tp_color, structure and chain A and resi 20
color tp_color, structure and chain A and resi 22
color tp_color, structure and chain A and resi 49
color tp_color, structure and chain A and resi 50
color tp_color, structure and chain A and resi 51
color tp_color, structure and chain A and resi 53
color tp_color, structure and chain A and resi 57
color tp_color, structure and chain A and resi 63
color tp_color, structure and chain A and resi 65
color tp_color, structure and chain A and resi 66
color tp_color, structure and chain A and resi 68
color tp_color, structure and chain A and resi 69
color tp_color, structure and chain A and resi 71
color tp_color, structure and chain A and resi 72
color tp_color, structure and chain A and resi 73
color tp_color, structure and chain A and resi 74
color tp_color, structure and chain A and resi 75
color tp_color, structure and chain A and resi 76
color tp_color, structure and chain A and resi 77
color tp_color, structure and chain A and resi 84
# FP (Sig. CSP -- Allosteric): 13 residues
color fp_color, structure and chain A and resi 4
color fp_color, structure and chain A and resi 6
color fp_color, structure and chain A and resi 16
color fp_color, structure and chain A and resi 19
color fp_color, structure and chain A and resi 28
color fp_color, structure and chain A and resi 29
color fp_color, structure and chain A and resi 30
color fp_color, structure and chain A and resi 31
color fp_color, structure and chain A and resi 48
color fp_color, structure and chain A and resi 60
color fp_color, structure and chain A and resi 67
color fp_color, structure and chain A and resi 70
color fp_color, structure and chain A and resi 89
# TN (low CSP -- Allosteric): 24 residues
color tn_color, structure and chain A and resi 10
color tn_color, structure and chain A and resi 13
color tn_color, structure and chain A and resi 25
color tn_color, structure and chain A and resi 26
color tn_color, structure and chain A and resi 32
color tn_color, structure and chain A and resi 33
color tn_color, structure and chain A and resi 34
color tn_color, structure and chain A and resi 42
color tn_color, structure and chain A and resi 45
color tn_color, structure and chain A and resi 46
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 56
color tn_color, structure and chain A and resi 58
color tn_color, structure and chain A and resi 61
color tn_color, structure and chain A and resi 78
color tn_color, structure and chain A and resi 81
color tn_color, structure and chain A and resi 86
color tn_color, structure and chain A and resi 90
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 93
color tn_color, structure and chain A and resi 94
color tn_color, structure and chain A and resi 95
# FN (low CSP in Binding Site): 16 residues
color fn_color, structure and chain A and resi 5
color fn_color, structure and chain A and resi 7
color fn_color, structure and chain A and resi 8
color fn_color, structure and chain A and resi 11
color fn_color, structure and chain A and resi 24
color fn_color, structure and chain A and resi 35
color fn_color, structure and chain A and resi 43
color fn_color, structure and chain A and resi 44
color fn_color, structure and chain A and resi 47
color fn_color, structure and chain A and resi 59
color fn_color, structure and chain A and resi 62
color fn_color, structure and chain A and resi 64
color fn_color, structure and chain A and resi 82
color fn_color, structure and chain A and resi 85
color fn_color, structure and chain A and resi 88
color fn_color, structure and chain A and resi 91
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
# TP (Sig. CSP in Union Site): 23
# FP (Sig. CSP -- Allosteric): 13
# TN (low CSP -- Allosteric): 24
# FN (low CSP in Union Site): 16
# Residues without CSP data: 19
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/1l8c/1l8c_csp.pdb
