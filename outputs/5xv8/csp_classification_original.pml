reinitialize
load ./outputs/5xv8/5xv8_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 13 residues
color tp_color, structure and chain B and resi 18
color tp_color, structure and chain B and resi 19
color tp_color, structure and chain B and resi 53
color tp_color, structure and chain B and resi 54
color tp_color, structure and chain B and resi 55
color tp_color, structure and chain B and resi 56
color tp_color, structure and chain B and resi 60
color tp_color, structure and chain B and resi 62
color tp_color, structure and chain B and resi 64
color tp_color, structure and chain B and resi 65
color tp_color, structure and chain B and resi 66
color tp_color, structure and chain B and resi 78
color tp_color, structure and chain B and resi 94
# FP (Sig. CSP -- Allosteric): 35 residues
color fp_color, structure and chain B and resi 4
color fp_color, structure and chain B and resi 7
color fp_color, structure and chain B and resi 8
color fp_color, structure and chain B and resi 15
color fp_color, structure and chain B and resi 17
color fp_color, structure and chain B and resi 23
color fp_color, structure and chain B and resi 29
color fp_color, structure and chain B and resi 36
color fp_color, structure and chain B and resi 37
color fp_color, structure and chain B and resi 38
color fp_color, structure and chain B and resi 41
color fp_color, structure and chain B and resi 42
color fp_color, structure and chain B and resi 44
color fp_color, structure and chain B and resi 45
color fp_color, structure and chain B and resi 48
color fp_color, structure and chain B and resi 49
color fp_color, structure and chain B and resi 59
color fp_color, structure and chain B and resi 61
color fp_color, structure and chain B and resi 63
color fp_color, structure and chain B and resi 67
color fp_color, structure and chain B and resi 68
color fp_color, structure and chain B and resi 72
color fp_color, structure and chain B and resi 73
color fp_color, structure and chain B and resi 74
color fp_color, structure and chain B and resi 75
color fp_color, structure and chain B and resi 77
color fp_color, structure and chain B and resi 79
color fp_color, structure and chain B and resi 84
color fp_color, structure and chain B and resi 87
color fp_color, structure and chain B and resi 89
color fp_color, structure and chain B and resi 92
color fp_color, structure and chain B and resi 95
color fp_color, structure and chain B and resi 100
color fp_color, structure and chain B and resi 106
color fp_color, structure and chain B and resi 107
# TN (low CSP -- Allosteric): 42 residues
color tn_color, structure and chain B and resi 2
color tn_color, structure and chain B and resi 3
color tn_color, structure and chain B and resi 5
color tn_color, structure and chain B and resi 6
color tn_color, structure and chain B and resi 9
color tn_color, structure and chain B and resi 10
color tn_color, structure and chain B and resi 11
color tn_color, structure and chain B and resi 12
color tn_color, structure and chain B and resi 13
color tn_color, structure and chain B and resi 14
color tn_color, structure and chain B and resi 21
color tn_color, structure and chain B and resi 22
color tn_color, structure and chain B and resi 24
color tn_color, structure and chain B and resi 25
color tn_color, structure and chain B and resi 26
color tn_color, structure and chain B and resi 27
color tn_color, structure and chain B and resi 28
color tn_color, structure and chain B and resi 30
color tn_color, structure and chain B and resi 31
color tn_color, structure and chain B and resi 32
color tn_color, structure and chain B and resi 33
color tn_color, structure and chain B and resi 34
color tn_color, structure and chain B and resi 39
color tn_color, structure and chain B and resi 40
color tn_color, structure and chain B and resi 43
color tn_color, structure and chain B and resi 46
color tn_color, structure and chain B and resi 47
color tn_color, structure and chain B and resi 50
color tn_color, structure and chain B and resi 69
color tn_color, structure and chain B and resi 70
color tn_color, structure and chain B and resi 71
color tn_color, structure and chain B and resi 80
color tn_color, structure and chain B and resi 83
color tn_color, structure and chain B and resi 85
color tn_color, structure and chain B and resi 86
color tn_color, structure and chain B and resi 88
color tn_color, structure and chain B and resi 90
color tn_color, structure and chain B and resi 91
color tn_color, structure and chain B and resi 96
color tn_color, structure and chain B and resi 99
color tn_color, structure and chain B and resi 103
color tn_color, structure and chain B and resi 105
# FN (low CSP in Binding Site): 10 residues
color fn_color, structure and chain B and resi 16
color fn_color, structure and chain B and resi 20
color fn_color, structure and chain B and resi 51
color fn_color, structure and chain B and resi 52
color fn_color, structure and chain B and resi 76
color fn_color, structure and chain B and resi 93
color fn_color, structure and chain B and resi 97
color fn_color, structure and chain B and resi 98
color fn_color, structure and chain B and resi 102
color fn_color, structure and chain B and resi 104
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
# TP (Sig. CSP in Union Site): 13
# FP (Sig. CSP -- Allosteric): 35
# TN (low CSP -- Allosteric): 42
# FN (low CSP in Union Site): 10
# Residues without CSP data: 8
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/5xv8/5xv8_csp.pdb
