reinitialize
load ./outputs/6U4N/6U4N_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 2 residues
color tp_color, structure and chain A and resi 42
color tp_color, structure and chain A and resi 46
# FP (Sig. CSP -- Allosteric): 28 residues
color fp_color, structure and chain A and resi 7
color fp_color, structure and chain A and resi 10
color fp_color, structure and chain A and resi 14
color fp_color, structure and chain A and resi 15
color fp_color, structure and chain A and resi 21
color fp_color, structure and chain A and resi 23
color fp_color, structure and chain A and resi 24
color fp_color, structure and chain A and resi 25
color fp_color, structure and chain A and resi 28
color fp_color, structure and chain A and resi 30
color fp_color, structure and chain A and resi 31
color fp_color, structure and chain A and resi 32
color fp_color, structure and chain A and resi 33
color fp_color, structure and chain A and resi 40
color fp_color, structure and chain A and resi 41
color fp_color, structure and chain A and resi 44
color fp_color, structure and chain A and resi 48
color fp_color, structure and chain A and resi 51
color fp_color, structure and chain A and resi 60
color fp_color, structure and chain A and resi 62
color fp_color, structure and chain A and resi 63
color fp_color, structure and chain A and resi 76
color fp_color, structure and chain A and resi 77
color fp_color, structure and chain A and resi 82
color fp_color, structure and chain A and resi 84
color fp_color, structure and chain A and resi 86
color fp_color, structure and chain A and resi 89
color fp_color, structure and chain A and resi 93
# TN (low CSP -- Allosteric): 43 residues
color tn_color, structure and chain A and resi 1
color tn_color, structure and chain A and resi 2
color tn_color, structure and chain A and resi 3
color tn_color, structure and chain A and resi 4
color tn_color, structure and chain A and resi 5
color tn_color, structure and chain A and resi 6
color tn_color, structure and chain A and resi 8
color tn_color, structure and chain A and resi 11
color tn_color, structure and chain A and resi 12
color tn_color, structure and chain A and resi 13
color tn_color, structure and chain A and resi 16
color tn_color, structure and chain A and resi 17
color tn_color, structure and chain A and resi 18
color tn_color, structure and chain A and resi 19
color tn_color, structure and chain A and resi 22
color tn_color, structure and chain A and resi 29
color tn_color, structure and chain A and resi 34
color tn_color, structure and chain A and resi 35
color tn_color, structure and chain A and resi 36
color tn_color, structure and chain A and resi 37
color tn_color, structure and chain A and resi 38
color tn_color, structure and chain A and resi 39
color tn_color, structure and chain A and resi 61
color tn_color, structure and chain A and resi 64
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 66
color tn_color, structure and chain A and resi 67
color tn_color, structure and chain A and resi 68
color tn_color, structure and chain A and resi 69
color tn_color, structure and chain A and resi 70
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 72
color tn_color, structure and chain A and resi 74
color tn_color, structure and chain A and resi 78
color tn_color, structure and chain A and resi 79
color tn_color, structure and chain A and resi 80
color tn_color, structure and chain A and resi 81
color tn_color, structure and chain A and resi 83
color tn_color, structure and chain A and resi 85
color tn_color, structure and chain A and resi 88
color tn_color, structure and chain A and resi 90
color tn_color, structure and chain A and resi 91
color tn_color, structure and chain A and resi 92
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain A and resi 57
color fn_color, structure and chain A and resi 59
color fn_color, structure and chain A and resi 73
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
# TP (Sig. CSP in Union Site): 2
# FP (Sig. CSP -- Allosteric): 28
# TN (low CSP -- Allosteric): 43
# FN (low CSP in Union Site): 3
# Residues without CSP data: 18
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/6U4N/6U4N_csp.pdb
