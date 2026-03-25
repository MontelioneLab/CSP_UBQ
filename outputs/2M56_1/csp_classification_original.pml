reinitialize
load ./outputs/2M56_1/2M56_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 2 residues
color tp_color, structure and chain B and resi 28
color tp_color, structure and chain B and resi 66
# FP (Sig. CSP -- Allosteric): 28 residues
color fp_color, structure and chain B and resi 2
color fp_color, structure and chain B and resi 3
color fp_color, structure and chain B and resi 5
color fp_color, structure and chain B and resi 6
color fp_color, structure and chain B and resi 10
color fp_color, structure and chain B and resi 11
color fp_color, structure and chain B and resi 15
color fp_color, structure and chain B and resi 19
color fp_color, structure and chain B and resi 26
color fp_color, structure and chain B and resi 27
color fp_color, structure and chain B and resi 29
color fp_color, structure and chain B and resi 30
color fp_color, structure and chain B and resi 32
color fp_color, structure and chain B and resi 33
color fp_color, structure and chain B and resi 52
color fp_color, structure and chain B and resi 53
color fp_color, structure and chain B and resi 57
color fp_color, structure and chain B and resi 59
color fp_color, structure and chain B and resi 60
color fp_color, structure and chain B and resi 70
color fp_color, structure and chain B and resi 72
color fp_color, structure and chain B and resi 81
color fp_color, structure and chain B and resi 82
color fp_color, structure and chain B and resi 88
color fp_color, structure and chain B and resi 91
color fp_color, structure and chain B and resi 94
color fp_color, structure and chain B and resi 97
color fp_color, structure and chain B and resi 99
# TN (low CSP -- Allosteric): 40 residues
color tn_color, structure and chain B and resi 4
color tn_color, structure and chain B and resi 7
color tn_color, structure and chain B and resi 9
color tn_color, structure and chain B and resi 12
color tn_color, structure and chain B and resi 13
color tn_color, structure and chain B and resi 14
color tn_color, structure and chain B and resi 18
color tn_color, structure and chain B and resi 20
color tn_color, structure and chain B and resi 21
color tn_color, structure and chain B and resi 22
color tn_color, structure and chain B and resi 23
color tn_color, structure and chain B and resi 31
color tn_color, structure and chain B and resi 35
color tn_color, structure and chain B and resi 51
color tn_color, structure and chain B and resi 54
color tn_color, structure and chain B and resi 55
color tn_color, structure and chain B and resi 56
color tn_color, structure and chain B and resi 58
color tn_color, structure and chain B and resi 62
color tn_color, structure and chain B and resi 63
color tn_color, structure and chain B and resi 64
color tn_color, structure and chain B and resi 65
color tn_color, structure and chain B and resi 67
color tn_color, structure and chain B and resi 68
color tn_color, structure and chain B and resi 69
color tn_color, structure and chain B and resi 74
color tn_color, structure and chain B and resi 75
color tn_color, structure and chain B and resi 76
color tn_color, structure and chain B and resi 77
color tn_color, structure and chain B and resi 78
color tn_color, structure and chain B and resi 79
color tn_color, structure and chain B and resi 83
color tn_color, structure and chain B and resi 89
color tn_color, structure and chain B and resi 90
color tn_color, structure and chain B and resi 93
color tn_color, structure and chain B and resi 95
color tn_color, structure and chain B and resi 96
color tn_color, structure and chain B and resi 100
color tn_color, structure and chain B and resi 101
color tn_color, structure and chain B and resi 103
# FN (low CSP in Binding Site): 2 residues
color fn_color, structure and chain B and resi 105
color fn_color, structure and chain B and resi 106
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
# TP (Sig. CSP in Union Site): 2
# FP (Sig. CSP -- Allosteric): 28
# TN (low CSP -- Allosteric): 40
# FN (low CSP in Union Site): 2
# Residues without CSP data: 34
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2M56_1/2M56_csp.pdb
