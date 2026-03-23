reinitialize
load ./outputs/2JQR/2JQR_csp.pdb, structure
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
color tp_color, structure and chain A and resi 48
color tp_color, structure and chain A and resi 49
color tp_color, structure and chain A and resi 59
color tp_color, structure and chain A and resi 78
color tp_color, structure and chain A and resi 80
# FP (Sig. CSP -- Allosteric): 24 residues
color fp_color, structure and chain A and resi -4
color fp_color, structure and chain A and resi -3
color fp_color, structure and chain A and resi 4
color fp_color, structure and chain A and resi 24
color fp_color, structure and chain A and resi 26
color fp_color, structure and chain A and resi 27
color fp_color, structure and chain A and resi 30
color fp_color, structure and chain A and resi 31
color fp_color, structure and chain A and resi 39
color fp_color, structure and chain A and resi 43
color fp_color, structure and chain A and resi 44
color fp_color, structure and chain A and resi 50
color fp_color, structure and chain A and resi 61
color fp_color, structure and chain A and resi 66
color fp_color, structure and chain A and resi 70
color fp_color, structure and chain A and resi 73
color fp_color, structure and chain A and resi 77
color fp_color, structure and chain A and resi 83
color fp_color, structure and chain A and resi 84
color fp_color, structure and chain A and resi 89
color fp_color, structure and chain A and resi 96
color fp_color, structure and chain A and resi 102
color fp_color, structure and chain A and resi 103
color fp_color, structure and chain A and resi 104
# TN (low CSP -- Allosteric): 36 residues
color tn_color, structure and chain A and resi -2
color tn_color, structure and chain A and resi 2
color tn_color, structure and chain A and resi 3
color tn_color, structure and chain A and resi 5
color tn_color, structure and chain A and resi 6
color tn_color, structure and chain A and resi 7
color tn_color, structure and chain A and resi 8
color tn_color, structure and chain A and resi 9
color tn_color, structure and chain A and resi 12
color tn_color, structure and chain A and resi 19
color tn_color, structure and chain A and resi 20
color tn_color, structure and chain A and resi 21
color tn_color, structure and chain A and resi 22
color tn_color, structure and chain A and resi 23
color tn_color, structure and chain A and resi 25
color tn_color, structure and chain A and resi 29
color tn_color, structure and chain A and resi 32
color tn_color, structure and chain A and resi 33
color tn_color, structure and chain A and resi 34
color tn_color, structure and chain A and resi 36
color tn_color, structure and chain A and resi 37
color tn_color, structure and chain A and resi 38
color tn_color, structure and chain A and resi 42
color tn_color, structure and chain A and resi 51
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 60
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 72
color tn_color, structure and chain A and resi 74
color tn_color, structure and chain A and resi 85
color tn_color, structure and chain A and resi 86
color tn_color, structure and chain A and resi 87
color tn_color, structure and chain A and resi 88
color tn_color, structure and chain A and resi 90
color tn_color, structure and chain A and resi 92
color tn_color, structure and chain A and resi 95
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain A and resi 28
color fn_color, structure and chain A and resi 35
color fn_color, structure and chain A and resi 81
color fn_color, structure and chain A and resi 82
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
# FP (Sig. CSP -- Allosteric): 24
# TN (low CSP -- Allosteric): 36
# FN (low CSP in Union Site): 4
# Residues without CSP data: 40
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2JQR/2JQR_csp.pdb
