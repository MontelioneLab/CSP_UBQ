reinitialize
load ./outputs/6K5R/6K5R_csp.pdb, structure
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
color tp_color, structure and chain A and resi 17
color tp_color, structure and chain A and resi 30
color tp_color, structure and chain A and resi 31
color tp_color, structure and chain A and resi 32
color tp_color, structure and chain A and resi 34
color tp_color, structure and chain A and resi 37
color tp_color, structure and chain A and resi 38
color tp_color, structure and chain A and resi 42
color tp_color, structure and chain A and resi 47
color tp_color, structure and chain A and resi 50
# FP (Sig. CSP -- Allosteric): 24 residues
color fp_color, structure and chain A and resi 16
color fp_color, structure and chain A and resi 18
color fp_color, structure and chain A and resi 19
color fp_color, structure and chain A and resi 20
color fp_color, structure and chain A and resi 23
color fp_color, structure and chain A and resi 25
color fp_color, structure and chain A and resi 36
color fp_color, structure and chain A and resi 40
color fp_color, structure and chain A and resi 41
color fp_color, structure and chain A and resi 43
color fp_color, structure and chain A and resi 45
color fp_color, structure and chain A and resi 46
color fp_color, structure and chain A and resi 48
color fp_color, structure and chain A and resi 62
color fp_color, structure and chain A and resi 67
color fp_color, structure and chain A and resi 68
color fp_color, structure and chain A and resi 70
color fp_color, structure and chain A and resi 72
color fp_color, structure and chain A and resi 74
color fp_color, structure and chain A and resi 76
color fp_color, structure and chain A and resi 81
color fp_color, structure and chain A and resi 82
color fp_color, structure and chain A and resi 84
color fp_color, structure and chain A and resi 87
# TN (low CSP -- Allosteric): 32 residues
color tn_color, structure and chain A and resi 21
color tn_color, structure and chain A and resi 22
color tn_color, structure and chain A and resi 24
color tn_color, structure and chain A and resi 26
color tn_color, structure and chain A and resi 27
color tn_color, structure and chain A and resi 49
color tn_color, structure and chain A and resi 51
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 56
color tn_color, structure and chain A and resi 57
color tn_color, structure and chain A and resi 58
color tn_color, structure and chain A and resi 59
color tn_color, structure and chain A and resi 60
color tn_color, structure and chain A and resi 61
color tn_color, structure and chain A and resi 63
color tn_color, structure and chain A and resi 64
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 69
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 75
color tn_color, structure and chain A and resi 77
color tn_color, structure and chain A and resi 78
color tn_color, structure and chain A and resi 79
color tn_color, structure and chain A and resi 80
color tn_color, structure and chain A and resi 83
color tn_color, structure and chain A and resi 85
color tn_color, structure and chain A and resi 86
color tn_color, structure and chain A and resi 88
color tn_color, structure and chain A and resi 89
# FN (low CSP in Binding Site): 2 residues
color fn_color, structure and chain A and resi 28
color fn_color, structure and chain A and resi 29
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
# FP (Sig. CSP -- Allosteric): 24
# TN (low CSP -- Allosteric): 32
# FN (low CSP in Union Site): 2
# Residues without CSP data: 9
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/6K5R/6K5R_csp.pdb
