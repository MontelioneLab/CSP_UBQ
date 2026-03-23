reinitialize
load ./outputs/2MZP/2MZP_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain I
set cartoon_tube_radius, 0.45, structure and chain I
color gray, structure and chain C
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 4 residues
color tp_color, structure and chain C and resi 65
color tp_color, structure and chain C and resi 67
color tp_color, structure and chain C and resi 69
color tp_color, structure and chain C and resi 71
# FP (Sig. CSP -- Allosteric): 34 residues
color fp_color, structure and chain C and resi 3
color fp_color, structure and chain C and resi 5
color fp_color, structure and chain C and resi 20
color fp_color, structure and chain C and resi 25
color fp_color, structure and chain C and resi 27
color fp_color, structure and chain C and resi 28
color fp_color, structure and chain C and resi 30
color fp_color, structure and chain C and resi 31
color fp_color, structure and chain C and resi 32
color fp_color, structure and chain C and resi 36
color fp_color, structure and chain C and resi 38
color fp_color, structure and chain C and resi 41
color fp_color, structure and chain C and resi 43
color fp_color, structure and chain C and resi 45
color fp_color, structure and chain C and resi 50
color fp_color, structure and chain C and resi 51
color fp_color, structure and chain C and resi 53
color fp_color, structure and chain C and resi 57
color fp_color, structure and chain C and resi 58
color fp_color, structure and chain C and resi 60
color fp_color, structure and chain C and resi 61
color fp_color, structure and chain C and resi 64
color fp_color, structure and chain C and resi 66
color fp_color, structure and chain C and resi 68
color fp_color, structure and chain C and resi 70
color fp_color, structure and chain C and resi 72
color fp_color, structure and chain C and resi 73
color fp_color, structure and chain C and resi 78
color fp_color, structure and chain C and resi 80
color fp_color, structure and chain C and resi 81
color fp_color, structure and chain C and resi 82
color fp_color, structure and chain C and resi 84
color fp_color, structure and chain C and resi 87
color fp_color, structure and chain C and resi 89
# TN (low CSP -- Allosteric): 42 residues
color tn_color, structure and chain C and resi 4
color tn_color, structure and chain C and resi 6
color tn_color, structure and chain C and resi 7
color tn_color, structure and chain C and resi 8
color tn_color, structure and chain C and resi 9
color tn_color, structure and chain C and resi 10
color tn_color, structure and chain C and resi 11
color tn_color, structure and chain C and resi 12
color tn_color, structure and chain C and resi 13
color tn_color, structure and chain C and resi 14
color tn_color, structure and chain C and resi 15
color tn_color, structure and chain C and resi 16
color tn_color, structure and chain C and resi 17
color tn_color, structure and chain C and resi 18
color tn_color, structure and chain C and resi 19
color tn_color, structure and chain C and resi 21
color tn_color, structure and chain C and resi 22
color tn_color, structure and chain C and resi 24
color tn_color, structure and chain C and resi 29
color tn_color, structure and chain C and resi 33
color tn_color, structure and chain C and resi 34
color tn_color, structure and chain C and resi 35
color tn_color, structure and chain C and resi 37
color tn_color, structure and chain C and resi 39
color tn_color, structure and chain C and resi 40
color tn_color, structure and chain C and resi 42
color tn_color, structure and chain C and resi 44
color tn_color, structure and chain C and resi 46
color tn_color, structure and chain C and resi 47
color tn_color, structure and chain C and resi 49
color tn_color, structure and chain C and resi 55
color tn_color, structure and chain C and resi 56
color tn_color, structure and chain C and resi 59
color tn_color, structure and chain C and resi 62
color tn_color, structure and chain C and resi 63
color tn_color, structure and chain C and resi 74
color tn_color, structure and chain C and resi 75
color tn_color, structure and chain C and resi 77
color tn_color, structure and chain C and resi 79
color tn_color, structure and chain C and resi 83
color tn_color, structure and chain C and resi 86
color tn_color, structure and chain C and resi 88
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain C and resi 26
color fn_color, structure and chain C and resi 48
color fn_color, structure and chain C and resi 76
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
# Receptor chain: C
# Ligand chain: I
# TP (Sig. CSP in Union Site): 4
# FP (Sig. CSP -- Allosteric): 34
# TN (low CSP -- Allosteric): 42
# FN (low CSP in Union Site): 3
# Residues without CSP data: 6
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2MZP/2MZP_csp.pdb
