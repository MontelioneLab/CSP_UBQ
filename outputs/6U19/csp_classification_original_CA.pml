reinitialize
load ./outputs/6U19/6U19_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain B and resi 29
color tp_color, structure and chain B and resi 30
color tp_color, structure and chain B and resi 33
color tp_color, structure and chain B and resi 66
color tp_color, structure and chain B and resi 69
color tp_color, structure and chain B and resi 70
color tp_color, structure and chain B and resi 76
color tp_color, structure and chain B and resi 77
# FP (Sig. CSP -- Allosteric): 12 residues
color fp_color, structure and chain B and resi 31
color fp_color, structure and chain B and resi 32
color fp_color, structure and chain B and resi 34
color fp_color, structure and chain B and resi 35
color fp_color, structure and chain B and resi 36
color fp_color, structure and chain B and resi 37
color fp_color, structure and chain B and resi 38
color fp_color, structure and chain B and resi 42
color fp_color, structure and chain B and resi 67
color fp_color, structure and chain B and resi 73
color fp_color, structure and chain B and resi 74
color fp_color, structure and chain B and resi 86
# TN (low CSP -- Allosteric): 30 residues
color tn_color, structure and chain B and resi 39
color tn_color, structure and chain B and resi 40
color tn_color, structure and chain B and resi 41
color tn_color, structure and chain B and resi 43
color tn_color, structure and chain B and resi 45
color tn_color, structure and chain B and resi 46
color tn_color, structure and chain B and resi 47
color tn_color, structure and chain B and resi 48
color tn_color, structure and chain B and resi 50
color tn_color, structure and chain B and resi 51
color tn_color, structure and chain B and resi 52
color tn_color, structure and chain B and resi 53
color tn_color, structure and chain B and resi 55
color tn_color, structure and chain B and resi 56
color tn_color, structure and chain B and resi 57
color tn_color, structure and chain B and resi 59
color tn_color, structure and chain B and resi 60
color tn_color, structure and chain B and resi 61
color tn_color, structure and chain B and resi 62
color tn_color, structure and chain B and resi 63
color tn_color, structure and chain B and resi 64
color tn_color, structure and chain B and resi 68
color tn_color, structure and chain B and resi 71
color tn_color, structure and chain B and resi 72
color tn_color, structure and chain B and resi 75
color tn_color, structure and chain B and resi 78
color tn_color, structure and chain B and resi 80
color tn_color, structure and chain B and resi 81
color tn_color, structure and chain B and resi 82
color tn_color, structure and chain B and resi 84
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain B and resi 44
color fn_color, structure and chain B and resi 49
color fn_color, structure and chain B and resi 54
color fn_color, structure and chain B and resi 83
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 12
# TN (low CSP -- Allosteric): 30
# FN (low CSP in Union Site): 4
# Residues without CSP data: 10
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/6U19/6U19_csp.pdb
