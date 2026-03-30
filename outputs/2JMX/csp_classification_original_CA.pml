reinitialize
load ./outputs/2JMX/2JMX_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# FP (Sig. CSP -- Allosteric): 15 residues
color fp_color, structure and chain A and resi 13
color fp_color, structure and chain A and resi 19
color fp_color, structure and chain A and resi 26
color fp_color, structure and chain A and resi 30
color fp_color, structure and chain A and resi 33
color fp_color, structure and chain A and resi 41
color fp_color, structure and chain A and resi 43
color fp_color, structure and chain A and resi 51
color fp_color, structure and chain A and resi 83
color fp_color, structure and chain A and resi 95
color fp_color, structure and chain A and resi 98
color fp_color, structure and chain A and resi 103
color fp_color, structure and chain A and resi 107
color fp_color, structure and chain A and resi 110
color fp_color, structure and chain A and resi 111
# TN (low CSP -- Allosteric): 27 residues
color tn_color, structure and chain A and resi 2
color tn_color, structure and chain A and resi 27
color tn_color, structure and chain A and resi 34
color tn_color, structure and chain A and resi 39
color tn_color, structure and chain A and resi 42
color tn_color, structure and chain A and resi 46
color tn_color, structure and chain A and resi 47
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 53
color tn_color, structure and chain A and resi 54
color tn_color, structure and chain A and resi 55
color tn_color, structure and chain A and resi 59
color tn_color, structure and chain A and resi 60
color tn_color, structure and chain A and resi 61
color tn_color, structure and chain A and resi 62
color tn_color, structure and chain A and resi 68
color tn_color, structure and chain A and resi 70
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 74
color tn_color, structure and chain A and resi 76
color tn_color, structure and chain A and resi 93
color tn_color, structure and chain A and resi 106
color tn_color, structure and chain A and resi 114
color tn_color, structure and chain A and resi 116
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 119
color tn_color, structure and chain A and resi 120
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
# TP (Sig. CSP in Union Site): 0
# FP (Sig. CSP -- Allosteric): 15
# TN (low CSP -- Allosteric): 27
# FN (low CSP in Union Site): 0
# Residues without CSP data: 78
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2JMX/2JMX_csp.pdb
