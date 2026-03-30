reinitialize
load ./outputs/2NBV/2NBV_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# FP (Sig. CSP -- Allosteric): 12 residues
color fp_color, structure and chain A and resi 24
color fp_color, structure and chain A and resi 34
color fp_color, structure and chain A and resi 37
color fp_color, structure and chain A and resi 53
color fp_color, structure and chain A and resi 61
color fp_color, structure and chain A and resi 66
color fp_color, structure and chain A and resi 82
color fp_color, structure and chain A and resi 84
color fp_color, structure and chain A and resi 101
color fp_color, structure and chain A and resi 105
color fp_color, structure and chain A and resi 114
color fp_color, structure and chain A and resi 119
# TN (low CSP -- Allosteric): 20 residues
color tn_color, structure and chain A and resi 23
color tn_color, structure and chain A and resi 35
color tn_color, structure and chain A and resi 36
color tn_color, structure and chain A and resi 49
color tn_color, structure and chain A and resi 51
color tn_color, structure and chain A and resi 52
color tn_color, structure and chain A and resi 64
color tn_color, structure and chain A and resi 65
color tn_color, structure and chain A and resi 68
color tn_color, structure and chain A and resi 69
color tn_color, structure and chain A and resi 71
color tn_color, structure and chain A and resi 80
color tn_color, structure and chain A and resi 81
color tn_color, structure and chain A and resi 99
color tn_color, structure and chain A and resi 103
color tn_color, structure and chain A and resi 118
color tn_color, structure and chain A and resi 122
color tn_color, structure and chain A and resi 124
color tn_color, structure and chain A and resi 125
color tn_color, structure and chain A and resi 130
# FN (low CSP in Binding Site): 3 residues
color fn_color, structure and chain A and resi 55
color fn_color, structure and chain A and resi 78
color fn_color, structure and chain A and resi 100
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
# FP (Sig. CSP -- Allosteric): 12
# TN (low CSP -- Allosteric): 20
# FN (low CSP in Union Site): 3
# Residues without CSP data: 75
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2NBV/2NBV_csp.pdb
