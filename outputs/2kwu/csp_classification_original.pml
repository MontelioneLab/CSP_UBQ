reinitialize
load ./outputs/2kwu/2kwu_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 9 residues
color tp_color, structure and chain B and resi 7
color tp_color, structure and chain B and resi 8
color tp_color, structure and chain B and resi 9
color tp_color, structure and chain B and resi 10
color tp_color, structure and chain B and resi 42
color tp_color, structure and chain B and resi 49
color tp_color, structure and chain B and resi 68
color tp_color, structure and chain B and resi 72
color tp_color, structure and chain B and resi 73
# FP (Sig. CSP -- Allosteric): 18 residues
color fp_color, structure and chain B and resi 2
color fp_color, structure and chain B and resi 3
color fp_color, structure and chain B and resi 6
color fp_color, structure and chain B and resi 13
color fp_color, structure and chain B and resi 18
color fp_color, structure and chain B and resi 20
color fp_color, structure and chain B and resi 28
color fp_color, structure and chain B and resi 31
color fp_color, structure and chain B and resi 34
color fp_color, structure and chain B and resi 40
color fp_color, structure and chain B and resi 43
color fp_color, structure and chain B and resi 45
color fp_color, structure and chain B and resi 47
color fp_color, structure and chain B and resi 48
color fp_color, structure and chain B and resi 50
color fp_color, structure and chain B and resi 60
color fp_color, structure and chain B and resi 69
color fp_color, structure and chain B and resi 71
# TN (low CSP -- Allosteric): 41 residues
color tn_color, structure and chain B and resi 4
color tn_color, structure and chain B and resi 5
color tn_color, structure and chain B and resi 11
color tn_color, structure and chain B and resi 12
color tn_color, structure and chain B and resi 14
color tn_color, structure and chain B and resi 15
color tn_color, structure and chain B and resi 16
color tn_color, structure and chain B and resi 17
color tn_color, structure and chain B and resi 21
color tn_color, structure and chain B and resi 22
color tn_color, structure and chain B and resi 23
color tn_color, structure and chain B and resi 25
color tn_color, structure and chain B and resi 26
color tn_color, structure and chain B and resi 27
color tn_color, structure and chain B and resi 29
color tn_color, structure and chain B and resi 30
color tn_color, structure and chain B and resi 32
color tn_color, structure and chain B and resi 33
color tn_color, structure and chain B and resi 35
color tn_color, structure and chain B and resi 36
color tn_color, structure and chain B and resi 39
color tn_color, structure and chain B and resi 41
color tn_color, structure and chain B and resi 46
color tn_color, structure and chain B and resi 51
color tn_color, structure and chain B and resi 52
color tn_color, structure and chain B and resi 54
color tn_color, structure and chain B and resi 55
color tn_color, structure and chain B and resi 56
color tn_color, structure and chain B and resi 57
color tn_color, structure and chain B and resi 58
color tn_color, structure and chain B and resi 59
color tn_color, structure and chain B and resi 61
color tn_color, structure and chain B and resi 62
color tn_color, structure and chain B and resi 63
color tn_color, structure and chain B and resi 64
color tn_color, structure and chain B and resi 65
color tn_color, structure and chain B and resi 66
color tn_color, structure and chain B and resi 67
color tn_color, structure and chain B and resi 70
color tn_color, structure and chain B and resi 74
color tn_color, structure and chain B and resi 76
# FN (low CSP in Binding Site): 2 residues
color fn_color, structure and chain B and resi 44
color fn_color, structure and chain B and resi 75
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
# TP (Sig. CSP in Union Site): 9
# FP (Sig. CSP -- Allosteric): 18
# TN (low CSP -- Allosteric): 41
# FN (low CSP in Union Site): 2
# Residues without CSP data: 6
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2kwu/2kwu_csp.pdb
