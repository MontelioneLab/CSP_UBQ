reinitialize
load ./outputs/2mur_2/2mur_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 11 residues
color tp_color, structure and chain A and resi 28
color tp_color, structure and chain A and resi 32
color tp_color, structure and chain A and resi 33
color tp_color, structure and chain A and resi 34
color tp_color, structure and chain A and resi 35
color tp_color, structure and chain A and resi 36
color tp_color, structure and chain A and resi 38
color tp_color, structure and chain A and resi 39
color tp_color, structure and chain A and resi 40
color tp_color, structure and chain A and resi 41
color tp_color, structure and chain A and resi 44
# FP (Sig. CSP -- Allosteric): 5 residues
color fp_color, structure and chain A and resi 10
color fp_color, structure and chain A and resi 13
color fp_color, structure and chain A and resi 37
color fp_color, structure and chain A and resi 42
color fp_color, structure and chain A and resi 43
# TN (low CSP -- Allosteric): 14 residues
color tn_color, structure and chain A and resi 4
color tn_color, structure and chain A and resi 5
color tn_color, structure and chain A and resi 6
color tn_color, structure and chain A and resi 9
color tn_color, structure and chain A and resi 15
color tn_color, structure and chain A and resi 16
color tn_color, structure and chain A and resi 17
color tn_color, structure and chain A and resi 18
color tn_color, structure and chain A and resi 19
color tn_color, structure and chain A and resi 21
color tn_color, structure and chain A and resi 22
color tn_color, structure and chain A and resi 23
color tn_color, structure and chain A and resi 26
color tn_color, structure and chain A and resi 27
# FN (low CSP in Binding Site): 8 residues
color fn_color, structure and chain A and resi 7
color fn_color, structure and chain A and resi 8
color fn_color, structure and chain A and resi 11
color fn_color, structure and chain A and resi 14
color fn_color, structure and chain A and resi 24
color fn_color, structure and chain A and resi 25
color fn_color, structure and chain A and resi 30
color fn_color, structure and chain A and resi 31
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
# TP (Sig. CSP in Union Site): 11
# FP (Sig. CSP -- Allosteric): 5
# TN (low CSP -- Allosteric): 14
# FN (low CSP in Union Site): 8
# Residues without CSP data: 6
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2mur_2/2mur_csp.pdb
