reinitialize
load ./outputs/7ZEY_2/7ZEY_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 11 residues
color tp_color, structure and chain B and resi 1588
color tp_color, structure and chain B and resi 1596
color tp_color, structure and chain B and resi 1601
color tp_color, structure and chain B and resi 1605
color tp_color, structure and chain B and resi 1606
color tp_color, structure and chain B and resi 1609
color tp_color, structure and chain B and resi 1612
color tp_color, structure and chain B and resi 1616
color tp_color, structure and chain B and resi 1617
color tp_color, structure and chain B and resi 1618
color tp_color, structure and chain B and resi 1620
# FP (Sig. CSP -- Allosteric): 12 residues
color fp_color, structure and chain B and resi 1574
color fp_color, structure and chain B and resi 1576
color fp_color, structure and chain B and resi 1577
color fp_color, structure and chain B and resi 1580
color fp_color, structure and chain B and resi 1582
color fp_color, structure and chain B and resi 1586
color fp_color, structure and chain B and resi 1603
color fp_color, structure and chain B and resi 1607
color fp_color, structure and chain B and resi 1610
color fp_color, structure and chain B and resi 1613
color fp_color, structure and chain B and resi 1622
color fp_color, structure and chain B and resi 1625
# TN (low CSP -- Allosteric): 20 residues
color tn_color, structure and chain B and resi 1565
color tn_color, structure and chain B and resi 1567
color tn_color, structure and chain B and resi 1575
color tn_color, structure and chain B and resi 1578
color tn_color, structure and chain B and resi 1579
color tn_color, structure and chain B and resi 1581
color tn_color, structure and chain B and resi 1583
color tn_color, structure and chain B and resi 1584
color tn_color, structure and chain B and resi 1585
color tn_color, structure and chain B and resi 1587
color tn_color, structure and chain B and resi 1592
color tn_color, structure and chain B and resi 1593
color tn_color, structure and chain B and resi 1594
color tn_color, structure and chain B and resi 1595
color tn_color, structure and chain B and resi 1597
color tn_color, structure and chain B and resi 1600
color tn_color, structure and chain B and resi 1611
color tn_color, structure and chain B and resi 1623
color tn_color, structure and chain B and resi 1626
color tn_color, structure and chain B and resi 1702
# FN (low CSP in Binding Site): 13 residues
color fn_color, structure and chain B and resi 1569
color fn_color, structure and chain B and resi 1589
color fn_color, structure and chain B and resi 1590
color fn_color, structure and chain B and resi 1591
color fn_color, structure and chain B and resi 1599
color fn_color, structure and chain B and resi 1602
color fn_color, structure and chain B and resi 1604
color fn_color, structure and chain B and resi 1608
color fn_color, structure and chain B and resi 1615
color fn_color, structure and chain B and resi 1619
color fn_color, structure and chain B and resi 1621
color fn_color, structure and chain B and resi 1624
color fn_color, structure and chain B and resi 1627
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
# TP (Sig. CSP in Union Site): 11
# FP (Sig. CSP -- Allosteric): 12
# TN (low CSP -- Allosteric): 20
# FN (low CSP in Union Site): 13
# Residues without CSP data: 9
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/7ZEY_2/7ZEY_csp.pdb
