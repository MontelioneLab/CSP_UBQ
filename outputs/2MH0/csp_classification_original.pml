reinitialize
load ./outputs/2MH0/2MH0_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain A
set cartoon_tube_radius, 0.45, structure and chain A
color gray, structure and chain B
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 10 residues
color tp_color, structure and chain B and resi 1735
color tp_color, structure and chain B and resi 1738
color tp_color, structure and chain B and resi 1760
color tp_color, structure and chain B and resi 1761
color tp_color, structure and chain B and resi 1763
color tp_color, structure and chain B and resi 1764
color tp_color, structure and chain B and resi 1781
color tp_color, structure and chain B and resi 1783
color tp_color, structure and chain B and resi 1784
color tp_color, structure and chain B and resi 1788
# FP (Sig. CSP -- Allosteric): 32 residues
color fp_color, structure and chain B and resi 1723
color fp_color, structure and chain B and resi 1724
color fp_color, structure and chain B and resi 1725
color fp_color, structure and chain B and resi 1726
color fp_color, structure and chain B and resi 1728
color fp_color, structure and chain B and resi 1730
color fp_color, structure and chain B and resi 1740
color fp_color, structure and chain B and resi 1742
color fp_color, structure and chain B and resi 1745
color fp_color, structure and chain B and resi 1747
color fp_color, structure and chain B and resi 1748
color fp_color, structure and chain B and resi 1749
color fp_color, structure and chain B and resi 1750
color fp_color, structure and chain B and resi 1751
color fp_color, structure and chain B and resi 1754
color fp_color, structure and chain B and resi 1757
color fp_color, structure and chain B and resi 1762
color fp_color, structure and chain B and resi 1766
color fp_color, structure and chain B and resi 1768
color fp_color, structure and chain B and resi 1771
color fp_color, structure and chain B and resi 1772
color fp_color, structure and chain B and resi 1775
color fp_color, structure and chain B and resi 1786
color fp_color, structure and chain B and resi 1790
color fp_color, structure and chain B and resi 1793
color fp_color, structure and chain B and resi 1799
color fp_color, structure and chain B and resi 1801
color fp_color, structure and chain B and resi 1807
color fp_color, structure and chain B and resi 1809
color fp_color, structure and chain B and resi 1810
color fp_color, structure and chain B and resi 1811
color fp_color, structure and chain B and resi 1812
# TN (low CSP -- Allosteric): 37 residues
color tn_color, structure and chain B and resi 1729
color tn_color, structure and chain B and resi 1732
color tn_color, structure and chain B and resi 1733
color tn_color, structure and chain B and resi 1736
color tn_color, structure and chain B and resi 1737
color tn_color, structure and chain B and resi 1739
color tn_color, structure and chain B and resi 1741
color tn_color, structure and chain B and resi 1743
color tn_color, structure and chain B and resi 1744
color tn_color, structure and chain B and resi 1746
color tn_color, structure and chain B and resi 1752
color tn_color, structure and chain B and resi 1753
color tn_color, structure and chain B and resi 1755
color tn_color, structure and chain B and resi 1758
color tn_color, structure and chain B and resi 1759
color tn_color, structure and chain B and resi 1765
color tn_color, structure and chain B and resi 1769
color tn_color, structure and chain B and resi 1770
color tn_color, structure and chain B and resi 1773
color tn_color, structure and chain B and resi 1774
color tn_color, structure and chain B and resi 1777
color tn_color, structure and chain B and resi 1778
color tn_color, structure and chain B and resi 1779
color tn_color, structure and chain B and resi 1782
color tn_color, structure and chain B and resi 1785
color tn_color, structure and chain B and resi 1789
color tn_color, structure and chain B and resi 1791
color tn_color, structure and chain B and resi 1792
color tn_color, structure and chain B and resi 1794
color tn_color, structure and chain B and resi 1795
color tn_color, structure and chain B and resi 1796
color tn_color, structure and chain B and resi 1797
color tn_color, structure and chain B and resi 1798
color tn_color, structure and chain B and resi 1800
color tn_color, structure and chain B and resi 1803
color tn_color, structure and chain B and resi 1806
color tn_color, structure and chain B and resi 1808
# FN (low CSP in Binding Site): 4 residues
color fn_color, structure and chain B and resi 1731
color fn_color, structure and chain B and resi 1734
color fn_color, structure and chain B and resi 1767
color fn_color, structure and chain B and resi 1787
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
# TP (Sig. CSP in Union Site): 10
# FP (Sig. CSP -- Allosteric): 32
# TN (low CSP -- Allosteric): 37
# FN (low CSP in Union Site): 4
# Residues without CSP data: 9
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2MH0/2MH0_csp.pdb
