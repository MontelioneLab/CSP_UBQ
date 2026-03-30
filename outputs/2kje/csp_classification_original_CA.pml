reinitialize
load ./outputs/2kje/2kje_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 14 residues
color tp_color, structure and chain A and resi 1766
color tp_color, structure and chain A and resi 1772
color tp_color, structure and chain A and resi 1775
color tp_color, structure and chain A and resi 1778
color tp_color, structure and chain A and resi 1793
color tp_color, structure and chain A and resi 1796
color tp_color, structure and chain A and resi 1798
color tp_color, structure and chain A and resi 1801
color tp_color, structure and chain A and resi 1802
color tp_color, structure and chain A and resi 1817
color tp_color, structure and chain A and resi 1826
color tp_color, structure and chain A and resi 1834
color tp_color, structure and chain A and resi 1839
color tp_color, structure and chain A and resi 1844
# FP (Sig. CSP -- Allosteric): 22 residues
color fp_color, structure and chain A and resi 1767
color fp_color, structure and chain A and resi 1768
color fp_color, structure and chain A and resi 1780
color fp_color, structure and chain A and resi 1783
color fp_color, structure and chain A and resi 1788
color fp_color, structure and chain A and resi 1790
color fp_color, structure and chain A and resi 1797
color fp_color, structure and chain A and resi 1803
color fp_color, structure and chain A and resi 1808
color fp_color, structure and chain A and resi 1812
color fp_color, structure and chain A and resi 1815
color fp_color, structure and chain A and resi 1816
color fp_color, structure and chain A and resi 1827
color fp_color, structure and chain A and resi 1831
color fp_color, structure and chain A and resi 1833
color fp_color, structure and chain A and resi 1836
color fp_color, structure and chain A and resi 1845
color fp_color, structure and chain A and resi 1846
color fp_color, structure and chain A and resi 1847
color fp_color, structure and chain A and resi 1848
color fp_color, structure and chain A and resi 1849
color fp_color, structure and chain A and resi 1850
# TN (low CSP -- Allosteric): 21 residues
color tn_color, structure and chain A and resi 1771
color tn_color, structure and chain A and resi 1774
color tn_color, structure and chain A and resi 1777
color tn_color, structure and chain A and resi 1781
color tn_color, structure and chain A and resi 1784
color tn_color, structure and chain A and resi 1785
color tn_color, structure and chain A and resi 1787
color tn_color, structure and chain A and resi 1789
color tn_color, structure and chain A and resi 1800
color tn_color, structure and chain A and resi 1804
color tn_color, structure and chain A and resi 1806
color tn_color, structure and chain A and resi 1807
color tn_color, structure and chain A and resi 1811
color tn_color, structure and chain A and resi 1821
color tn_color, structure and chain A and resi 1824
color tn_color, structure and chain A and resi 1828
color tn_color, structure and chain A and resi 1829
color tn_color, structure and chain A and resi 1832
color tn_color, structure and chain A and resi 1835
color tn_color, structure and chain A and resi 1838
color tn_color, structure and chain A and resi 1841
# FN (low CSP in Binding Site): 14 residues
color fn_color, structure and chain A and resi 1773
color fn_color, structure and chain A and resi 1779
color fn_color, structure and chain A and resi 1782
color fn_color, structure and chain A and resi 1786
color fn_color, structure and chain A and resi 1791
color fn_color, structure and chain A and resi 1792
color fn_color, structure and chain A and resi 1795
color fn_color, structure and chain A and resi 1799
color fn_color, structure and chain A and resi 1805
color fn_color, structure and chain A and resi 1809
color fn_color, structure and chain A and resi 1819
color fn_color, structure and chain A and resi 1820
color fn_color, structure and chain A and resi 1825
color fn_color, structure and chain A and resi 1830
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
# TP (Sig. CSP in Union Site): 14
# FP (Sig. CSP -- Allosteric): 22
# TN (low CSP -- Allosteric): 21
# FN (low CSP in Union Site): 14
# Residues without CSP data: 21
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2kje/2kje_csp.pdb
