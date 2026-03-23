reinitialize
load ./outputs/2mwp/2mwp_csp.pdb, structure
hide everything, structure
show cartoon, structure
color cyan, structure and chain B
set cartoon_tube_radius, 0.45, structure and chain B
color gray, structure and chain A
set_color tp_color, [0.1804, 0.8000, 0.4431]
set_color tn_color, [0.2039, 0.5961, 0.8588]
set_color fp_color, [0.6078, 0.3490, 0.7137]
set_color fn_color, [0.9529, 0.6118, 0.0706]
# TP (Sig. CSP in Binding Site): 8 residues
color tp_color, structure and chain A and resi 1500
color tp_color, structure and chain A and resi 1502
color tp_color, structure and chain A and resi 1521
color tp_color, structure and chain A and resi 1523
color tp_color, structure and chain A and resi 1547
color tp_color, structure and chain A and resi 1553
color tp_color, structure and chain A and resi 1584
color tp_color, structure and chain A and resi 1587
# FP (Sig. CSP -- Allosteric): 47 residues
color fp_color, structure and chain A and resi 1485
color fp_color, structure and chain A and resi 1486
color fp_color, structure and chain A and resi 1492
color fp_color, structure and chain A and resi 1493
color fp_color, structure and chain A and resi 1494
color fp_color, structure and chain A and resi 1501
color fp_color, structure and chain A and resi 1503
color fp_color, structure and chain A and resi 1504
color fp_color, structure and chain A and resi 1505
color fp_color, structure and chain A and resi 1506
color fp_color, structure and chain A and resi 1507
color fp_color, structure and chain A and resi 1509
color fp_color, structure and chain A and resi 1517
color fp_color, structure and chain A and resi 1518
color fp_color, structure and chain A and resi 1519
color fp_color, structure and chain A and resi 1520
color fp_color, structure and chain A and resi 1522
color fp_color, structure and chain A and resi 1524
color fp_color, structure and chain A and resi 1525
color fp_color, structure and chain A and resi 1526
color fp_color, structure and chain A and resi 1527
color fp_color, structure and chain A and resi 1529
color fp_color, structure and chain A and resi 1530
color fp_color, structure and chain A and resi 1531
color fp_color, structure and chain A and resi 1532
color fp_color, structure and chain A and resi 1533
color fp_color, structure and chain A and resi 1535
color fp_color, structure and chain A and resi 1536
color fp_color, structure and chain A and resi 1543
color fp_color, structure and chain A and resi 1546
color fp_color, structure and chain A and resi 1548
color fp_color, structure and chain A and resi 1550
color fp_color, structure and chain A and resi 1552
color fp_color, structure and chain A and resi 1554
color fp_color, structure and chain A and resi 1555
color fp_color, structure and chain A and resi 1568
color fp_color, structure and chain A and resi 1589
color fp_color, structure and chain A and resi 1590
color fp_color, structure and chain A and resi 1593
color fp_color, structure and chain A and resi 1594
color fp_color, structure and chain A and resi 1596
color fp_color, structure and chain A and resi 1597
color fp_color, structure and chain A and resi 1598
color fp_color, structure and chain A and resi 1599
color fp_color, structure and chain A and resi 1601
color fp_color, structure and chain A and resi 1602
color fp_color, structure and chain A and resi 1603
# TN (low CSP -- Allosteric): 52 residues
color tn_color, structure and chain A and resi 1487
color tn_color, structure and chain A and resi 1488
color tn_color, structure and chain A and resi 1489
color tn_color, structure and chain A and resi 1490
color tn_color, structure and chain A and resi 1491
color tn_color, structure and chain A and resi 1508
color tn_color, structure and chain A and resi 1510
color tn_color, structure and chain A and resi 1511
color tn_color, structure and chain A and resi 1512
color tn_color, structure and chain A and resi 1513
color tn_color, structure and chain A and resi 1514
color tn_color, structure and chain A and resi 1515
color tn_color, structure and chain A and resi 1516
color tn_color, structure and chain A and resi 1528
color tn_color, structure and chain A and resi 1534
color tn_color, structure and chain A and resi 1538
color tn_color, structure and chain A and resi 1540
color tn_color, structure and chain A and resi 1541
color tn_color, structure and chain A and resi 1542
color tn_color, structure and chain A and resi 1544
color tn_color, structure and chain A and resi 1545
color tn_color, structure and chain A and resi 1556
color tn_color, structure and chain A and resi 1557
color tn_color, structure and chain A and resi 1558
color tn_color, structure and chain A and resi 1559
color tn_color, structure and chain A and resi 1560
color tn_color, structure and chain A and resi 1561
color tn_color, structure and chain A and resi 1562
color tn_color, structure and chain A and resi 1563
color tn_color, structure and chain A and resi 1564
color tn_color, structure and chain A and resi 1567
color tn_color, structure and chain A and resi 1569
color tn_color, structure and chain A and resi 1570
color tn_color, structure and chain A and resi 1571
color tn_color, structure and chain A and resi 1572
color tn_color, structure and chain A and resi 1573
color tn_color, structure and chain A and resi 1574
color tn_color, structure and chain A and resi 1576
color tn_color, structure and chain A and resi 1577
color tn_color, structure and chain A and resi 1578
color tn_color, structure and chain A and resi 1579
color tn_color, structure and chain A and resi 1580
color tn_color, structure and chain A and resi 1581
color tn_color, structure and chain A and resi 1582
color tn_color, structure and chain A and resi 1583
color tn_color, structure and chain A and resi 1585
color tn_color, structure and chain A and resi 1586
color tn_color, structure and chain A and resi 1588
color tn_color, structure and chain A and resi 1591
color tn_color, structure and chain A and resi 1592
color tn_color, structure and chain A and resi 1595
color tn_color, structure and chain A and resi 1600
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
# TP (Sig. CSP in Union Site): 8
# FP (Sig. CSP -- Allosteric): 47
# TN (low CSP -- Allosteric): 52
# FN (low CSP in Union Site): 0
# Residues without CSP data: 16
# Structure colored by CSP classification (TP/FP/TN/FN)
# Modified PDB file: /Users/tiburon/Desktop/CSP_UBQ/outputs/2mwp/2mwp_csp.pdb
