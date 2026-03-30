reinitialize
load ./outputs/2JMF/2JMF_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color green, structure and chain B
# Delta SASA range: 0.00 to 4.92 Å²
spectrum b, blue_white_red, structure and chain A, minimum=0.00, maximum=4.92
set cartoon_transparency, 0.2, structure
set cartoon_fancy_helices, 1
set cartoon_ring_mode, 1
# Colorbar information:
# Blue regions: Low delta SASA (minimal occlusion)
# White regions: Medium delta SASA (moderate occlusion)
# Red regions: High delta SASA (strong occlusion)
# Delta SASA Heatmap Analysis Summary:
# Receptor chain: A
# Ligand chain: B
# Occluded residues: 5
# Fraction occluded: 11.63%
# SASA threshold: 0.0 Å²
# Average percent burial: 0.8%
# Heatmap: Blue (low ΔSASA) -> White -> Red (high ΔSASA)
# Backbone colored by delta SASA values (protein-peptide occlusion)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2JMF/2JMF_delta_sasa.pdb
