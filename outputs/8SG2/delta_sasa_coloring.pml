reinitialize
load ./outputs/8SG2/8SG2_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color green, structure and chain B
# Delta SASA range: 0.00 to 10.85 Å²
spectrum b, blue_white_red, structure and chain A, minimum=0.00, maximum=10.85
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
# Occluded residues: 34
# Fraction occluded: 20.86%
# SASA threshold: 0.0 Å²
# Average percent burial: 3.5%
# Heatmap: Blue (low ΔSASA) -> White -> Red (high ΔSASA)
# Backbone colored by delta SASA values (protein-peptide occlusion)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/8SG2/8SG2_delta_sasa.pdb
