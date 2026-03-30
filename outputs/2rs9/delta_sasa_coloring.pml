reinitialize
load ./outputs/2rs9/2rs9_delta_sasa.pdb, structure
hide everything, structure
show cartoon, structure
color green, structure and chain A
# Delta SASA range: 0.00 to 3.60 Å²
spectrum b, blue_white_red, structure and chain B, minimum=0.00, maximum=3.60
set cartoon_transparency, 0.2, structure
set cartoon_fancy_helices, 1
set cartoon_ring_mode, 1
# Colorbar information:
# Blue regions: Low delta SASA (minimal occlusion)
# White regions: Medium delta SASA (moderate occlusion)
# Red regions: High delta SASA (strong occlusion)
# Delta SASA Heatmap Analysis Summary:
# Receptor chain: B
# Ligand chain: A
# Occluded residues: 7
# Fraction occluded: 5.79%
# SASA threshold: 0.0 Å²
# Average percent burial: 0.4%
# Heatmap: Blue (low ΔSASA) -> White -> Red (high ΔSASA)
# Backbone colored by delta SASA values (protein-peptide occlusion)
# Modified PDB file: /Users/tiburon/Desktop/new_CSP_UBQ/CSP_UBQ/outputs/2rs9/2rs9_delta_sasa.pdb
