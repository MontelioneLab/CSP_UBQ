"""
Central configuration for the CSP pipeline.

This module centralizes paths, thresholds, and IO settings. It is imported by
all other pipeline modules. Values can be overridden via environment variables
or CLI flags in `pipeline.py`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class Paths:
    workspace_root: str = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    input_csv: str = os.path.join(workspace_root, "CSP_UBQ.csv")
    cs_cache_dir: str = os.path.join(workspace_root, "CS_Lists")
    pdb_cache_dir: str = os.path.join(workspace_root, "PDB_FILES")
    outputs_dir: str = os.path.join(workspace_root, "outputs")
    pymol_views_dir: str = os.path.join(workspace_root, "pymol_views")


@dataclass(frozen=True)
class Thresholds:
    # Outlier detection settings for iterative method
    outlier_z_score: float = 3.0  # Z-score for outlier detection
    max_outlier_iterations: int = 10  # Max iterations for outlier removal
    max_outlier_fraction: float = 0.2  # Max 20% of residues removed
    
    # Final significance threshold settings
    significance_z_score: float = 0.0  # Z-score for final threshold (mean + z*SD)
    
    # Legacy option: absolute cutoff overrides everything
    absolute_cutoff: Optional[float] = None


@dataclass(frozen=True)
class Network:
    # Timeouts in seconds
    connect_timeout: float = 10.0
    read_timeout: float = 60.0
    retries: int = 3

    # BMRB FTP server URL templates for NMR-STAR files
    bmrb_ftp_url_template: str = "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{id}/bmr{id}_21.str"
    bmrb_ftp_fallback_url_template: str = "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{id}/bmr{id}_3.str"
    rcsb_pdb_url_template: str = "https://files.rcsb.org/download/{pdb_id}.pdb"


@dataclass(frozen=True)
class Alignment:
    # Global alignment scoring (manual NW fallback if Biopython not available)
    match_score: int = 2
    mismatch_penalty: int = -1
    gap_open_penalty: int = -2
    gap_extend_penalty: int = -1


@dataclass(frozen=True)
class Compute:
    # Combined CSP uses multipliers on Δδ (ppm) before squaring, same form as legacy (1/7, 1/4):
    # N/H: sqrt(1/2 * (Δδ_H^2 + (wN*Δδ_N)^2))
    # N/H/CA: sqrt(1/3 * (Δδ_H^2 + (wN*Δδ_N)^2 + (wCA*Δδ_Cα)^2))
    csp_delta_n_scale: float = 0.14
    csp_delta_ca_scale: float = 0.3


@dataclass(frozen=True)
class Referencing:
    # Referencing method: 'mean' (legacy) or 'grid' (new)
    method: str = "grid"

    # Grid search parameters for offsets
    grid_h_min: float = -0.12
    grid_h_max: float = 0.12
    grid_h_step: float = 0.01

    grid_n_min: float = -1.2
    grid_n_max: float = 1.2
    grid_n_step: float = 0.05

    # CA grid search parameters (scaled similarly to N range)
    grid_ca_min: float = -6.0
    grid_ca_max: float = 6.0
    grid_ca_step: float = 0.2

    # CSP cutoff used during grid search (count CSP < cutoff)
    grid_cutoff: float = 0.05

    # Persistence and visualization
    cache_results: bool = True
    save_heatmap: bool = True

@dataclass(frozen=True)
class Concurrency:
    workers: int = max(1, int(os.environ.get("CSP_WORKERS", "4")))


@dataclass(frozen=True)
class SASAAnalysis:
    # SASA threshold for determining occluded residues (Å²)
    sasa_threshold: float = 0.0


@dataclass(frozen=True)
class CADistanceAnalysis:
    # CA-CA distance threshold for proximity filter (Å)
    ca_distance_threshold: float = 6.0
    # Min inter-chain atom-atom distance for direct binding site (Å)
    direct_contact_threshold: float = 2.0


@dataclass(frozen=True)
class ClassificationColors:
    """Hex codes for confusion matrix visualization (TP/FP/TN/FN)."""
    TP: str = "#2ecc71"   # Green
    FP: str = "#9b59b6"   # Purple
    TN: str = "#3498db"   # Blue
    FN: str = "#f39c12"   # Orange


def hex_to_rgb01(hex_color: str) -> Tuple[float, float, float]:
    """Convert hex color to RGB tuple in 0–1 range (e.g. for PyMOL `set_color`)."""
    hex_clean = hex_color.lstrip("#")
    if len(hex_clean) != 6:
        return (0.5, 0.5, 0.5)
    return (
        int(hex_clean[0:2], 16) / 255.0,
        int(hex_clean[2:4], 16) / 255.0,
        int(hex_clean[4:6], 16) / 255.0,
    )


# Export a single config object for convenience
paths = Paths()
thresholds = Thresholds()
network = Network()
alignment = Alignment()
compute = Compute()
concurrency = Concurrency()
sasa_analysis = SASAAnalysis()
ca_distance_analysis = CADistanceAnalysis()
classification_colors = ClassificationColors()


def ensure_directories(base_output: Optional[str] = None) -> None:
    """Create required directories if they do not exist.

    Parameters
    ----------
    base_output: Optional[str]
        If provided, overrides the default outputs directory for this run.
    """
    os.makedirs(paths.cs_cache_dir, exist_ok=True)
    os.makedirs(paths.pdb_cache_dir, exist_ok=True)
    os.makedirs(paths.pymol_views_dir, exist_ok=True)
    os.makedirs(base_output or paths.outputs_dir, exist_ok=True)


