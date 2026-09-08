from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from csp_reusable_snippets import (
    compute_atom_deltas_with_offset as compute_atom_deltas_with_offset_bundle,
    compute_ca_inclusive_csp,
    compute_ha_ca_csp,
    compute_nh_csp,
    compute_threshold_with_outlier_removal,
    run_offset_grid_search,
)

try:
    from .csp import (
        compute_atom_deltas_with_offset as compute_atom_deltas_with_offset_pipeline,
        compute_csp_A,
        compute_csp_from_aligned_sequences_ca,
        compute_csp_from_aligned_sequences_ha_ca,
        compute_threshold_with_outlier_removal as compute_threshold_with_outlier_removal_pipeline,
        collect_large_dn_exclusions,
    )
except Exception:
    from scripts.csp import (
        compute_atom_deltas_with_offset as compute_atom_deltas_with_offset_pipeline,
        compute_csp_A,
        compute_csp_from_aligned_sequences_ca,
        compute_csp_from_aligned_sequences_ha_ca,
        compute_threshold_with_outlier_removal as compute_threshold_with_outlier_removal_pipeline,
        collect_large_dn_exclusions,
    )


def test_compute_threshold_with_outlier_removal_matches_pipeline():
    values = [0.10, 0.12, 0.11, 0.09, 1.80]
    bundle = compute_threshold_with_outlier_removal(values, 1.5, 0.0, 10, 0.5)
    pipeline = compute_threshold_with_outlier_removal_pipeline(values, 1.5, 0.0, 10, 0.5)
    assert bundle.threshold == pytest.approx(pipeline.threshold)
    assert bundle.cleaned_values == pytest.approx(pipeline.cleaned_values)
    assert bundle.outliers_removed == pipeline.outliers_removed


def test_nh_formula_and_absolute_cutoff_override():
    result = compute_nh_csp(
        mapping=[(1, 1)],
        apo_seq="M",
        holo_seq="M",
        H_apo={1: 8.40},
        N_apo={1: 121.0},
        H_holo={1: 8.30},
        N_holo={1: 119.0},
        threshold_config={"absolute_cutoff": 0.15},
        enable_referencing=False,
    )
    expected = math.sqrt(0.5 * ((-0.10) ** 2 + (0.14 * -2.0) ** 2))
    assert result.residues[0].csp_A == pytest.approx(expected)
    assert result.cutoff == pytest.approx(0.15)
    assert result.residues[0].significant is True


def test_large_abs_delta_n_excludes_csp():
    """Residues with |ΔN_raw| > 15 ppm get no CSP; smaller |ΔN| still recorded."""
    results = compute_csp_A(
        mapping=[(1, 1), (2, 2)],
        apo_seq="MG",
        holo_seq="MG",
        H_apo={1: 8.40, 2: 8.20},
        N_apo={1: 121.0, 2: 110.0},
        H_holo={1: 8.30, 2: 8.25},
        N_holo={1: 101.0, 2: 115.0},  # |ΔN|=20 and 5
        threshold_config={"absolute_cutoff": 0.05},
        enable_referencing=False,
    )
    assert results[0].excluded_large_dn is True
    assert results[0].csp_A is None
    assert results[0].significant is None
    assert results[1].excluded_large_dn is False
    assert results[1].csp_A is not None
    assert results[1].significant is True

    excl = collect_large_dn_exclusions(results)
    assert len(excl) == 1
    assert excl[0]["apo_resi"] == 1
    assert excl[0]["abs_dN_raw"] == pytest.approx(20.0)


def test_aa_mismatch_excludes_csp():
    """Misaligned AA pairs (e.g. G↔E) get no CSP; matched pairs still recorded."""
    results = compute_csp_A(
        mapping=[(1, 1), (2, 2)],
        apo_seq="GL",
        holo_seq="EL",
        H_apo={1: 8.24, 2: 7.67},
        N_apo={1: 108.3, 2: 120.4},
        H_holo={1: 8.42, 2: 8.09},
        N_holo={1: 119.2, 2: 121.7},
        threshold_config={"absolute_cutoff": 0.05},
        enable_referencing=False,
    )
    # G↔E mismatch: no CSP, not in significance value list
    assert results[0].excluded_aa_mismatch is True
    assert results[0].csp_A is None
    assert results[0].significant is None
    # L↔L match: CSP computed
    assert results[1].excluded_aa_mismatch is False
    assert results[1].csp_A is not None
    assert results[1].significant is True
    # Only the matched CSP contributes to thresholding
    assert sum(1 for r in results if r.csp_A is not None) == 1


def test_ca_inclusive_formula():
    result = compute_ca_inclusive_csp(
        aligned_apo="MG",
        aligned_holo="MG",
        H_apo={1: 8.40},
        N_apo={1: 121.0},
        CA_apo={1: 55.0},
        H_holo={1: 8.30},
        N_holo={1: 119.0},
        CA_holo={1: 56.0},
        threshold_config={"absolute_cutoff": 0.2},
        enable_referencing=False,
    )
    expected = math.sqrt((1.0 / 3.0) * ((-0.10) ** 2 + (0.14 * -2.0) ** 2 + (0.3 * 1.0) ** 2))
    assert result.residues[0].csp_A == pytest.approx(expected)


def test_ha_ca_formula():
    result = compute_ha_ca_csp(
        aligned_apo="MG",
        aligned_holo="MG",
        HA_apo={1: 4.20},
        CA_apo={1: 55.0},
        HA_holo={1: 4.35},
        CA_holo={1: 55.7},
        threshold_config={"absolute_cutoff": 0.2},
        enable_referencing=False,
    )
    expected = math.sqrt(0.5 * ((1.0 * 0.15) ** 2 + (0.3 * 0.7) ** 2))
    assert result.residues[0].csp_A == pytest.approx(expected)


def test_grid_search_picks_expected_best_offsets():
    bundle_result = compute_nh_csp(
        mapping=[(1, 1), (2, 2)],
        apo_seq="MG",
        holo_seq="MG",
        H_apo={1: 8.40, 2: 8.20},
        N_apo={1: 121.0, 2: 110.0},
        H_holo={1: 8.50, 2: 8.30},
        N_holo={1: 122.0, 2: 111.0},
        enable_referencing=False,
        threshold_config={"absolute_cutoff": 0.0},
    )
    grid = run_offset_grid_search(
        bundle_result.residues,
        h_min=-0.2,
        h_max=0.0,
        h_step=0.1,
        n_min=-1.5,
        n_max=0.0,
        n_step=0.5,
        cutoff=0.05,
    )
    assert grid["best_h_offset"] == pytest.approx(-0.1)
    assert grid["best_n_offset"] == pytest.approx(-0.5)


def test_atom_specific_1d_csp_all_atom_types():
    base = compute_nh_csp(
        mapping=[(1, 1), (2, 2)],
        apo_seq="MG",
        holo_seq="MG",
        H_apo={1: 8.40, 2: 8.10},
        N_apo={1: 121.0, 2: 110.0},
        H_holo={1: 8.30, 2: 8.20},
        N_holo={1: 119.7, 2: 110.4},
        threshold_config={"absolute_cutoff": 0.05},
        enable_referencing=False,
    )
    for residue, ca_apo, ca_holo, ha_apo, ha_holo in zip(
        base.residues,
        [55.0, 44.0],
        [55.4, 44.3],
        [4.20, 3.95],
        [4.35, 4.10],
    ):
        residue.CA_apo = ca_apo
        residue.CA_holo_original = ca_holo
        residue.HA_apo = ha_apo
        residue.HA_holo_original = ha_holo

    h_result = compute_atom_deltas_with_offset_bundle(base.residues, "H", {"absolute_cutoff": 0.05}, {"cutoff": 10.0})
    n_result = compute_atom_deltas_with_offset_bundle(base.residues, "N", {"absolute_cutoff": 0.2}, {"cutoff": 10.0})
    ca_result = compute_atom_deltas_with_offset_bundle(base.residues, "CA", {"absolute_cutoff": 0.2}, {"cutoff": 10.0})
    ha_result = compute_atom_deltas_with_offset_bundle(base.residues, "HA", {"absolute_cutoff": 0.1}, {"cutoff": 10.0})

    assert h_result.values_abs == pytest.approx([0.1, 0.1])
    assert n_result.values_abs == pytest.approx([1.3, 0.4])
    assert ca_result.values_abs == pytest.approx([0.4, 0.3])
    assert ha_result.values_abs == pytest.approx([0.15, 0.15])


def test_bundle_matches_pipeline_nh_and_1d_behavior():
    mapping = [(1, 1), (2, 2), (3, 3)]
    apo_seq = "MGA"
    holo_seq = "MGA"
    H_apo = {1: 8.40, 2: 8.10, 3: 8.20}
    N_apo = {1: 121.0, 2: 110.0, 3: 118.0}
    H_holo = {1: 8.50, 2: 8.18, 3: 8.10}
    N_holo = {1: 122.2, 2: 109.8, 3: 116.9}

    bundle = compute_nh_csp(
        mapping,
        apo_seq,
        holo_seq,
        H_apo,
        N_apo,
        H_holo,
        N_holo,
        enable_referencing=True,
        grid_params={"h_min": -0.2, "h_max": 0.2, "h_step": 0.1, "n_min": -1.5, "n_max": 1.5, "n_step": 0.1, "cutoff": 0.2},
    )
    pipeline = compute_csp_A(
        mapping,
        apo_seq,
        holo_seq,
        H_apo,
        N_apo,
        H_holo,
        N_holo,
        enable_referencing=True,
        grid_params={"h_min": -0.2, "h_max": 0.2, "h_step": 0.1, "n_min": -1.5, "n_max": 1.5, "n_step": 0.1, "cutoff": 0.2},
    )

    assert bundle.cutoff == pytest.approx(
        compute_threshold_with_outlier_removal_pipeline(
            [row.csp_A for row in pipeline if row.csp_A is not None],
            3.0,
            0.0,
            10,
            0.2,
        ).threshold
    )
    for bundle_row, pipeline_row in zip(bundle.residues, pipeline):
        assert bundle_row.csp_A == pytest.approx(pipeline_row.csp_A)
        assert bundle_row.H_offset == pytest.approx(pipeline_row.H_offset)
        assert bundle_row.N_offset == pytest.approx(pipeline_row.N_offset)
        assert bundle_row.significant == pipeline_row.significant

    bundle_atom = compute_atom_deltas_with_offset_bundle(bundle.residues, "N", grid_params={"cutoff": 0.2})
    pipeline_atom = compute_atom_deltas_with_offset_pipeline(pipeline, "N", grid_params={"cutoff": 0.2})
    assert bundle_atom.offset == pytest.approx(pipeline_atom.offset)
    assert bundle_atom.cutoff == pytest.approx(pipeline_atom.cutoff)
    assert bundle_atom.values_abs == pytest.approx(pipeline_atom.values_abs)


def test_bundle_matches_pipeline_ca_and_ha_ca_behavior():
    ca_bundle = compute_ca_inclusive_csp(
        aligned_apo="MGA",
        aligned_holo="MGA",
        H_apo={1: 8.40, 2: 8.10, 3: 8.20},
        N_apo={1: 121.0, 2: 110.0, 3: 118.0},
        CA_apo={1: 55.0, 2: 44.0, 3: 52.0},
        H_holo={1: 8.50, 2: 8.15, 3: 8.15},
        N_holo={1: 122.2, 2: 109.7, 3: 117.2},
        CA_holo={1: 55.8, 2: 44.3, 3: 52.6},
        enable_referencing=False,
    )
    ca_pipeline = compute_csp_from_aligned_sequences_ca(
        "MGA",
        "MGA",
        {1: 8.40, 2: 8.10, 3: 8.20},
        {1: 121.0, 2: 110.0, 3: 118.0},
        {1: 55.0, 2: 44.0, 3: 52.0},
        {1: 8.50, 2: 8.15, 3: 8.15},
        {1: 122.2, 2: 109.7, 3: 117.2},
        {1: 55.8, 2: 44.3, 3: 52.6},
        enable_referencing=False,
    )
    for bundle_row, pipeline_row in zip(ca_bundle.residues, ca_pipeline):
        assert bundle_row.csp_A == pytest.approx(pipeline_row.csp_A)
        assert bundle_row.significant == pipeline_row.significant

    ha_bundle = compute_ha_ca_csp(
        aligned_apo="MGA",
        aligned_holo="MGA",
        HA_apo={1: 4.20, 2: 3.90, 3: 4.10},
        CA_apo={1: 55.0, 2: 44.0, 3: 52.0},
        HA_holo={1: 4.35, 2: 4.00, 3: 4.25},
        CA_holo={1: 55.8, 2: 44.3, 3: 52.6},
        enable_referencing=False,
    )
    ha_pipeline = compute_csp_from_aligned_sequences_ha_ca(
        "MGA",
        "MGA",
        {1: 4.20, 2: 3.90, 3: 4.10},
        {1: 55.0, 2: 44.0, 3: 52.0},
        {1: 4.35, 2: 4.00, 3: 4.25},
        {1: 55.8, 2: 44.3, 3: 52.6},
        enable_referencing=False,
    )
    for bundle_row, pipeline_row in zip(ha_bundle.residues, ha_pipeline):
        assert bundle_row.csp_A == pytest.approx(pipeline_row.csp_A)
        assert bundle_row.significant == pipeline_row.significant


def test_bundle_1d_results_align_with_generated_1d_analysis_csv(tmp_path: Path):
    try:
        from scripts.analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    except Exception as exc:
        pytest.skip(f"Optional 1D CSV cross-check unavailable in this environment: {exc}")

    csp_path = tmp_path / "csp_table.csv"
    with open(csp_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "apo_bmrb", "holo_bmrb", "holo_pdb", "chain", "apo_resi", "apo_aa", "holo_resi", "holo_aa",
            "H_apo", "N_apo", "CA_apo", "H_holo", "N_holo", "CA_holo",
            "H_offset", "N_offset", "CA_offset", "dH", "dN", "csp_A", "significant",
        ])
        writer.writerow(["18251", "4700", "1cf4", "A", "1", "M", "1", "M", "8.41", "122.0", "55.6", "8.26", "115.7", "55.8", "-0.05", "1.1", "0", "0.04", "-1.5", "0.31", "1"])
        writer.writerow(["18251", "4700", "1cf4", "A", "2", "G", "2", "G", "8.10", "110.0", "45.0", "8.12", "109.8", "45.3", "0", "0", "0", "0.02", "-0.2", "0.03", "0"])

    compute_1d_metrics_for_target(tmp_path)
    nh_result = compute_nh_csp(
        mapping=[(1, 1), (2, 2)],
        apo_seq="MG",
        holo_seq="MG",
        H_apo={1: 8.41, 2: 8.10},
        N_apo={1: 122.0, 2: 110.0},
        H_holo={1: 8.26, 2: 8.12},
        N_holo={1: 115.7, 2: 109.8},
        enable_referencing=False,
    )
    for residue, ca_apo, ca_holo in zip(nh_result.residues, [55.6, 45.0], [55.8, 45.3]):
        residue.CA_apo = ca_apo
        residue.CA_holo_original = ca_holo

    ca_result = compute_atom_deltas_with_offset_bundle(
        nh_result.residues,
        "CA",
        {"outlier_z_score": 3.0, "significance_z_score": 0.0, "max_outlier_iterations": 10, "max_outlier_fraction": 0.2},
        {"cutoff": 10.0},
    )

    with open(tmp_path / "1d_analysis.csv", "r", newline="") as handle:
        rows = list(csv.DictReader(handle))

    csv_ca_values = [float(row["CSP_CA_1d"]) for row in rows]
    csv_ca_significance = [row["csp_CA_1d_significant"] == "True" for row in rows]
    assert ca_result.values_abs == pytest.approx(csv_ca_values)
    assert [entry["significant"] for entry in ca_result.per_residue] == csv_ca_significance
