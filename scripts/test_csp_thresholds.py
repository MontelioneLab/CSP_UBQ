"""Scientific invariants for the independent, fixed-eligibility CSP sweep."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.analyze_csp_thresholds import (
    STRUCTURAL_FLAGS, analyze, at_thresholds, eligibility_reason,
    roc_table, sha256, structural_label,
)
from scripts.create_si_fig_s27 import log_density, validate_analysis
from scripts.export_target_roc_curves import export_target_rocs


def row(**updates):
    result = {k: "False" for k in STRUCTURAL_FLAGS}
    result.update(csp_A="0.1", csp_z="0", significant="1", holo_aa="A",
                  passes_sub_2A_filter_any_atom="False", min_any_atom_distance="12")
    result.update(updates)
    return result


def test_structural_union_and_unknowns():
    assert structural_label(row()) is False
    for key in (*STRUCTURAL_FLAGS, "passes_sub_2A_filter_any_atom"):
        assert structural_label(row(**{key: "True"})) is True
    assert structural_label(row(passes_sub_2A_filter_any_atom="", min_any_atom_distance="1.999")) is True
    assert structural_label(row(passes_sub_2A_filter_any_atom="", min_any_atom_distance="2")) is False
    assert structural_label({}) is None
    assert structural_label(row(has_hbond_interaction="")) is None
    assert structural_label(row(has_hbond_interaction="", is_occluded_occlusion="True")) is True


@pytest.mark.parametrize("updates,expected", [
    ({}, "eligible"), ({"csp_A": ""}, "missing_or_nonfinite_csp"),
    ({"csp_A": "inf"}, "missing_or_nonfinite_csp"), ({"csp_A": "nan"}, "missing_or_nonfinite_csp"),
    ({"csp_A": "-1"}, "negative_csp"), ({"holo_aa": "P"}, "proline"),
    ({"csp_z": "100.01"}, "extreme_csp_z"), ({"csp_z": "100"}, "eligible"),
    ({"significant": ""}, "missing_baseline_mask"),
    ({"has_hbond_interaction": ""}, "unknown_structural_label"),
])
def test_fixed_eligibility(updates, expected):
    assert eligibility_reason(row(**updates)) == expected


@pytest.mark.parametrize("labels,scores,expected", [
    ([1, 0], [2, 1], 1), ([1, 0], [1, 2], 0),
    ([1, 0, 1, 0], [2, 2, 2, 2], 0.5),
    ([1, 0, 1, 0], [3, 2, 2, 1], 0.875),
])
def test_auc_endpoints_ties_and_rank_identity(labels, scores, expected):
    curve, auc = roc_table(labels, scores)
    assert auc == pytest.approx(expected)
    labels, scores = np.asarray(labels, dtype=bool), np.asarray(scores)
    u = mannwhitneyu(scores[labels], scores[~labels]).statistic
    assert auc == pytest.approx(u / (labels.sum() * (~labels).sum()))
    assert curve.iloc[0][["tp", "fp"]].tolist() == [0, 0]
    assert curve.iloc[-1][["tpr", "fpr"]].tolist() == [1, 1]
    assert curve.tp.is_monotonic_increasing and curve.fp.is_monotonic_increasing
    assert (curve.tp + curve.fn == labels.sum()).all()
    assert (curve.fp + curve.tn == (~labels).sum()).all()
    assert np.isnan(curve.iloc[0].precision) and np.isnan(curve.iloc[0].tp_fp_ratio)
    requested = np.r_[np.inf, np.unique(scores), 0]
    looked_up = at_thresholds(curve, requested)
    for i, threshold in enumerate(requested):
        prediction = scores >= threshold
        assert looked_up.iloc[i].tp == (prediction & labels).sum()
        assert looked_up.iloc[i].fp == (prediction & ~labels).sum()


def test_single_class_and_invalid_scores():
    curve, auc = roc_table([1, 1], [1, 2])
    assert np.isnan(auc) and curve.fpr.isna().all()
    with pytest.raises(ValueError):
        roc_table([1, 0], [1, np.inf])


def test_bounded_density_and_sparse_groups():
    from scipy.integrate import trapezoid
    assert log_density([0.1, 0.2]) is None
    assert log_density([0.1] * 6) is None
    x, density = log_density([0.05, 0.08, 0.1, 0.2, 0.4])
    assert x[[0, -1]] == pytest.approx(np.log10([0.05, 0.4]))
    assert trapezoid(density, x) == pytest.approx(1)


def make_pair_fixture(tmp_path):
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    targets = tmp_path / "targets.csv"
    with targets.open("w") as f:
        f.write("holo_pdb,apo_bmrb,holo_bmrb\n1ABC,1,2\n2DEF,3,4\n")
    for name, apo, holo, values in [("1ABC_1", "1", "2", [0.2, 0.05]), ("2DEF_3", "3", "4", [0.04, 0.3])]:
        path = outputs / name
        path.mkdir()
        records = []
        for i, value in enumerate(values):
            records.append(row(csp_A=str(value), significant=str(int(value > 0.1)),
                               passes_filter_distance=str(i == 0), sequential_position=str(i + 1),
                               chain="A", apo_bmrb=apo, holo_bmrb=holo))
        with (path / "master_alignment.csv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(records)
    return targets, outputs


def test_full_analysis_bootstrap_and_immutability(tmp_path):
    targets, outputs = make_pair_fixture(tmp_path)
    run = tmp_path / "2026-11-04_test_001"
    inputs = {p: sha256(p) for p in outputs.glob("*/master_alignment.csv")}
    analyze(targets, outputs, run, bootstrap=40, seed=123)
    assert all(sha256(p) == digest for p, digest in inputs.items())
    metrics = pd.read_csv(run / "aggregate_metrics.csv")
    assert metrics.loc[metrics.aggregation == "equal_target_mean", "auc"].iloc[0] == pytest.approx(0.5)
    assert metrics.loc[metrics.score_mode == "absolute_ppm", "auc"].iloc[0] == pytest.approx(0.25)
    boot = pd.read_csv(run / "bootstrap_auc.csv")
    weights = np.random.default_rng(123).multinomial(2, [0.5, 0.5], size=40)
    expected = []
    for w in weights:
        y, scores = [], []
        for count, values in zip(w, [[0.2, 0.05], [0.04, 0.3]]):
            y.extend([1, 0] * count)
            scores.extend(values * count)
        expected.append(roc_table(y, scores)[1])
    assert boot.pooled_auc_absolute_ppm.to_numpy() == pytest.approx(expected)
    manifest = validate_analysis(run, check_upstream=True)
    assert manifest["n_targets"] == 2 and manifest["n_eligible_residues"] == 4
    for entry in manifest["sources"]:
        assert sha256(run / entry["snapshot_path"]) == entry["sha256"]
    with pytest.raises(FileExistsError):
        analyze(targets, outputs, run, bootstrap=40)
    with (run / "target_summary.csv").open("a") as stream:
        stream.write("changed\n")
    with pytest.raises(ValueError, match="changed"):
        validate_analysis(run)


def test_missing_or_mismatched_pairs_fail_before_creating_run(tmp_path):
    targets, outputs = make_pair_fixture(tmp_path)
    targets.write_text("holo_pdb,apo_bmrb,holo_bmrb\n1ABC,1,999\n")
    run = tmp_path / "2026-11-04_invalid_001"
    with pytest.raises(ValueError, match="pair mismatch"):
        analyze(targets, outputs, run, bootstrap=2)
    assert not run.exists()


def test_target_roc_exports_route_correctly_and_reuse_without_writes(tmp_path):
    targets, outputs = make_pair_fixture(tmp_path)
    run = tmp_path / "2026-11-04_test_001"
    analyze(targets, outputs, run, bootstrap=8, seed=123)
    export_id = "2026-11-04_csp-roc_001"
    index_dir = export_target_rocs(run, outputs, export_id=export_id, dpi=72)
    index = pd.read_csv(index_dir / "index.csv")
    assert len(index) == 2
    original_curves = pd.read_csv(run / "roc_curves.csv")
    for entry in index.to_dict("records"):
        directory = outputs / entry["target_id"] / export_id
        receipt = json.loads((directory / "provenance.json").read_text())
        assert receipt["signature"]["target_id"] == entry["target_id"]
        assert receipt["summary"]["auc"] == pytest.approx(entry["auc"])
        exported = pd.read_csv(directory / "roc_curve.csv")
        expected = original_curves.loc[original_curves.target_id == entry["target_id"]].reset_index(drop=True)
        pd.testing.assert_frame_equal(exported, expected, check_dtype=False, rtol=1e-11, atol=1e-12)
        from pypdf import PdfReader
        assert len(PdfReader(directory / "roc_curve.pdf").pages) == 1
        assert (directory / "roc_curve.png").is_file()
    before = {p: (sha256(p), p.stat().st_mtime_ns) for p in outputs.rglob("*") if p.is_file()}
    assert export_target_rocs(run, outputs, export_id=export_id, dpi=72) == index_dir
    assert all((sha256(p), p.stat().st_mtime_ns) == state for p, state in before.items())
    relocated = outputs.rename(tmp_path / "relocated")
    validate_analysis(run, check_upstream=True, targets_csv=targets, outputs_dir=relocated)
    assert export_target_rocs(run, relocated, export_id=export_id, dpi=72) == relocated / export_id
    outputs = relocated
    changed = outputs / "2DEF_3/master_alignment.csv"
    changed.write_text(changed.read_text() + "\n")
    with pytest.raises(ValueError, match="Target input changed"):
        export_target_rocs(run, outputs, export_id="2026-11-04_csp-roc_002", dpi=72)
    assert not list(outputs.glob("*/2026-11-04_csp-roc_002"))
    assert not (outputs / "2026-11-04_csp-roc_002").exists()


def test_target_roc_export_preserves_damaged_existing_artifacts(tmp_path):
    targets, outputs = make_pair_fixture(tmp_path)
    run = tmp_path / "2026-11-04_test_001"
    analyze(targets, outputs, run, bootstrap=8, seed=123)
    export_id = "2026-11-04_csp-roc_001"
    export_target_rocs(run, outputs, export_id=export_id, dpi=72)
    damaged = outputs / "1ABC_1" / export_id / "roc_curve.png"
    damaged.write_bytes(b"intentionally changed")
    with pytest.raises(ValueError, match="ROC artifact changed"):
        export_target_rocs(run, outputs, export_id=export_id, dpi=72)
    assert damaged.read_bytes() == b"intentionally changed"
